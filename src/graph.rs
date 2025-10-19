use std::cmp::Ordering;

use mzdeisotope::DeconvolvedSolutionPeak;
use rustyms::{molecular_formula, prelude::{CheckedAminoAcid, Chemical, MolecularFormula, Peptidoform, SequenceElement}, quantities::Tolerance, sequence::UnAmbiguous, system::{dalton, Mass, OrderedMass}};
use serde::Deserialize;
use mzdeisotope::solution::Envelope;

use crate::deconvolution::{AmbiguousPeak, DeconvolutedPeak};

/// Describes a modification that may occur to a peak.
/// Examples include loss of water or ammonia.
#[derive(Debug, Clone, Deserialize)]
pub struct Modification {
    /// Name of the modification.
    pub name: String,
    /// Mass change due to the modification.
    /// Can be negative (water loss) or positive (proton addition).
    pub mass_change: f64,
    /// How the modification affects node score. Typically negative.
    pub score: f64,
}

/// Represents a set of modifications applied to a peak.
#[derive(Debug, Clone)]
pub struct ModificationSet {
    /// List of modifications in the set.
    pub modifications: Vec<Modification>,
    /// Total mass change due to the modifications.
    pub total_mass: Mass,
    /// Total score change due to the modifications.
    pub total_score: f64,
}

/// Represents a node in the graph. Each node is a unambiguous
/// while many nodes can correspond to the same ambiguous peak
/// due to, for example, charge errors or modifications.
#[derive(Debug, Clone)]
pub struct Node {
    /// Neutral mass of the node.
    pub mass: f64,
    /// Tolerance of the mass.
    /// The real mass is in the range [mass - mass_tolerance, mass + mass_tolerance].
    pub mass_tolerance: f64,
    /// Intensity of the (deconvoluted) peak.
    pub intensity: f64,
    /// Set of modifications applied to the peak.
    pub modification_set: ModificationSet,
    /// Score represents confidence of such ion actually existing.
    /// Higher score means higher confidence.
    pub score: f64,
    /// Each node corresponds to a (possibly ambiguous) peak in the deconvoluted spectrum.
    /// This field stores original index of that peak.
    pub peak_index: usize,
}

/// An edge in the graph represents possible amino acid compositions
/// such that the mass difference between two nodes corresponds to
/// the mass of the amino acid composition (with some error).
#[derive(Debug, Clone)]
pub struct Edge {
    /// The index of the starting node of this edge.
    pub from_index: usize,
    /// The index of the ending node of this edge.
    pub to_index: usize,
    /// List of the sequences that can be read along this edge
    /// together with mass error and score.
    pub sequences: Vec<(Peptidoform<UnAmbiguous>, Mass, f64)>,
    /// Score of the edge.
    /// Higher score means higher confidence that this edge is correct.
    pub score: f64,
}

/// Represents the graph built from deconvoluted peaks.
#[derive(Debug, Clone)]
pub struct Graph {
    /// Nodes of the graph.
    pub nodes: Vec<Node>,

    /// Edges of the graph.
    pub edges: Vec<Edge>,

    /// Outgoing edges for each node.
    pub out_edges: Vec<Vec<Edge>>,
    /// Incoming edges for each node.
    pub in_edges: Vec<Vec<Edge>>,
}

/// Various options for building the graph.
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "type", content = "weight", rename_all = "lowercase")]

pub enum MassErrorModel {
    /// Mass error score is -weight * |mass_error|
    Linear(f64),
    /// Mass error score is -weight * ln(max(|mass_error|, 0.01))
    Log(f64),
}

impl MassErrorModel {
    pub fn score(&self, mass_error: Mass) -> f64 {
        match self {
            MassErrorModel::Linear(weight) => -weight * mass_error.value.abs(),
            MassErrorModel::Log(weight) => {
                let error = mass_error.value.abs().max(0.01);
                -weight * error.ln()
            }
        }
    }
}

/// Various options for building the graph.
#[derive(Debug, Clone, Deserialize)]
pub struct GraphBuildingConfig {
    /// Maximum length of fragments (in amino acids) to consider when building the graph.
    #[serde(default = "default_max_fragment_length")]
    pub max_fragment_length: usize,

    /// List of peak modifications to consider when building the graph.
    #[serde(default)]
    pub modifications: Vec<Modification>,

    /// Maximum number of modifications to apply per peak.
    #[serde(default = "default_max_modifications_per_peak")]
    pub max_modifications_per_peak: usize,

    /// Mass error model to use when scoring edges.
    #[serde(default = "default_mass_error_model")]
    pub mass_error_model: MassErrorModel,

    /// Length penalty for edges. Score of an edge is reduced by length_penalty * (length - 1).
    #[serde(default = "default_length_penalty")]
    pub length_penalty: f64,
}

impl GraphBuildingConfig {
    pub fn default() -> Self {
        GraphBuildingConfig {
            max_fragment_length: default_max_fragment_length(),
            modifications: vec![],
            max_modifications_per_peak: default_max_modifications_per_peak(),
            mass_error_model: default_mass_error_model(),
            length_penalty: default_length_penalty(),
        }
    }
}

fn default_mass_error_model() -> MassErrorModel {
    MassErrorModel::Linear(5.0)
}

fn default_max_fragment_length() -> usize {
    1
}

fn default_max_modifications_per_peak() -> usize {
    1
}

fn default_length_penalty() -> f64 {
    5.0
}

fn generate_fragments(max_length: usize) -> Vec<(Peptidoform<UnAmbiguous>, Mass)> {
    let amino_acids = vec![
        CheckedAminoAcid::Alanine,
        CheckedAminoAcid::Arginine,
        CheckedAminoAcid::Asparagine,
        CheckedAminoAcid::AsparticAcid,
        CheckedAminoAcid::Cysteine,
        CheckedAminoAcid::GlutamicAcid,
        CheckedAminoAcid::Glutamine,
        CheckedAminoAcid::Glycine,
        CheckedAminoAcid::Histidine,
        CheckedAminoAcid::AmbiguousLeucine,
        CheckedAminoAcid::Lysine,
        CheckedAminoAcid::Methionine,
        CheckedAminoAcid::Phenylalanine,
        CheckedAminoAcid::Proline,
        CheckedAminoAcid::Serine,
        CheckedAminoAcid::Threonine,
        CheckedAminoAcid::Tryptophan,
        CheckedAminoAcid::Tyrosine,
        CheckedAminoAcid::Valine,
    ];
    // There are 20 standard amino acids, but Leucine and Isoleucine are isobaric.
    assert!(amino_acids.len() == 19);

    let mut fragments: Vec<(Peptidoform<UnAmbiguous>, Mass)> = vec![];
    let mut prefixes: Vec<Vec<SequenceElement<UnAmbiguous>>> = vec![vec![]];
    for _ in 0..max_length {
        let mut current_fragments = vec![];
        for prefix in &prefixes {
            for aa in &amino_acids {
                if !prefix.is_empty() && prefix.last().unwrap().aminoacid < *aa {
                    continue;
                }

                let mut new_fragment = prefix.clone();
                let element = SequenceElement::new(*aa, None);
                new_fragment.push(element);
                current_fragments.push(new_fragment);
            }
        }
        for fragment in &current_fragments {
            let mut mass = Mass::new::<dalton>(0.0);
            for element in fragment {
                mass += element.aminoacid.formula().monoisotopic_mass();
            }
            fragments.push((Peptidoform::new(fragment.clone()), mass));
        }
        prefixes = current_fragments
    }

    fragments
}

fn generate_modification_sets(config: &GraphBuildingConfig) -> Vec<ModificationSet> {
    let mut sets = vec![ModificationSet {
        modifications: vec![],
        total_mass: Mass::new::<dalton>(0.0),
        total_score: 0.0,
    }];

    for _ in 1..=config.max_modifications_per_peak {
        let mut new_sets = vec![];
        for set in &sets {
            if set.modifications.len() + 1 > config.max_modifications_per_peak {
                continue;
            }
            for modification in &config.modifications {
                if !set.modifications.is_empty() && set.modifications.last().unwrap().name > modification.name {
                    continue;
                }

                let mut new_mods = set.modifications.clone();
                new_mods.push(modification.clone());
                let total_mass = set.total_mass + Mass::new::<dalton>(modification.mass_change);
                let total_score = set.total_score + modification.score;
                new_sets.push(ModificationSet { modifications: new_mods, total_mass, total_score });
            }
        }

        sets.extend(new_sets);
    }

    sets
}

pub fn build_graph(peaks: &Vec<AmbiguousPeak>, config: &GraphBuildingConfig) -> Graph {
    let modification_sets = generate_modification_sets(config);

    let mut nodes = vec![];

    // Generate nodes from peaks.
    for (peak_index, peak) in peaks.iter().enumerate() {
        // For each peak consider all possible variants.
        for variant in &peak.variants {
            // For each variant consider all possible modification sets.
            for modification_set in &modification_sets {
                // "Demodified" mass that (hopefully) corresponds to the mass of b/y ion.
                // Note that modification mass is subtracted here. Suppose the modification
                // is -18 Da (loss of water). Then the mass should be increased compared
                // to observed peak.
                let mass = variant.mass - modification_set.total_mass.value;
                if mass > 0.0 {
                    nodes.push(Node {
                        mass: mass,
                        mass_tolerance: variant.mass_tolerance,
                        intensity: variant.intensity,
                        modification_set: modification_set.clone(),
                        score: variant.score + modification_set.total_score,
                        peak_index: peak_index,
                    });
                }
            }
        }
    }

    // Some nodes may have the same or very close masses (for example, +1 proton vs +1 neutron).
    // We keep only the node with the highest score among such nodes.
    // TODO: sometimes this is an evidence of having ions with different charges
    // but the same mass which should be encouraged with additional score for the node.
    nodes.sort_unstable_by(|a, b| a.mass.partial_cmp(&b.mass).unwrap_or(Ordering::Equal));

    let mut filtered_nodes: Vec<Node> = vec![];
    for node in nodes {
        if !filtered_nodes.is_empty() {
            let node_min_mass = node.mass - node.mass_tolerance;
            let last_node = filtered_nodes.last().unwrap();
            let last_node_max_mass = last_node.mass + last_node.mass_tolerance;
            if node_min_mass > last_node_max_mass {
                filtered_nodes.push(node);
            } else {
                // Merge nodes by keeping the one with the highest score.
                if node.score > last_node.score {
                    filtered_nodes.pop();
                    filtered_nodes.push(node);
                }
            }
        } else {
            filtered_nodes.push(node);
        }
    }

    filtered_nodes.sort_unstable_by(|a, b| a.mass.partial_cmp(&b.mass).unwrap_or(Ordering::Equal));

    let mut fragments = generate_fragments(config.max_fragment_length);
    fragments.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut graph = Graph {
        nodes: filtered_nodes.clone(),
        edges: Vec::new(),
        out_edges: vec![Vec::new(); filtered_nodes.len()],
        in_edges: vec![Vec::new(); filtered_nodes.len()],
    };

    // Now generate edges.
    // Iterate over all pairs of nodes and their variants.
    for i in 0..graph.nodes.len() {
        for j in i + 1..graph.nodes.len() {
            // Since each mass is known approximately, we can find the range of possible mass differences.
            let tolerance = graph.nodes[i].mass_tolerance + graph.nodes[j].mass_tolerance;
            let mass = graph.nodes[j].mass - graph.nodes[i].mass;
            let min_mass = Mass::new::<dalton>(mass - tolerance);
            let max_mass = Mass::new::<dalton>(mass + tolerance);

            // Perform binary search to find the first fragment that has mass >= min_mass.
            let mut sequences = vec![];
            let start_index = fragments.binary_search_by(|f| f.1.partial_cmp(&min_mass).unwrap_or(Ordering::Equal));
            let mut start_index = match start_index {
                Ok(idx) => idx,
                Err(idx) => idx,
            };

            // Iterate over fragments until mass > max_mass.
            while start_index < fragments.len() && fragments[start_index].1 <= max_mass {
                // The fragment fits the mass difference, so it can be read along this edge.
                let fragment = &fragments[start_index];
                let mass_error = Mass::new::<dalton>( (graph.nodes[j].mass - graph.nodes[i].mass) - fragment.1.value ).abs();
                // The score of the fragment is based on mass error and length.
                let score = config.mass_error_model.score(mass_error) - config.length_penalty * (fragment.0.sequence().len() as f64 - 1.0);
                sequences.push((fragment.0.clone(), mass_error, score));
                start_index += 1;
            }

            // Sort sequences by score.
            sequences.sort_unstable_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(Ordering::Equal));

            if !sequences.is_empty() {
                // The score of the edge is the score of the best sequence.
                let score = sequences[0].2;
                let edge = Edge {
                    from_index: i,
                    to_index: j,
                    sequences,
                    score,
                };
                graph.out_edges[i].push(edge.clone());
                graph.in_edges[j].push(edge.clone());
                graph.edges.push(edge);
            }
        }
    }

    graph
}
