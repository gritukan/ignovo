mod graph;
mod peptide_search;
mod deconvolution;
mod novoviz;
mod align;

use std::fs;
use std::path::PathBuf;

use clap::{arg, Parser};

use mzdata::io::MZFileReader;
use mzdata::MZReader;
use mzdata::spectrum::MultiLayerSpectrum;

use serde::de::DeserializeOwned;
use serde::Deserialize;

use crate::deconvolution::{deconvolute, AmbiguousPeak, DeconvolutionConfig};
use crate::peptide_search::PeptideSearchConfig;
use crate::graph::{build_graph, GraphBuildingConfig};

#[derive(Debug, Deserialize)]
pub struct InputConfig {
    mzml_path: String,
    scan_index: usize,
}

fn get_spectrum(config: &InputConfig) -> Result<MultiLayerSpectrum, Box<dyn std::error::Error>> {
    let reader = MZReader::open_path(&config.mzml_path)?;
    let mut spectrum: Option<MultiLayerSpectrum> = None;
    for (i, s) in reader.enumerate() {
        if i == config.scan_index {
            spectrum = Some(s);
            break;
        }
    }

    if spectrum.is_none() {
        return Err("Spectrum not found".into());
    }

    Ok(spectrum.unwrap())
}

#[derive(Debug, Deserialize)]
pub struct GraphVisualizationConfig {
    output_html_path: String,
}

fn visualize_graph(graph: &graph::Graph, deconvoluted_spectrum: &[AmbiguousPeak], config: &GraphVisualizationConfig) {
    let mut viz_peaks = vec![];
    let mut viz_edges = vec![];
    for node in &graph.nodes {
        let mut label_lines = vec![
            format!("Mass: {:.4}", node.mass),
            format!("Intensity: {:.2}", node.intensity),
            format!("Original m/z: {:.4}", deconvoluted_spectrum[node.peak_index].original_peak.mz),
            format!("Score: {:.2}", node.score),
        ];
        if node.modification_set.modifications.len() > 0 {
            let mut mods = "Modifications:".to_string();
            for modification in &node.modification_set.modifications {
                mods += &format!(" {}", modification.name);
            }
            label_lines.push(mods);
        }

        let peak = novoviz::Peak {
            mass: node.mass,
            intensity: node.intensity,
            label: label_lines.join("\n"),
        };
        viz_peaks.push(peak);
    }
    for edge in &graph.edges {
        let mut fragment_seqs = vec![];
        for (seq, _, _) in &edge.sequences {
            fragment_seqs.push(seq.to_string());
        }

        let edge = novoviz::Edge {
            from_mass: graph.nodes[edge.from_index].mass,
            to_mass: graph.nodes[edge.to_index].mass,
            value: edge.score,
            label: fragment_seqs.join("|"),
            _idx: 0,
        };
        viz_edges.push(edge);
    }

    let html = novoviz::graph_html(viz_peaks, viz_edges);
    std::fs::write(&config.output_html_path, html).unwrap();
}

#[derive(Debug, Deserialize)]
pub struct PeptidePrintingConfig { }

fn print_peptides(peptide_paths: &[peptide_search::Path], _config: &PeptidePrintingConfig) {
    for path in peptide_paths {
        let mut seq = String::new();
        for edge in &path.edges {
            if seq.len() > 0 {
                seq += "|";
            }
            seq += edge.sequences[0].0.to_string().as_str();
        }

        println!("{}", seq);
    }
}

#[derive(Debug, Deserialize)]
pub struct AlignmentsPrintingConfig { }

fn print_alignments(alignments: &align::AlignmentResults, _: &AlignmentsPrintingConfig) {
    println!("V Gene Alignments:");
    for (allele, result) in &alignments.v_genes {
        println!("  {}: Score = {}", allele.fancy_name(), result.score);
    }

    println!("\nJ Gene Alignments:");
    for (allele, result) in &alignments.j_genes {
        println!("  {}: Score = {}", allele.fancy_name(), result.score);
    }
}

#[derive(Debug, Deserialize)]
pub struct AppConfig {
    input: InputConfig,

    #[serde(default = "default_deconvolution_config")]
    deconvolution: DeconvolutionConfig,

    #[serde(default = "default_graph_building_config")]
    build_graph: GraphBuildingConfig,

    visualize_graph: Option<GraphVisualizationConfig>,

    peptide_search: Option<PeptideSearchConfig>,

    print_peptides: Option<PeptidePrintingConfig>,

    align_peptides: Option<align::AlignmentConfig>,

    print_alignments: Option<AlignmentsPrintingConfig>,
}

fn default_deconvolution_config() -> DeconvolutionConfig {
    DeconvolutionConfig::None {
        cfg: deconvolution::NoneDeconvolutionConfig::default(),
    }
}

fn default_graph_building_config() -> GraphBuildingConfig {
    GraphBuildingConfig::default()
}

#[derive(Parser, Debug)]
struct Cli {
    /// Path to YAML config
    #[arg(long)]
    config: Option<PathBuf>,
}

pub fn from_yaml<T: DeserializeOwned>(s: &str) -> Result<T, String> {
    let de = serde_yaml::Deserializer::from_str(s); // note: not &mut
    match serde_path_to_error::deserialize::<_, T>(de) {
        Ok(v) => Ok(v),
        Err(err) => {
            let p = err.path().to_string();
            let loc = err.inner().location()
                .map(|l| format!(" at line {}, column {}", l.line(), l.column()))
                .unwrap_or_default();
            Err(format!("Failed to parse at `{}`{}: {}", p, loc, err))
        }
    }
}

fn load_config() -> AppConfig {
    let cli = Cli::parse();

    let config_path = cli.config.unwrap_or(PathBuf::from("config.yaml"));
    let config_str = fs::read_to_string(config_path).expect("Failed to read config file");
    let config = from_yaml::<AppConfig>(&config_str).expect("Failed to parse config file");

    config
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = load_config();

    let spectrum = get_spectrum(&config.input)?;
    let deconvoluted_spectrum = deconvolute(&spectrum, &config.deconvolution);

    let graph = build_graph(&deconvoluted_spectrum, &config.build_graph);

    if let Some(viz_config) = &config.visualize_graph {
        visualize_graph(&graph, &deconvoluted_spectrum, viz_config);
    }

    let mut peptide_paths = None;
    if let Some(peptide_search_config) = &config.peptide_search {
        peptide_paths = Some(peptide_search::peptide_search(&graph, peptide_search_config));
    }

    if let Some(printing_config) = &config.print_peptides {
        if let Some(peptide_paths) = peptide_paths.as_ref() {
            print_peptides(peptide_paths, printing_config);
        } else {
            eprintln!("No peptide search performed, cannot print peptides");
        }
    }

    let mut alignment_results = None;
    if let Some(align_config) = &config.align_peptides {
        if let Some(peptide_paths) = peptide_paths.as_ref() {
            alignment_results = Some(align::align(peptide_paths, align_config));
        } else {
            eprintln!("No peptide search performed, cannot align peptides");
        }
    }

    if let Some(printing_config) = &config.print_alignments {
        if let Some(alignments) = alignment_results.as_ref() {
            print_alignments(alignments, printing_config);
        } else {
            eprintln!("No alignments performed, cannot print alignments");
        }
    }

    Ok(())
}
