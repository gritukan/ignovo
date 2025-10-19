use serde::Deserialize;

use crate::graph;

/// A path in the graph representing some peptide sequence.
#[derive(Debug, Clone)]
pub struct Path {
    /// Nodes in the path.
    pub nodes: Vec<graph::Node>,
    /// Edges in the path. It's guaranteed that edges[i]
    /// connects some variants in nodes[i] and nodes[i+1].
    pub edges: Vec<graph::Edge>,
}

/// Peptide search configuration.
#[derive(Debug, Deserialize)]
pub struct PeptideSearchConfig {
    #[serde(default = "default_min_peptide_length")]
    /// Minimum length of the peptide (in number of amino acids).
    pub min_peptide_length: usize,

    #[serde(default)]
    /// Minimum score of the peptide path to search.
    pub min_score: Option<f64>,

    #[serde(default)]
    /// Maximum number of peptides to search for.
    pub max_peptides: Option<usize>,

    /// Do not look for peptides with score less than
    /// x * first found score.
    #[serde(default)]
    pub relative_score_threshold: Option<f64>,
}

fn default_min_peptide_length() -> usize {
    5
}

pub fn peptide_search(
    graph: &graph::Graph,
    config: &PeptideSearchConfig,
) -> Vec<Path> {
    let mut paths = Vec::new();

    let max_peak_index = graph
        .nodes
        .iter()
        .map(|node| node.peak_index)
        .max()
        .unwrap_or(0);

    // Many nodes can correspond to the same peak but we allow using each peak
    // only in a single peptide. This is done to avoid duplicate peptides that differ
    // only by water loss, for example. In the future this might be changed but will
    // require more sophisticated path search algorithms.
    let mut peak_used = vec![false; max_peak_index + 1];

    let mut first_found_score = None;

    loop {
        let mut best_score = vec![0.0f64; graph.nodes.len()];
        let mut next_edge = vec![None; graph.nodes.len()];
        for i in (0..graph.nodes.len()).rev() {
            let node = &graph.nodes[i];
            if peak_used[node.peak_index] {
                continue;
            }

            best_score[i] = node.score;

            for edge in &graph.out_edges[i] {
                let j = edge.to_index;
                if peak_used[graph.nodes[j].peak_index] {
                    continue;
                }

                let score = node.score + edge.score + best_score[j];
                if score > best_score[i] {
                    best_score[i] = score;
                    next_edge[i] = Some(edge.clone());
                }
            }
        }

        let mut start = None;
        for i in 0..graph.nodes.len() {
            if peak_used[graph.nodes[i].peak_index] {
                continue;
            }
            if start.is_none() || best_score[i] > best_score[start.unwrap()] {
                start = Some(i);
            }
        }

        if start.is_none() {
            break;
        }

        let path_score = best_score[start.unwrap_or(0)];
        if let Some(min_score) = config.min_score {
            if path_score < min_score {
                break;
            }
        }

        let mut start = start.unwrap();

        let mut path = Path {
            nodes: Vec::new(),
            edges: Vec::new(),
        };

        loop {
            path.nodes.push(graph.nodes[start].clone());

            peak_used[graph.nodes[start].peak_index] = true;

            let next_edge = &next_edge[start];

            if next_edge.is_none() {
                break;
            }

            path.edges.push(next_edge.as_ref().unwrap().clone());
            start = next_edge.as_ref().unwrap().to_index;
        }

        if path.nodes.len() - 1 < config.min_peptide_length {
            break;
        }

        if let Some(relative_threshold) = config.relative_score_threshold {
            if let Some(first_score) = first_found_score {
                if path_score < first_score * relative_threshold {
                    break;
                }
            } else {
                first_found_score = Some(path_score);
            }
        }

        paths.push(path);

        if let Some(max_peptides) = config.max_peptides {
            if paths.len() >= max_peptides {
                break;
            }
        }
    }

    paths
}
