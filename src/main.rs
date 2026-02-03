mod graph;
mod peptide_search;
mod deconvolution;
mod novoviz;
mod align;

use std::{fs, io};
use std::path::PathBuf;

use clap::{arg, Parser};

use mzdata::io::{MZFileReader, SpectrumWriter};
use mzdata::mzpeaks::{CentroidPeak, MZPeakSetType};
use mzdata::prelude::SpectrumLike;
use mzdata::{MZReader, MzMLWriter};
use mzdata::spectrum::MultiLayerSpectrum;

use serde::de::{self, DeserializeOwned};
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
pub struct DeconvolutionVisualizationConfig {
    output_html_path: String,
}

fn visualize_deconvolution(original_spectrum: &MultiLayerSpectrum, deconvoluted_spectrum: &[AmbiguousPeak], config: &DeconvolutionVisualizationConfig) {
    let mut viz_original_peaks = vec![];
    let mut viz_deconvoluted_peaks = vec![];
    for peak in original_spectrum.peaks().iter() {
        let viz_peak = novoviz::Peak {
            mass: peak.mz as f64,
            intensity: peak.intensity as f64,
            label: format!("m/z: {:.4}\nIntensity: {:.2}", peak.mz, peak.intensity),
        };
        viz_original_peaks.push(viz_peak);
    }
    for peak in deconvoluted_spectrum {
        for variant in &peak.variants {
            let viz_peak = novoviz::DeconvolutedPeak {
                mass: variant.mass,
                intensity: variant.intensity,
                label: format!("Mass: {:.4}\nIntensity: {:.2}\nCharge: {}\nNeutrons: {}", variant.mass, variant.intensity, variant.charge, variant.neutrons),
                sources_mz: peak.original_peaks.iter().map(|p| p.mz).collect(),
            };
            viz_deconvoluted_peaks.push(viz_peak);
        }
    }

    let html = novoviz::deconvolution_html(viz_original_peaks, viz_deconvoluted_peaks);
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
pub struct SaveDeconvolutedSpectrumConfig {
    output_path: String,
}

#[derive(Debug, Deserialize)]
pub struct AppConfig {
    input: InputConfig,

    #[serde(default = "default_deconvolution_config")]
    deconvolution: DeconvolutionConfig,

    visualize_deconvolution: Option<DeconvolutionVisualizationConfig>,

    save_deconvoluted_spectrum: Option<SaveDeconvolutedSpectrumConfig>,

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

    let mut spectrum = get_spectrum(&config.input)?;
    // TODO: move it somewhere else, it is here just for deconvolution visualization
    if spectrum.signal_continuity() != mzdata::spectrum::SignalContinuity::Centroid {
        spectrum.pick_peaks(3.0).unwrap();
    }

    let deconvoluted_spectrum = deconvolute(&spectrum, &config.deconvolution);

    if let Some(viz_config) = &config.visualize_deconvolution {
        visualize_deconvolution(&spectrum, &deconvoluted_spectrum, viz_config);
    }

    let graph = build_graph(&deconvoluted_spectrum, &config.build_graph);

    if let Some(viz_config) = &config.visualize_graph {
        visualize_graph(&graph, &deconvoluted_spectrum, viz_config);
    }

    if let Some(save_config) = &config.save_deconvoluted_spectrum {
        let mut writer = MzMLWriter::new(io::BufWriter::new(fs::File::create(&save_config.output_path)?));

        let spectrum_description = mzdata::spectrum::SpectrumDescription {
            id: "BlahBlahBlah".to_string(),
            index: 0,
            ms_level: 1,
            polarity: mzdata::spectrum::ScanPolarity::Positive,
            signal_continuity: mzdata::spectrum::SignalContinuity::Centroid,
            params: vec![],
            acquisition: mzdata::spectrum::Acquisition {
                scans: vec![],
                combination: mzdata::spectrum::ScanCombination::NoCombination,
                params: None,
            },
            precursor: vec![],
        };

        let peaks: Vec<CentroidPeak> = deconvoluted_spectrum
            .iter()
            .enumerate()
            .map(|(i, peak)| CentroidPeak {
                mz: peak.variants[0].mass,
                intensity: peak.variants[0].intensity as f32,
                index: i as u32,
            })
            .collect();
        let peaks = MZPeakSetType::<CentroidPeak>::new(peaks);

        let spectrum = mzdata::Spectrum {
            description: spectrum_description,
            arrays: None,
            peaks: Some(peaks),
            deconvoluted_peaks: None,
        };

        writer.write_spectrum(&spectrum)?;

        writer.close()?;
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
