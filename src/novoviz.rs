use anyhow::{anyhow, Result};
use serde::{Deserialize, Serialize};
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Peak in the spectrum
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Peak {
    pub mass: f64,
    pub intensity: f64,
    pub label: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeconvolutedPeak {
    pub mass: f64,
    pub intensity: f64,
    pub label: String,
    /// list of original m/z that contributed to this deconvoluted peak
    pub sources_mz: Vec<f64>,
}

/// Edge between two peaks (mass -> mass), `value` is an arbitrary score/weight
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge {
    pub from_mass: f64,
    pub to_mass: f64,
    pub value: f64,
    pub label: String, // shown on-arc; can be overridden by annotation in the HTML UI
    #[serde(default)]
    pub _idx: usize, // internal for lane assignment; not serialized from JSON
}

/// Edge after lane routing
#[derive(Debug, Clone, Serialize)]
struct EdgeOut {
    from_mass: f64,
    to_mass: f64,
    value: f64,
    label: String,
    lane: usize,
}

#[derive(Debug, Clone, Serialize)]
struct GraphOut {
    peaks: Vec<Peak>,
    edges: Vec<EdgeOut>,
}

#[derive(Debug, Clone, Serialize)]
struct DeconvolutionOut {
    original: Vec<Peak>,
    deconvoluted: Vec<DeconvolutedPeak>,
}

// compile-time embed the HTML file:
const GRAPH_TEMPLATE: &str = include_str!("static/graph.html");
const DECONVOLUTION_TEMPLATE: &str = include_str!("static/deconvolution.html");

/// Generate the graph visualization HTML.
pub fn graph_html(mut peaks: Vec<Peak>, mut edges: Vec<Edge>) -> String {
    peaks.sort_by(|a, b| a.mass.partial_cmp(&b.mass).unwrap());

    // Normalize each edge so from_mass <= to_mass
    for e in edges.iter_mut() {
        if e.to_mass < e.from_mass {
            std::mem::swap(&mut e.from_mass, &mut e.to_mass);
        }
    }

    // Assign non-overlapping lanes (interval graph coloring)
    let lanes = assign_lanes(&edges);

    let edges_out: Vec<EdgeOut> = edges
        .iter()
        .enumerate()
        .map(|(i, e)| EdgeOut {
            from_mass: e.from_mass,
            to_mass: e.to_mass,
            value: e.value,
            label: e.label.clone(),
            lane: lanes[i],
        })
        .collect();

    let graph = GraphOut { peaks, edges: edges_out };

    // Build HTML by injecting JSON into template
    let tpl = GRAPH_TEMPLATE;
    let json = serde_json::to_string_pretty(&graph).unwrap();
    let html = tpl.replace("__GRAPH_JSON__", &json);
    html
}

fn assign_lanes(edges: &Vec<Edge>) -> Vec<usize> {
    #[derive(Clone)]
    struct Interval {
        start: f64,
        end: f64,
        original_idx: usize,
        span: f64,
    }
    let mut intervals: Vec<Interval> = edges
        .iter()
        .map(|e| {
            let start = e.from_mass.min(e.to_mass);
            let end = e.from_mass.max(e.to_mass);
            Interval {
                start,
                end,
                span: end - start,
                original_idx: e._idx,
            }
        })
        .collect();

    // Sort by start, then by longer span first (nesting looks nicer)
    intervals.sort_by(|a, b| {
        let c = a.start.partial_cmp(&b.start).unwrap();
        if c == std::cmp::Ordering::Equal {
            b.span.partial_cmp(&a.span).unwrap()
        } else {
            c
        }
    });

    let mut lanes: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut idx_to_lane = vec![0usize; edges.len()];

    'place: for iv in intervals {
        for (lane_idx, lane) in lanes.iter_mut().enumerate() {
            if !overlaps_any(iv.start, iv.end, lane) {
                lane.push((iv.start, iv.end));
                idx_to_lane[iv.original_idx] = lane_idx;
                continue 'place;
            }
        }
        lanes.push(vec![(iv.start, iv.end)]);
        idx_to_lane[iv.original_idx] = lanes.len() - 1;
    }
    idx_to_lane
}

fn overlaps_any(start: f64, end: f64, lane: &Vec<(f64, f64)>) -> bool {
    for (s, e) in lane {
        if !(end <= *s || start >= *e) {
            return true;
        }
    }
    false
}

pub fn deconvolution_html(
    original_peaks: Vec<Peak>,
    deconvoluted_peaks: Vec<DeconvolutedPeak>,
) -> String {
    let deconv = DeconvolutionOut {
        original: original_peaks,
        deconvoluted: deconvoluted_peaks,
    };
    let tpl = DECONVOLUTION_TEMPLATE;
    let json = serde_json::to_string_pretty(&deconv).unwrap();
    let html = tpl.replace("__DATA_JSON__", &json);
    html
}

