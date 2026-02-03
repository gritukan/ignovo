use mzdata::{mzpeaks::peak::MZPoint, spectrum::MultiLayerSpectrum};

use mzdata::prelude::SpectrumLike;
use serde::Deserialize;
use std::collections::HashSet;

/// Represents isotopic distribution for a given molecule.
#[derive(Clone, Debug)]
struct IsotopicDistribution {
    /// Minimum number of neutrons in the distribution.
    /// Small abundance peaks before this are ignored.
    min_neutrons: i32,
    /// Abundances of isotopic peaks.
    /// i-th element corresponds to (min_neutrons + i) neutrons.
    abundances: Vec<f64>,
    /// Index of the most abundant peak in the distribution.
    /// That is, the most abundant peak corresponds to
    /// (min_neutrons + most_abundant_index) neutrons.
    most_abundant_index: usize,
}

impl IsotopicDistribution {
    /// When dealing with large molecules, given the peaks corresponding to
    /// isotopic distribution, it is hard to distinguish between +-1 neutron errors
    /// trying to fit the distribution. We use pretty simple heuristic: first of all,
    /// we consider that the largest peak is the most abundant one. Then, we extend
    /// the range of considered neutrons on both sides as long as the abundance
    /// is within relative_tolerance of the most abundant peak that makes algorithm
    /// more robust to intensity measurement errors.
    /// This function returns range of additional neutrons for the most abundant peak
    /// to consider given the relative tolerance of the intensity measurements.
    fn get_neutrons_range(&self, relative_tolerance: f64) -> (i32, i32) {
        let mut min_index = self.most_abundant_index;
        while min_index > 0 && self.abundances[min_index] >= self.abundances[self.most_abundant_index] * relative_tolerance {
            min_index -= 1;
        }
        let mut max_index = self.most_abundant_index;
        while max_index < (self.abundances.len() - 1) && self.abundances[max_index + 1] >= self.abundances[self.most_abundant_index] * relative_tolerance {
            max_index += 1;
        }

        (min_index as i32 + self.min_neutrons, max_index as i32 + self.min_neutrons)
    }
}

trait IsotopicModel {
    fn get_isotopic_distribution(&self, mass: f64) -> IsotopicDistribution;
}

struct AveragineModel {
    // (mass, distribution) pairs for different lengths of averagine units
    distributions: Vec<(f64, IsotopicDistribution)>,
}

impl AveragineModel {
    pub fn new() -> Self {
        // python3 generate_averagine.py --min_aas 1 --max_aas 400 --max_neutrons 60 --output averagine.txt
        const AVERAGINE_DATA: &str = include_str!("data/averagine.txt");
    
        #[derive(Debug, Deserialize)]
        struct AveragineEntry {
            mass: f64,
            distribution: Vec<f64>,
        }
        let entries: Vec<AveragineEntry> = serde_json::from_str(AVERAGINE_DATA).unwrap();

        let mut distributions = Vec::with_capacity(entries.len());

        for entry in &entries {
            let mut max_index = 0;
            for (i, &abundance) in entry.distribution.iter().enumerate() {
                if abundance > entry.distribution[max_index] {
                    max_index = i;
                }
            }

            let mut min_neutrons = max_index;
            let mut max_neutrons = max_index;
            let mut explained_abundance = entry.distribution[max_index];

            loop {
                // Corner case: for very small molecules monoisotopic peak can explain
                // almost all abundance, but we still want to rely on additional peaks.
                if explained_abundance >= 0.95 && max_neutrons != min_neutrons {
                    break;
                }

                // Select the side with larger abundance
                let left = if min_neutrons == 0 {
                    false
                } else if max_neutrons == entry.distribution.len() - 1 {
                    true
                } else {
                    entry.distribution[min_neutrons - 1] > entry.distribution[max_neutrons + 1]
                };

                let next_abundance = if left {
                    entry.distribution[min_neutrons - 1]
                } else {
                    entry.distribution[max_neutrons + 1]
                };

                // Peak is too small to consider, so we stop here
                if next_abundance < entry.distribution[max_index] * 0.1 {
                    break;
                }

                explained_abundance += next_abundance;
                if left {
                    min_neutrons -= 1;
                } else {
                    max_neutrons += 1;
                }
            }

            let isotopic_distribution = IsotopicDistribution {
                min_neutrons: min_neutrons as i32,
                abundances: entry.distribution[min_neutrons..=max_neutrons].to_vec(),
                most_abundant_index: max_index - min_neutrons,
            };

            distributions.push((entry.mass, isotopic_distribution));
        }

        Self { distributions }
    }
}

impl IsotopicModel for AveragineModel {
    fn get_isotopic_distribution(&self, mass: f64) -> IsotopicDistribution {
        let idx = self.distributions.binary_search_by(|entry| {
            entry.0.partial_cmp(&mass).unwrap_or(std::cmp::Ordering::Equal)
        });

        let closest_idx = match idx {
            Ok(i) => i,
            Err(i) => {
                if i == 0 {
                    0
                } else if i == self.distributions.len() {
                    self.distributions.len() - 1
                } else {
                    let dist_prev = (self.distributions[i - 1].0 - mass).abs();
                    let dist_next = (self.distributions[i].0 - mass).abs();
                    if dist_prev < dist_next {
                        i - 1
                    } else {
                        i
                    }
                }
            }
        };

        self.distributions[closest_idx].1.clone()
    }
}

/// Represents a deconvoluted peak with definite mass and intensity.
#[derive(Clone, Debug)]
pub struct DeconvolutedPeak {
    /// Neutral mass of the deconvoluted peak.
    pub mass: f64,
    /// Tolerance of the deconvoluted peak.
    /// The real mass is in the range [mass - mass_tolerance, mass + mass_tolerance].
    pub mass_tolerance: f64,
    /// Intensity of the deconvoluted peak.
    pub intensity: f64,
    /// Score represents confidence of the deconvolution of the peak.
    pub score: f64,
    /// Charge of the original peak
    pub charge: i32,
    /// How many additional neutrons does original peak have
    /// compared to the monoisotopic mass.
    pub neutrons: i32,
}

/// Represents a deconvoluted peak with possible ambiguities.
#[derive(Clone, Debug)]
pub struct AmbiguousPeak {
    /// List of possible deconvoluted peak variants.
    /// Variants are sorted by score in descending order,
    /// so the first variant is the most likely one.
    pub variants: Vec<DeconvolutedPeak>,
    /// The peak from the original spectrum that this deconvoluted peak corresponds to.
    pub original_peak: MZPoint,
    /// All peaks (including isotopic peaks) from the original spectrum that
    /// this deconvoluted peak corresponds to.
    pub original_peaks: Vec<MZPoint>,
}

/// Tolerance that is used during deconvolution.
#[derive(Debug, Deserialize, Clone)]
#[serde(tag = "type", content = "value", rename_all = "lowercase")]
pub enum DeconvolutionTolerance {
    /// Parts-per-million tolerance.
    PPM(f64),
    /// Absolute Dalton tolerance.
    Da(f64),
    /// Orbitrap-specific tolerance. Decreases as sqrt(mz).
    Orbitrap { calibration_mz: f64, calibration_resolution: f64 },
}

impl DeconvolutionTolerance {
    pub fn resolution_at(&self, mz: f64) -> f64 {
        match self {
            DeconvolutionTolerance::PPM(ppm) => 1e6 / ppm,
            DeconvolutionTolerance::Da(da) => mz / da,
            DeconvolutionTolerance::Orbitrap { calibration_mz, calibration_resolution } => {
                calibration_resolution * (mz / calibration_mz).sqrt()
            }
        }
    }

    pub fn tolerance_at(&self, mz: f64) -> f64 {
        mz / self.resolution_at(mz)
    }
}

/// Model for intensity scoring during deconvolution.
#[derive(Debug, Deserialize, Clone)]
#[serde(tag = "type", content = "weight", rename_all = "lowercase")]
pub enum IntensityScoringModel {
    /// Absolute intensity scoring.
    /// Largest peak gets score of "weight", other peaks are scaled accordingly.
    Absolute(f64),
    /// Logarithmic intensity scoring.
    /// Largest peak gets score of "weight", other peaks are scaled logarithmically.
    Log(f64),
}

impl IntensityScoringModel {
    pub fn score(&self, intensity: f64, max_intensity: f64) -> f64 {
        match self {
            IntensityScoringModel::Absolute(weight) => {
                if max_intensity == 0.0 {
                    0.0
                } else {
                    intensity / max_intensity * weight
                }
            }
            IntensityScoringModel::Log(weight) => {
                if intensity <= 0.0 || max_intensity <= 0.0 {
                    0.0
                } else {
                    (intensity.ln() / max_intensity.ln()) * weight
                }
            }
        }
    }
}

/// Configuration for none deconvolution.
#[derive(Debug, Deserialize, Clone)]
pub struct NoneDeconvolutionConfig {
    /// If set, returns at most max_peaks with highest intensity.
    #[serde(default)]
    max_peaks: Option<usize>,

    /// If set, returns only peaks with intensity above the threshold.
    #[serde(default)]
    intensity_threshold: Option<f64>,

    /// If spectrum is not centroided, pick peaks with given SNR.
    #[serde(default = "default_pick_peaks_snr")]
    pick_peaks_snr: f64,

    /// Intensity scoring model to use during deconvolution.
    #[serde(default = "default_intensity_scoring_model")]
    intensity_scoring: IntensityScoringModel,

    /// Tolerance that is used during deconvolution.
    #[serde(default = "default_deconvolution_tolerance")]
    tolerance: DeconvolutionTolerance,
}

impl NoneDeconvolutionConfig {
    pub fn default() -> Self {
        Self {
            max_peaks: None,
            intensity_threshold: None,
            pick_peaks_snr: default_pick_peaks_snr(),
            intensity_scoring: default_intensity_scoring_model(),
            tolerance: default_deconvolution_tolerance(),
        }
    }
}

fn default_deconvolution_tolerance() -> DeconvolutionTolerance {
    DeconvolutionTolerance::PPM(50.0)
}

fn default_pick_peaks_snr() -> f64 {
    3.0
}

fn default_intensity_scoring_model() -> IntensityScoringModel {
    IntensityScoringModel::Log(10.0)
}

pub fn none_deconvolute(spectrum: &MultiLayerSpectrum, config: &NoneDeconvolutionConfig) -> Vec<AmbiguousPeak> {
    // If spectrum is not centroided, centroid it first.
    let mut spectrum = spectrum.clone();
    if spectrum.signal_continuity() != mzdata::spectrum::SignalContinuity::Centroid {
        spectrum.pick_peaks(config.pick_peaks_snr as f32).unwrap();
    }

    let mut peaks = spectrum.peaks().iter().collect::<Vec<_>>();

    // Apply intensity threshold if set.
    if let Some(threshold) = config.intensity_threshold {
        peaks.retain(|peak| peak.intensity as f64 >= threshold);
    }

    // Sort peaks by intensity in descending order.
    peaks.sort_unstable_by(|a, b| b.intensity.partial_cmp(&a.intensity).unwrap_or(std::cmp::Ordering::Equal));

    // Retain only top max_peaks if set.
    if let Some(max_peaks) = config.max_peaks {
        peaks.truncate(max_peaks);
    }

    let max_intensity = peaks.iter().map(|peak| peak.intensity as f64).fold(0.0, f64::max);

    peaks.into_iter().map(|peak| AmbiguousPeak {
        variants: vec![ DeconvolutedPeak {
            mass: peak.mz,
            mass_tolerance: config.tolerance.tolerance_at(peak.mz),
            intensity: peak.intensity as f64,
            score: config.intensity_scoring.score(peak.intensity as f64, max_intensity),
            charge: 1,
            neutrons: 0,
        }, ],
        original_peak: MZPoint {
            mz: peak.mz,
            intensity: peak.intensity,
        },
        original_peaks: vec![ MZPoint {
            mz: peak.mz,
            intensity: peak.intensity,
        }, ],
    }).collect::<Vec<_>>()
}

/// Configuration for greedy deconvolution.
#[derive(Debug, Deserialize, Clone)]
pub struct GreedyDeconvolutionConfig {
    /// If set, returns at most max_peaks with highest intensity.
    #[serde(default)]
    max_peaks: Option<usize>,

    /// If set, returns not more than max_variants per peak.
    #[serde(default)]
    max_variants: Option<usize>,

    /// If set, returns only peaks with intensity above the threshold.
    #[serde(default)]
    intensity_threshold: Option<f64>,

    /// Charge range to consider during deconvolution.
    #[serde(default = "default_charge_range")]
    charge_range: (i32, i32),

    /// If spectrum is not centroided, pick peaks with given SNR.
    #[serde(default = "default_pick_peaks_snr")]
    pick_peaks_snr: f64,

    /// Intensity scoring model to use during deconvolution.
    #[serde(default = "default_intensity_scoring_model")]
    intensity_scoring: IntensityScoringModel,

    /// Penalty for isotopic error during deconvolution.
    #[serde(default = "default_isotope_error_penalty")]
    isotope_error_penalty: f64,

    /// See get_neutrons_range for details.
    #[serde(default = "default_isotope_tolerance")]
    isotope_tolerance: f64,

    /// This value is added to the score for each ion
    /// with the same mass and another charge.
    #[serde(default = "default_charge_ladder_score")]
    charge_ladder_score: f64,

    /// Tolerance that is used during deconvolution.
    #[serde(default = "default_deconvolution_tolerance")]
    tolerance: DeconvolutionTolerance,
}

fn default_charge_range() -> (i32, i32) {
    (5, 15)
}

fn default_isotope_error_penalty() -> f64 {
    0.1
}

fn default_charge_ladder_score() -> f64 {
    0.1
}

fn default_isotope_tolerance() -> f64 {
    0.8
}

const NEUTRON_MASS: f64 = 1.0033548378;
const PROTON_MASS: f64 = 1.007276466812;

fn greedy_deconvolute(spectrum: &MultiLayerSpectrum, config: &GreedyDeconvolutionConfig) -> Vec<AmbiguousPeak> {
    // If spectrum is not centroided, centroid it first.
    let mut spectrum = spectrum.clone();
    spectrum.pick_peaks(config.pick_peaks_snr as f32).unwrap();

    // Main data structure over original spectra: all peaks sorted by mz.
    let mut peaks = spectrum.peaks().iter().collect::<Vec<_>>();
    peaks.sort_unstable_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));

    let get_resolution = |mz: f64| -> f64 {
        config.tolerance.resolution_at(mz)
    };
    let get_bounds = |mz: f64| -> (f64, f64) {
        let resolution = get_resolution(mz);
        let delta_mz = mz / resolution;
        (mz - delta_mz, mz + delta_mz)
    };

    // Given a mz, finds the index of the peak with maximum intensity
    // within the tolerance bounds around the mz.
    let get_max_peak_index = |mz: f64| -> Option<usize> {
        let (min_mz, max_mz) = get_bounds(mz);
        let idx = peaks.binary_search_by(
            |p| p.mz.partial_cmp(&min_mz)
                .unwrap());
        let mut idx = match idx {
            Ok(idx) => idx,
            Err(idx) => idx,
        };

        let mut result: Option<usize> = None;

        while idx < peaks.len() && peaks[idx].mz <= max_mz {
            if result.is_none() || peaks[idx].intensity > peaks[result.unwrap()].intensity {
                result = Some(idx);
            }
            idx += 1;
        }

        if result.is_some() && peaks[result.unwrap()].intensity == 0.0 {
            return None;
        }

        result
    };

    // A candidate for deconvoluted peak.
    #[derive(Clone, Debug)]
    struct Candidate {
        // The index of the peak in the original spectrum.
        // As every candidate is associated with a number of peaks
        // due to isotopic distribution, this index corresponds
        // to the most abundant peak in the distribution.
        peak_index: usize,
        // For two peaks mz1 and mz2 within the same isotopic distribution,
        // charge can be easilty determined as z = round(NEUTRON_MASS / (mz2 - mz1)).
        // However, note that to reliably distinguish between z and z + 1
        // we need a resolution of order 1 / z^2 which is not possible for
        // large z that occurs in antibodies. That's why we store a range of
        // possible charges for the candidate that are indistinguishable.
        charge_range: (i32, i32),
        // Which peaks in the original spectrum match the isotopic distribution.
        // Note, that for different charges the profile is the same due to
        // indistinguishability described above.
        profile: Vec<Option<usize>>,
        // Predicted intensity of the candidate.
        intensity: f64,
    }

    let averagine_model = AveragineModel::new();

    // The very first pass: identify potential (mz, z) candidates
    // that fit isotopic patterns well.
    let mut candidates = vec![];

    for peak in &peaks {
        if peak.intensity == 0.0 {
            continue;
        }

        let mut charges = vec![];

        for charge in config.charge_range.0..=config.charge_range.1 {
            let bounds = get_bounds(peak.mz);
            let tolerance = bounds.1 - bounds.0;
            if tolerance > NEUTRON_MASS / charge as f64 {
                dbg!(peak.mz, &config.tolerance, tolerance, charge);
                panic!("Charge too high for peak at mz {}", peak.mz);
            }

            let neutral_mass = peak.mz * charge as f64 - charge as f64 * PROTON_MASS;
            let isotopic_distribution = averagine_model.get_isotopic_distribution(neutral_mass);

            let mut matched_peaks = 0;
            let mut profile = vec![];
            for (i, _) in isotopic_distribution.abundances.iter().enumerate() {
                let expected_mz = (neutral_mass + (i as i32 - isotopic_distribution.most_abundant_index as i32) as f64 * NEUTRON_MASS) / charge as f64 + PROTON_MASS;
                let max_peak_index = get_max_peak_index(expected_mz);
                if max_peak_index.is_some() {
                    matched_peaks += 1;
                }
                profile.push(max_peak_index);
            }

            let max_missing = if isotopic_distribution.abundances.len() > 4 { 1 } else { 0 };
            if matched_peaks + max_missing >= isotopic_distribution.abundances.len() {
                charges.push((charge, profile.clone()));
            }
        }

        let mut idx = 0;
        while idx < charges.len() {
            let start_charge = charges[idx].0;
            let mut end_charge = start_charge;

            while idx + 1 < charges.len() && charges[idx + 1].0 == end_charge + 1 && charges[idx + 1].1 == charges[idx].1 {
                end_charge = charges[idx + 1].0;
                idx += 1;
            }

            candidates.push(Candidate {
                peak_index: peaks.iter().position(|p| p.mz == peak.mz).unwrap(),
                charge_range: (start_charge, end_charge),
                profile: charges[idx].1.clone(),
                intensity: 0.0,
            });

            idx += 1;
        }
    }

    // Gets an intensity for a candidate based on how well its profile
    // fits the expected isotopic distribution.
    let candidate_intensity = |candidate: &Candidate, peaks: &[MZPoint]| -> f64 {
        let mut intensity = f64::MAX;

        let charge = (candidate.charge_range.0 + candidate.charge_range.1) / 2;
        let neutral_mass = peaks[candidate.peak_index].mz * charge as f64 - charge as f64 * PROTON_MASS;
        let isotopic_distribution = averagine_model.get_isotopic_distribution(neutral_mass);

        for (i, &peak_idx_opt) in candidate.profile.iter().enumerate() {
            if let Some(peak_idx) = peak_idx_opt {
                intensity = intensity.min(
                    peaks[peak_idx].intensity as f64 / isotopic_distribution.abundances[i]
                );
            }
        }

        intensity
    };

    // For each peak of original spectrum, the list of candidates that
    // have it in the profile.
    let mut peak_to_candidates: Vec<Vec<usize>> = vec![vec![]; peaks.len()];
    for (i, candidate) in candidates.iter().enumerate() {
        for &peak_idx_opt in &candidate.profile {
            if let Some(peak_idx) = peak_idx_opt {
                peak_to_candidates[peak_idx].push(i);
            }
        }
    }

    // Compute intensities for all candidates.
    for candidate in &mut candidates {
        candidate.intensity = candidate_intensity(candidate, &peaks);
    }

    // Now let's greedily evaluate candidate intensities.
    // At each step we extract the current best candidate,
    // call its intesity final, and then reduce the intensities
    // of all other candidates accordingly.
    let mut final_candidates = vec![];
    let mut bst: std::collections::BTreeSet<(ordered_float::OrderedFloat<f64>, usize)> = std::collections::BTreeSet::new();
    for (i, candidate) in candidates.iter().enumerate() {
        bst.insert((ordered_float::OrderedFloat(candidate.intensity), i));
    }

    // Whether or not a candidate was already processed and final intensity assigned.
    let mut candidate_processed = vec![false; candidates.len()];

    while let Some((ordered_float::OrderedFloat(intensity), candidate_idx)) = bst.iter().next_back().cloned() {
        bst.remove(&(ordered_float::OrderedFloat(intensity), candidate_idx));

        let candidate = &candidates[candidate_idx];
        candidate_processed[candidate_idx] = true;

        if intensity < 1e-3 {
            continue;
        }

        final_candidates.push(candidate.clone());

        // Now decrease intensity of peaks in the original spectrum.
        // Also, remember which candidates are affected.
        let charge = (candidate.charge_range.0 + candidate.charge_range.1) / 2;
        let neutral_mass = peaks[candidate.peak_index].mz * charge as f64 - charge as f64 * PROTON_MASS;
        let isotopic_distribution = averagine_model.get_isotopic_distribution(neutral_mass);

        let mut affected_candidates = vec![];

        for (i, &peak_idx_opt) in candidate.profile.iter().enumerate() {
            if let Some(peak_idx) = peak_idx_opt {
                peaks[peak_idx].intensity -= (intensity * isotopic_distribution.abundances[i]) as f32;
                affected_candidates.extend(peak_to_candidates[peak_idx].clone());
            }
        }

        affected_candidates.sort_unstable();
        affected_candidates.dedup();
        affected_candidates.retain(|&candidate| !candidate_processed[candidate]);

        for candidate in affected_candidates {
            // remove old BST entry (if present), recompute intensity and re-insert
            let old_key = (ordered_float::OrderedFloat(candidates[candidate].intensity), candidate);
            bst.remove(&old_key);

            let new_intensity = candidate_intensity(&candidates[candidate], &peaks);
            candidates[candidate].intensity = new_intensity;

            bst.insert((ordered_float::OrderedFloat(new_intensity), candidate));
        }
    }

    // Now convert final candidates to deconvoluted peaks.
    // For each candidate, consider all possible charges
    // as well as possible isotopic errors in case of large
    // masses where it is hard to distinguish beween +-1 neutron
    // isotopic distributions.
    let mut final_peaks = vec![];

    for candidate in final_candidates {
        let mut candidate_peaks = vec![];

        let peak = &peaks[candidate.peak_index];
        for charge in candidate.charge_range.0..=candidate.charge_range.1 {
            let neutral_mass = peak.mz * charge as f64 - charge as f64 * PROTON_MASS;
            let isotopic_distribution = averagine_model.get_isotopic_distribution(neutral_mass);

            let (min_neutrons, max_neutrons) = isotopic_distribution.get_neutrons_range(0.8);
            let abundant_neutrons = isotopic_distribution.most_abundant_index as i32 + isotopic_distribution.min_neutrons;

            for neutrons in min_neutrons..=max_neutrons {
                let adjusted_mass = neutral_mass - neutrons as f64 * NEUTRON_MASS;
                let (min_mz, max_mz) = get_bounds(peak.mz);
                let mass_tolerance = (max_mz - min_mz) / 2.0 * charge as f64;
                let neutrons_error = (neutrons - abundant_neutrons).abs() as f64;
                candidate_peaks.push(DeconvolutedPeak {
                    mass: adjusted_mass,
                    mass_tolerance,
                    intensity: candidate.intensity,
                    score: -neutrons_error * config.isotope_error_penalty,
                    charge,
                    neutrons,
                });
            }
        }

        let original_peaks = candidate.profile.iter()
            .filter_map(|&idx_opt| idx_opt.map(|idx| MZPoint {
                mz: peaks[idx].mz,
                intensity: peaks[idx].intensity,
            }))
            .collect::<Vec<_>>();

        final_peaks.push(AmbiguousPeak {
            variants: candidate_peaks,
            original_peak: peak.clone(),
            original_peaks,
        });
    }

    // Add intensity scores.
    let max_intensity = final_peaks.iter()
        .flat_map(|ap| ap.variants.iter())
        .map(|dp| dp.intensity)
        .fold(0.0, f64::max);
    for peak in &mut final_peaks {
        for variant in &mut peak.variants {
            variant.score += config.intensity_scoring.score(variant.intensity, max_intensity);
        }
    }

    // Add charge-ladder bonus.
    {
        let mut variants: Vec<(f64, f64, i32)> = Vec::new();
        for peak in &final_peaks {
            for variant in &peak.variants {
                // (m + z*PROTON) / z = mz => z = m / (mz - PROTON)
                let charge = ((peak.original_peak.mz - PROTON_MASS) / variant.mass).round() as i32;
                variants.push((variant.mass, variant.mass_tolerance, charge));
            }
        }

        variants.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let max_other_tol = variants.iter().map(|(_, t, _)| *t).fold(0.0_f64, f64::max);

        for peak in &mut final_peaks {
            for variant in &mut peak.variants {
                let mut charges: HashSet<i32> = HashSet::new();

                let search_radius = variant.mass_tolerance + max_other_tol;
                let lower = variant.mass - search_radius;
                let upper = variant.mass + search_radius;

                let start_idx = match variants.binary_search_by(|(m, _, _)| {
                    m.partial_cmp(&lower).unwrap_or(std::cmp::Ordering::Equal)
                }) {
                    Ok(i) => i,
                    Err(i) => i,
                };

                let self_charge = ((peak.original_peak.mz - PROTON_MASS) / variant.mass).round() as i32;

                let mut idx = start_idx;
                while idx < variants.len() && variants[idx].0 <= upper {
                    let (variant_mass, variant_tolerance, variant_charge) = variants[idx];
                    if variant_charge != self_charge && (variant_mass - variant.mass).abs() <= (variant.mass_tolerance + variant_tolerance) {
                        charges.insert(variant_charge);
                    }
                    idx += 1;
                }

                let count = charges.len() as f64;
                variant.score += config.charge_ladder_score * count;
            }
        }
    }

    // For each peak, sort variants by score and retain only top max_variants if set.
    for peak in &mut final_peaks {
        peak.variants.sort_unstable_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
        if let Some(max_variants) = config.max_variants {
            peak.variants.truncate(max_variants);
        }
    }

    // Sort final peaks by score and retain only top max_peaks if set.
    final_peaks.sort_unstable_by(|a, b| b.variants[0].score.partial_cmp(&a.variants[0].score).unwrap_or(std::cmp::Ordering::Equal));
    if let Some(max_peaks) = config.max_peaks {
        final_peaks.truncate(max_peaks);
    }

    dbg!(final_peaks.len());

    final_peaks
}

/// Configuration for deconvolution.
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "method", rename_all = "lowercase")]
pub enum DeconvolutionConfig {
    None {
        #[serde(flatten)]
        cfg: NoneDeconvolutionConfig,
    },
    Greedy {
        #[serde(flatten)]
        cfg: GreedyDeconvolutionConfig,
    },
}

pub fn deconvolute(spectrum: &MultiLayerSpectrum, config: &DeconvolutionConfig) -> Vec<AmbiguousPeak> {
    match config {
        DeconvolutionConfig::None { cfg } => none_deconvolute(spectrum, cfg),
        DeconvolutionConfig::Greedy { cfg } => greedy_deconvolute(spectrum, cfg),
    }
}
