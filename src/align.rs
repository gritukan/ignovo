use std::{collections::HashSet, str::FromStr};

use serde::Deserialize;

use rustyms::{align::{AlignScoring, AlignType}, imgt::{Allele, AlleleSelection, ChainType, GeneType, Selection, Species}, prelude::Peptidoform, sequence::SimpleLinear};

use crate::peptide_search::Path;

#[derive(Debug, Clone, Deserialize)]
pub struct AlignmentConfig {
    /// Species to consider during alignment.
    /// If not specified, all species will be considered.
    #[serde(default)]
    pub species: Option<Vec<String>>,

    /// Maximum number of alignments to return per gene type (V/J).
    #[serde(default = "default_max_alignments")]
    pub max_alignments: Option<usize>,
}

fn default_max_alignments() -> Option<usize> {
    Some(5)
}

pub struct AlignmentResult {
    pub score: f64,
}

pub struct AlignmentResults<'a> {
    pub v_genes: Vec<(Allele<'a>, AlignmentResult)>,
    pub j_genes: Vec<(Allele<'a>, AlignmentResult)>,
}

pub fn align_allele<'a>(
    peptides: &Vec<Path>,
    allele: &Allele<'a>,
) -> AlignmentResult {
    let mut peptidoforms = vec![];
    for peptide in peptides {
        let peptide_string = peptide
            .edges
            .iter()
            .map(|e| e.sequences[0].0.to_string())
            .collect::<Vec<String>>()
            .join("");
        peptidoforms.push(Peptidoform::pro_forma(&peptide_string, None).unwrap());

        // We do not know whether it is b or y ion, so we check both peptide
        // and its reverse.
        let rev_peptide_string = peptide_string.chars().rev().collect::<String>();
        peptidoforms.push(Peptidoform::pro_forma(&rev_peptide_string, None).unwrap());
    }

    let mut total_score = 0.0;
    for peptidoform in &peptidoforms {
        let peptidoform = peptidoform.as_simple_linear().unwrap();
        let allele = allele.sequence.as_simple_linear().unwrap();
        let alignment = rustyms::align::align::<4, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(
            peptidoform,
            allele,
            AlignScoring::default(),
            AlignType::LOCAL);
        total_score += alignment.score().absolute as f64;
    }

    AlignmentResult { score: total_score }
}

pub fn align_impl<'a>(
    peptides: &Vec<Path>,
    config: &AlignmentConfig,
    gene_type: GeneType
) -> Vec<(Allele<'a>, AlignmentResult)> {
    let species = if let Some(species_names) = &config.species {
        let mut s = HashSet::new();
        for species in species_names {
            s.insert(Species::from_str(species).unwrap());
        }
        Some(s)
    } else {
        None
    };

    let selection = Selection {
        species: species,
        chains: Option::<HashSet::<ChainType>>::None,
        genes: Some(HashSet::from([gene_type])),
        allele: AlleleSelection::First,
    };

    let mut result = vec![];

    for germline in selection.germlines() {
        let alignment = align_allele(peptides, &germline);
        result.push((germline, alignment));
    }

    result.sort_by(|a, b| b.1.score.partial_cmp(&a.1.score).unwrap());
    if let Some(max_alignments) = config.max_alignments {
        result.truncate(max_alignments);
    }

    result
}

pub fn align<'a>(
    peptides: &Vec<Path>,
    config: &AlignmentConfig,
) -> AlignmentResults<'a> {
    let v_alignments = align_impl(peptides, config, GeneType::V);
    let j_alignments = align_impl(peptides, config, GeneType::J);

    AlignmentResults {
        v_genes: v_alignments,
        j_genes: j_alignments,
    }
}
