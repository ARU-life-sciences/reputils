//! Evidence-based semi-automatic curation of aligned transposable-element copies.
//!
//! Coordinates in Rust values and JSON are zero-based, half-open unless named otherwise.
//! TSV boundary positions are one-based and inclusive for curator-facing interoperability.

use bio::io::fasta;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use thiserror::Error;

pub const SCHEMA_VERSION: &str = "reputils.curation.v1";

#[derive(Debug, Error)]
pub enum CurationError {
    #[error("I/O error for {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("invalid FASTA: {0}")]
    Fasta(String),
    #[error("alignment is empty")]
    EmptyAlignment,
    #[error("alignment record {id} has length {observed}, expected {expected}")]
    UnequalLength {
        id: String,
        observed: usize,
        expected: usize,
    },
    #[error("alignment record {id} contains unsupported symbol '{symbol}' at column {position}")]
    InvalidSymbol {
        id: String,
        symbol: char,
        position: usize,
    },
    #[error("duplicate alignment record ID: {0}")]
    DuplicateId(String),
    #[error("only {observed} copies were supplied; at least {minimum} are required")]
    TooFewCopies { observed: usize, minimum: usize },
    #[error("invalid parameter: {0}")]
    InvalidParameter(String),
    #[error("output directory already exists: {0}; use --force to replace curation output files")]
    OutputExists(String),
    #[error("manifest error: {0}")]
    Manifest(String),
    #[error("external tool error: {0}")]
    ExternalTool(String),
    #[error("serialization error: {0}")]
    Serialization(#[from] serde_json::Error),
}

fn io_error(path: &Path, source: std::io::Error) -> CurationError {
    CurationError::Io {
        path: path.display().to_string(),
        source,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignedCopy {
    pub id: String,
    pub sequence: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct Alignment {
    copies: Vec<AlignedCopy>,
    width: usize,
}

impl Alignment {
    pub fn from_fasta(path: &Path, minimum_copies: usize) -> Result<Self, CurationError> {
        let reader =
            fasta::Reader::from_file(path).map_err(|e| CurationError::Fasta(e.to_string()))?;
        let mut records = Vec::new();
        for item in reader.records() {
            let record = item.map_err(|e| CurationError::Fasta(e.to_string()))?;
            records.push(AlignedCopy {
                id: record.id().to_string(),
                sequence: record.seq().to_vec(),
            });
        }
        Self::from_records(records, minimum_copies)
    }

    pub fn from_records(
        mut copies: Vec<AlignedCopy>,
        minimum_copies: usize,
    ) -> Result<Self, CurationError> {
        if copies.is_empty() {
            return Err(CurationError::EmptyAlignment);
        }
        if copies.len() < minimum_copies {
            return Err(CurationError::TooFewCopies {
                observed: copies.len(),
                minimum: minimum_copies,
            });
        }
        let width = copies[0].sequence.len();
        if width == 0 {
            return Err(CurationError::EmptyAlignment);
        }
        let mut ids = HashSet::new();
        for copy in &mut copies {
            if !ids.insert(copy.id.clone()) {
                return Err(CurationError::DuplicateId(copy.id.clone()));
            }
            if copy.sequence.len() != width {
                return Err(CurationError::UnequalLength {
                    id: copy.id.clone(),
                    observed: copy.sequence.len(),
                    expected: width,
                });
            }
            for (i, base) in copy.sequence.iter_mut().enumerate() {
                *base = base.to_ascii_uppercase();
                if !is_accepted(*base) {
                    return Err(CurationError::InvalidSymbol {
                        id: copy.id.clone(),
                        symbol: char::from(*base),
                        position: i + 1,
                    });
                }
            }
        }
        Ok(Self { copies, width })
    }

    pub fn copies(&self) -> &[AlignedCopy] {
        &self.copies
    }
    pub fn width(&self) -> usize {
        self.width
    }
    pub fn len(&self) -> usize {
        self.copies.len()
    }
    pub fn is_empty(&self) -> bool {
        self.copies.is_empty()
    }
    pub fn trim(&self, left: usize, right: usize) -> Self {
        Self {
            copies: self
                .copies
                .iter()
                .map(|c| AlignedCopy {
                    id: c.id.clone(),
                    sequence: c.sequence[left..right].to_vec(),
                })
                .collect(),
            width: right - left,
        }
    }
}

fn is_accepted(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'N'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'-'
            | b'.'
            | b'?'
    )
}
fn is_called(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'U')
}
fn is_gap(b: u8) -> bool {
    matches!(b, b'-' | b'.')
}
fn is_unknown(b: u8) -> bool {
    matches!(b, b'N' | b'?')
}
fn is_ambiguous(b: u8) -> bool {
    matches!(
        b,
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V'
    )
}
fn normalized(b: u8) -> u8 {
    if b == b'U' {
        b'T'
    } else {
        b
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentColumnStats {
    pub position: usize,
    pub non_gap_count: usize,
    pub called_base_count: usize,
    pub ambiguous_count: usize,
    pub unknown_count: usize,
    pub occupancy: f64,
    pub gap_fraction: f64,
    pub n_fraction: f64,
    pub major_base: Option<char>,
    pub major_base_frequency_among_called_bases: f64,
    pub major_base_frequency_among_all_sequences: f64,
    pub entropy: f64,
    pub copy_start_count: usize,
    pub copy_end_count: usize,
    pub gap_open_count: usize,
    pub gap_close_count: usize,
    pub consensus_base: Option<char>,
    pub consensus_confidence: Option<f64>,
}

pub fn column_statistics(alignment: &Alignment) -> Vec<AlignmentColumnStats> {
    let n = alignment.len() as f64;
    let extents: Vec<_> = alignment
        .copies
        .iter()
        .map(|c| copy_extent(&c.sequence))
        .collect();
    (0..alignment.width)
        .map(|position| {
            let mut counts = [0usize; 4];
            let mut non_gap = 0;
            let mut unknown = 0;
            let mut ambiguous = 0;
            let mut starts = 0;
            let mut ends = 0;
            let mut opens = 0;
            let mut closes = 0;
            for (row, extent) in alignment.copies.iter().zip(&extents) {
                let b = normalized(row.sequence[position]);
                if !is_gap(b) {
                    non_gap += 1;
                }
                if is_unknown(b) {
                    unknown += 1;
                }
                if is_ambiguous(b) {
                    ambiguous += 1;
                }
                if let Some(i) = base_index(b) {
                    counts[i] += 1;
                }
                if extent.map(|x| x.0) == Some(position) {
                    starts += 1;
                }
                if extent.map(|x| x.1) == Some(position) {
                    ends += 1;
                }
                if is_gap(b) && position > 0 && !is_gap(row.sequence[position - 1]) {
                    opens += 1;
                }
                if !is_gap(b) && position > 0 && is_gap(row.sequence[position - 1]) {
                    closes += 1;
                }
            }
            let called: usize = counts.iter().sum();
            let (major_i, major_n) = counts
                .iter()
                .copied()
                .enumerate()
                .max_by_key(|x| x.1)
                .unwrap_or((0, 0));
            let entropy = if called == 0 {
                0.0
            } else {
                counts
                    .iter()
                    .filter(|&&x| x > 0)
                    .map(|&x| {
                        let p = x as f64 / called as f64;
                        -p * p.log2()
                    })
                    .sum()
            };
            AlignmentColumnStats {
                position,
                non_gap_count: non_gap,
                called_base_count: called,
                ambiguous_count: ambiguous,
                unknown_count: unknown,
                occupancy: non_gap as f64 / n,
                gap_fraction: (alignment.len() - non_gap) as f64 / n,
                n_fraction: unknown as f64 / n,
                major_base: if major_n == 0 {
                    None
                } else {
                    Some(['A', 'C', 'G', 'T'][major_i])
                },
                major_base_frequency_among_called_bases: if called == 0 {
                    0.0
                } else {
                    major_n as f64 / called as f64
                },
                major_base_frequency_among_all_sequences: major_n as f64 / n,
                entropy,
                copy_start_count: starts,
                copy_end_count: ends,
                gap_open_count: opens,
                gap_close_count: closes,
                consensus_base: None,
                consensus_confidence: None,
            }
        })
        .collect()
}

fn base_index(b: u8) -> Option<usize> {
    match normalized(b) {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}
fn copy_extent(seq: &[u8]) -> Option<(usize, usize)> {
    let first = seq.iter().position(|&b| !is_gap(b))?;
    let last = seq.iter().rposition(|&b| !is_gap(b))?;
    Some((first, last))
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CopyMetadata {
    pub subject_id: Option<String>,
    pub genomic_start: Option<u64>,
    pub genomic_end: Option<u64>,
    pub strand: Option<String>,
    pub locus_id: Option<String>,
    pub unextended_start: Option<u64>,
    pub unextended_end: Option<u64>,
    pub extended_start: Option<u64>,
    pub extended_end: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CopyMetrics {
    pub copy_id: String,
    pub first_non_gap_column: Option<usize>,
    pub last_non_gap_column: Option<usize>,
    pub aligned_span: usize,
    pub called_bases: usize,
    pub gap_bases: usize,
    pub unknown_bases: usize,
    pub ambiguous_bases: usize,
    pub internal_gap_runs: usize,
    pub terminal_gap_left: usize,
    pub terminal_gap_right: usize,
    pub identity_to_provisional_consensus: Option<f64>,
    pub identity_to_refined_consensus: Option<f64>,
    pub covered_fraction: f64,
    pub quality_score: f64,
    pub selected_for_consensus: bool,
    pub exclusion_reason: Option<String>,
    #[serde(flatten)]
    pub metadata: CopyMetadata,
}

pub fn copy_metrics(
    alignment: &Alignment,
    provisional: Option<&[u8]>,
    metadata: &HashMap<String, CopyMetadata>,
) -> Vec<CopyMetrics> {
    alignment
        .copies
        .iter()
        .map(|copy| {
            let extent = copy_extent(&copy.sequence);
            let left = extent.map_or(alignment.width, |x| x.0);
            let right = extent.map_or(0, |x| x.1);
            let span = extent.map_or(0, |x| x.1 - x.0 + 1);
            let called = copy.sequence.iter().filter(|&&b| is_called(b)).count();
            let gaps = copy.sequence.iter().filter(|&&b| is_gap(b)).count();
            let unknown = copy.sequence.iter().filter(|&&b| is_unknown(b)).count();
            let ambiguous = copy.sequence.iter().filter(|&&b| is_ambiguous(b)).count();
            let internal_runs = extent.map_or(0, |(a, b)| count_gap_runs(&copy.sequence[a..=b]));
            let coverage = called as f64 / alignment.width as f64;
            let internal_penalty = internal_runs as f64 / span.max(1) as f64;
            let quality = (0.65 * coverage
                + 0.35 * (1.0 - (unknown + ambiguous) as f64 / span.max(1) as f64)
                - internal_penalty)
                .clamp(0.0, 1.0);
            CopyMetrics {
                copy_id: copy.id.clone(),
                first_non_gap_column: extent.map(|x| x.0),
                last_non_gap_column: extent.map(|x| x.1),
                aligned_span: span,
                called_bases: called,
                gap_bases: gaps,
                unknown_bases: unknown,
                ambiguous_bases: ambiguous,
                internal_gap_runs: internal_runs,
                terminal_gap_left: left,
                terminal_gap_right: alignment
                    .width
                    .saturating_sub(right + usize::from(extent.is_some())),
                identity_to_provisional_consensus: provisional
                    .and_then(|p| sequence_identity(&copy.sequence, p, 1).map(|x| x.0)),
                identity_to_refined_consensus: None,
                covered_fraction: coverage,
                quality_score: quality,
                selected_for_consensus: false,
                exclusion_reason: None,
                metadata: metadata.get(&copy.id).cloned().unwrap_or_default(),
            }
        })
        .collect()
}

fn count_gap_runs(seq: &[u8]) -> usize {
    seq.iter()
        .enumerate()
        .filter(|(i, b)| is_gap(**b) && (*i == 0 || !is_gap(seq[*i - 1])))
        .count()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryConfig {
    pub smoothing_window: usize,
    pub minimum_body_length: usize,
    pub minimum_occupancy: f64,
    pub minimum_copies: usize,
    pub extension: usize,
}
impl Default for BoundaryConfig {
    fn default() -> Self {
        Self {
            smoothing_window: 11,
            minimum_body_length: 50,
            minimum_occupancy: 0.45,
            minimum_copies: 3,
            extension: 0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "kebab-case")]
pub enum ArchitectureHint {
    TirLike,
    LtrLike,
    LineLike,
    NonAutonomous,
    Asymmetric,
    Unknown,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryScoreComponents {
    pub alignment_body_support: f64,
    pub boundary_concentration: f64,
    pub occupancy_transition: f64,
    pub structural_support: f64,
    pub tsd_support: f64,
    pub domain_containment: f64,
    pub domain_cut_penalty: f64,
    pub conserved_flank_penalty: f64,
    pub insufficient_copy_penalty: f64,
    pub total: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryCandidate {
    pub left: usize,
    pub right: usize,
    pub score: f64,
    pub components: BoundaryScoreComponents,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryResult {
    pub left: usize,
    pub right: usize,
    pub left_confidence: f64,
    pub right_confidence: f64,
    pub architecture_hint: Option<ArchitectureHint>,
    pub alternatives: Vec<BoundaryCandidate>,
    pub winning_score: BoundaryScoreComponents,
    pub warnings: Vec<CurationWarning>,
}

/// Infer a supported interval. Per-column body support is
/// `0.50 occupancy + 0.35 called identity + 0.15 (1-normalized entropy)`.
/// Pair scores add mean internal support, recurrent start/end concentration, and
/// occupancy transitions. Confidence is evidence agreement, not a probability.
pub fn infer_boundaries(
    stats: &[AlignmentColumnStats],
    cfg: &BoundaryConfig,
    copies: usize,
) -> Result<BoundaryResult, CurationError> {
    if cfg.smoothing_window == 0 || cfg.minimum_body_length == 0 {
        return Err(CurationError::InvalidParameter(
            "windows and body length must be positive".into(),
        ));
    }
    if stats.is_empty() {
        return Err(CurationError::EmptyAlignment);
    }
    let raw: Vec<f64> = stats
        .iter()
        .map(|s| {
            0.50 * s.occupancy
                + 0.35 * s.major_base_frequency_among_called_bases
                + 0.15 * (1.0 - (s.entropy / 2.0).min(1.0))
        })
        .collect();
    let support = smooth(&raw, cfg.smoothing_window);
    let threshold = 0.50 * cfg.minimum_occupancy + 0.25;
    let mut intervals = sustained_intervals(&support, threshold);
    intervals.retain(|x| x.1 - x.0 >= cfg.minimum_body_length.min(stats.len()));
    if intervals.is_empty() {
        intervals.push((0, stats.len()));
    }
    let mut candidates = Vec::new();
    for (body_l, body_r) in intervals {
        let radius = cfg.smoothing_window.max(2);
        let left_range =
            body_l.saturating_sub(radius)..=(body_l + radius).min(body_r.saturating_sub(1));
        let right_start = body_r.saturating_sub(radius + 1).max(body_l);
        let right_range = right_start..=(body_r + radius).min(stats.len() - 1);
        let mut lefts: Vec<_> = left_range.collect();
        let mut rights: Vec<_> = right_range.map(|x| x + 1).collect();
        lefts.sort_by(|a, b| {
            boundary_left_score(*b, stats, &support, copies)
                .total_cmp(&boundary_left_score(*a, stats, &support, copies))
        });
        rights.sort_by(|a, b| {
            boundary_right_score(*b, stats, &support, copies)
                .total_cmp(&boundary_right_score(*a, stats, &support, copies))
        });
        for &left in lefts.iter().take(5) {
            for &right in rights.iter().take(5) {
                if right <= left || right - left < cfg.minimum_body_length.min(stats.len()) {
                    continue;
                }
                let body = support[left..right].iter().sum::<f64>() / (right - left) as f64;
                let concentration = (stats[left].copy_start_count as f64
                    + stats[right - 1].copy_end_count as f64)
                    / (2.0 * copies as f64);
                let transition = (support[left] - if left == 0 { 0.0 } else { support[left - 1] })
                    .max(0.0)
                    + (support[right - 1]
                        - if right == stats.len() {
                            0.0
                        } else {
                            support[right]
                        })
                    .max(0.0);
                let insuff = if copies < cfg.minimum_copies {
                    (cfg.minimum_copies - copies) as f64 / cfg.minimum_copies as f64
                } else {
                    0.0
                };
                let total = 0.65 * body + 0.20 * concentration + 0.15 * transition - 0.20 * insuff;
                let components = BoundaryScoreComponents {
                    alignment_body_support: body,
                    boundary_concentration: concentration,
                    occupancy_transition: transition,
                    structural_support: 0.0,
                    tsd_support: 0.0,
                    domain_containment: 0.0,
                    domain_cut_penalty: 0.0,
                    conserved_flank_penalty: 0.0,
                    insufficient_copy_penalty: insuff,
                    total,
                };
                candidates.push(BoundaryCandidate {
                    left: left.saturating_sub(cfg.extension),
                    right: (right + cfg.extension).min(stats.len()),
                    score: total,
                    components,
                });
            }
        }
    }
    candidates.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| a.left.cmp(&b.left))
            .then_with(|| a.right.cmp(&b.right))
    });
    candidates.dedup_by(|a, b| a.left == b.left && a.right == b.right);
    let best = candidates
        .first()
        .cloned()
        .ok_or_else(|| CurationError::InvalidParameter("no valid boundary pair".into()))?;
    let left_signal = boundary_left_score(best.left, stats, &support, copies).clamp(0.0, 1.0);
    let right_signal = boundary_right_score(best.right, stats, &support, copies).clamp(0.0, 1.0);
    let architecture = if left_signal > 0.55 && right_signal > 0.55 {
        ArchitectureHint::Unknown
    } else if right_signal > left_signal + 0.2 {
        ArchitectureHint::LineLike
    } else {
        ArchitectureHint::Asymmetric
    };
    let mut warnings = Vec::new();
    if left_signal < 0.45 {
        warnings.push(CurationWarning::new(
            "uncertain_left_boundary",
            "Left boundary evidence is diffuse or conflicting",
            "warning",
        ));
    }
    if right_signal < 0.45 {
        warnings.push(CurationWarning::new(
            "uncertain_right_boundary",
            "Right boundary evidence is diffuse or conflicting",
            "warning",
        ));
    }
    Ok(BoundaryResult {
        left: best.left,
        right: best.right,
        left_confidence: left_signal,
        right_confidence: right_signal,
        architecture_hint: Some(architecture),
        alternatives: candidates.into_iter().skip(1).take(9).collect(),
        winning_score: best.components,
        warnings,
    })
}
fn smooth(v: &[f64], window: usize) -> Vec<f64> {
    let half = window / 2;
    (0..v.len())
        .map(|i| {
            let a = i.saturating_sub(half);
            let b = (i + half + 1).min(v.len());
            v[a..b].iter().sum::<f64>() / (b - a) as f64
        })
        .collect()
}
fn sustained_intervals(v: &[f64], t: f64) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut start = None;
    for (i, &x) in v.iter().enumerate() {
        if x >= t && start.is_none() {
            start = Some(i)
        } else if x < t {
            if let Some(s) = start.take() {
                out.push((s, i))
            }
        }
    }
    if let Some(s) = start {
        out.push((s, v.len()))
    }
    out
}
fn boundary_left_score(i: usize, s: &[AlignmentColumnStats], sm: &[f64], n: usize) -> f64 {
    let prev = if i == 0 { 0.0 } else { sm[i - 1] };
    0.55 * (sm[i] - prev).max(0.0) + 0.45 * s[i].copy_start_count as f64 / n as f64
}
fn boundary_right_score(i: usize, s: &[AlignmentColumnStats], sm: &[f64], n: usize) -> f64 {
    let p = i.saturating_sub(1).min(s.len() - 1);
    let next = if i >= sm.len() { 0.0 } else { sm[i] };
    0.55 * (sm[p] - next).max(0.0) + 0.45 * s[p].copy_end_count as f64 / n as f64
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusPosition {
    pub alignment_position: usize,
    pub base: char,
    pub confidence: f64,
    pub called_weight: f64,
    pub gap_weight: f64,
    pub deletion_polymorphism: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefinedConsensus {
    pub sequence: String,
    pub profile: Vec<ConsensusPosition>,
}

pub fn select_representatives(
    alignment: &Alignment,
    metrics: &mut [CopyMetrics],
    min_quality: f64,
    near_identity: f64,
    min_overlap: usize,
    redundancy_control: bool,
) -> Vec<usize> {
    let mut order: Vec<usize> = (0..metrics.len()).collect();
    order.sort_by(|&a, &b| {
        metrics[b]
            .quality_score
            .total_cmp(&metrics[a].quality_score)
            .then_with(|| metrics[a].copy_id.cmp(&metrics[b].copy_id))
    });
    let mut selected: Vec<usize> = Vec::new();
    for i in order {
        if metrics[i].quality_score < min_quality {
            metrics[i].exclusion_reason = Some("low_alignment_quality".into());
            continue;
        }
        let redundant = redundancy_control
            && selected.iter().any(|&j| {
                sequence_identity(
                    &alignment.copies[i].sequence,
                    &alignment.copies[j].sequence,
                    min_overlap,
                )
                .is_some_and(|x| x.0 >= near_identity)
            });
        if redundant {
            metrics[i].exclusion_reason = Some("near_identical_redundancy".into())
        } else {
            metrics[i].selected_for_consensus = true;
            selected.push(i)
        }
    }
    if selected.is_empty() {
        if let Some((i, _)) = metrics
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.quality_score.total_cmp(&b.1.quality_score))
        {
            metrics[i].selected_for_consensus = true;
            metrics[i].exclusion_reason = None;
            selected.push(i)
        }
    }
    selected
}

pub fn refined_consensus(
    alignment: &Alignment,
    selected: &[usize],
    metrics: &[CopyMetrics],
    left: usize,
    right: usize,
) -> RefinedConsensus {
    let mut sequence = String::new();
    let mut profile = Vec::new();
    for pos in left..right {
        let mut counts = [0.0; 4];
        let mut gap = 0.0;
        for &i in selected {
            let w = metrics[i].quality_score.max(0.1);
            let b = normalized(alignment.copies[i].sequence[pos]);
            if let Some(j) = base_index(b) {
                counts[j] += w
            } else if is_gap(b) {
                let e = copy_extent(&alignment.copies[i].sequence);
                if e.is_some_and(|(a, z)| pos >= a && pos <= z) {
                    gap += w
                }
            }
        }
        let called = counts.iter().sum::<f64>();
        let total = called + gap;
        let max = counts.iter().copied().fold(0.0, f64::max);
        let deletion = total > 0.0 && gap / total >= 0.25 && gap / total <= 0.75;
        let mut present = 0u8;
        for (i, &x) in counts.iter().enumerate() {
            if called > 0.0 && x / called >= 0.20 {
                present |= 1 << i
            }
        }
        let base = if called < 1.0 || max / called < 0.40 {
            'N'
        } else {
            iupac(present)
        };
        let conf = if total == 0.0 {
            0.0
        } else {
            (max / total) * (selected.len().min(5) as f64 / 5.0)
        };
        sequence.push(base);
        profile.push(ConsensusPosition {
            alignment_position: pos,
            base,
            confidence: conf,
            called_weight: called,
            gap_weight: gap,
            deletion_polymorphism: deletion,
        })
    }
    RefinedConsensus { sequence, profile }
}
fn iupac(mask: u8) -> char {
    match mask {
        1 => 'A',
        2 => 'C',
        4 => 'G',
        8 => 'T',
        3 => 'M',
        5 => 'R',
        9 => 'W',
        6 => 'S',
        10 => 'Y',
        12 => 'K',
        7 => 'V',
        11 => 'H',
        13 => 'D',
        14 => 'B',
        15 => 'N',
        _ => 'N',
    }
}
fn sequence_identity(a: &[u8], b: &[u8], min_overlap: usize) -> Option<(f64, usize)> {
    let mut same = 0;
    let mut overlap = 0;
    for (&x, &y) in a.iter().zip(b) {
        if is_called(x) && is_called(y) {
            overlap += 1;
            if normalized(x) == normalized(y) {
                same += 1
            }
        }
    }
    if overlap < min_overlap {
        None
    } else {
        Some((same as f64 / overlap as f64, overlap))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructuralEvidence {
    pub repeat_type: String,
    pub candidate_length: usize,
    pub alignment_score: f64,
    pub identity: f64,
    pub coverage: f64,
    pub left_start: usize,
    pub left_end: usize,
    pub right_start: usize,
    pub right_end: usize,
    pub strand_relation: String,
    pub supporting_copies: usize,
}
pub fn terminal_repeat_evidence(
    seq: &str,
    min_len: usize,
    max_len: usize,
) -> Vec<StructuralEvidence> {
    let b = seq.as_bytes();
    let mut out = Vec::new();
    for len in min_len..=max_len.min(b.len() / 2) {
        let left = &b[..len];
        let right = &b[b.len() - len..];
        for (kind, target, rel) in [
            ("LTR", right.to_vec(), "direct"),
            ("TIR", revcomp(right), "inverted"),
        ] {
            let matches = left.iter().zip(&target).filter(|(x, y)| x == y).count();
            let identity = matches as f64 / len as f64;
            if identity >= 0.60 {
                out.push(StructuralEvidence {
                    repeat_type: kind.into(),
                    candidate_length: len,
                    alignment_score: matches as f64 - (len - matches) as f64,
                    identity,
                    coverage: 1.0,
                    left_start: 0,
                    left_end: len,
                    right_start: b.len() - len,
                    right_end: b.len(),
                    strand_relation: rel.into(),
                    supporting_copies: 1,
                })
            }
        }
    }
    out.sort_by(|a, b| b.alignment_score.total_cmp(&a.alignment_score));
    out.truncate(10);
    out
}
fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubfamilyEvidence {
    pub candidate_clusters: usize,
    pub cluster_membership: BTreeMap<String, usize>,
    pub mean_within_cluster_distance: Option<f64>,
    pub mean_between_cluster_distance: Option<f64>,
    pub support_for_splitting: f64,
    pub insufficient_pairs: usize,
}
#[allow(clippy::needless_range_loop)]
pub fn subfamily_analysis(alignment: &Alignment, minimum_overlap: usize) -> SubfamilyEvidence {
    let n = alignment.len();
    let mut d = vec![vec![None; n]; n];
    let mut vals = Vec::new();
    let mut insufficient = 0;
    for i in 0..n {
        for j in i + 1..n {
            if let Some((id, _)) = sequence_identity(
                &alignment.copies[i].sequence,
                &alignment.copies[j].sequence,
                minimum_overlap,
            ) {
                d[i][j] = Some(1.0 - id);
                d[j][i] = d[i][j];
                vals.push(1.0 - id)
            } else {
                insufficient += 1
            }
        }
    }
    let threshold = if vals.is_empty() {
        0.15
    } else {
        (vals.iter().sum::<f64>() / vals.len() as f64 + 0.08).min(0.25)
    };
    let mut cluster: Vec<usize> = (0..n).collect();
    fn root(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x]
        }
        x
    }
    for i in 0..n {
        for j in i + 1..n {
            if d[i][j].is_some_and(|x| x <= threshold) {
                let a = root(&mut cluster, i);
                let b = root(&mut cluster, j);
                cluster[a] = b
            }
        }
    }
    let mut ids = HashMap::new();
    let mut membership = BTreeMap::new();
    let mut next = 0;
    for (i, c) in alignment.copies.iter().enumerate() {
        let r = root(&mut cluster, i);
        let id = *ids.entry(r).or_insert_with(|| {
            let x = next;
            next += 1;
            x
        });
        membership.insert(c.id.clone(), id);
    }
    let mut within = Vec::new();
    let mut between = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            if let Some(x) = d[i][j] {
                if membership[&alignment.copies[i].id] == membership[&alignment.copies[j].id] {
                    within.push(x)
                } else {
                    between.push(x)
                }
            }
        }
    }
    let mean = |v: &[f64]| {
        if v.is_empty() {
            None
        } else {
            Some(v.iter().sum::<f64>() / v.len() as f64)
        }
    };
    let wi = mean(&within);
    let be = mean(&between);
    let support = match (wi, be) {
        (Some(a), Some(b)) => ((b - a) / b.max(0.001)).clamp(0.0, 1.0),
        _ => 0.0,
    };
    SubfamilyEvidence {
        candidate_clusters: next,
        cluster_membership: membership,
        mean_within_cluster_distance: wi,
        mean_between_cluster_distance: be,
        support_for_splitting: support,
        insufficient_pairs: insufficient,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainRuleSet {
    pub rules: Vec<DomainRule>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainRule {
    pub id: String,
    pub patterns: Vec<String>,
    #[serde(default)]
    pub te_associated: bool,
    pub architecture_hint: Option<ArchitectureHint>,
    pub warning_code: Option<String>,
    pub warning_message: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainRuleMatch {
    pub rule_id: String,
    pub sequence_id: String,
    pub domain_name: String,
}

impl Default for DomainRuleSet {
    fn default() -> Self {
        Self {
            rules: vec![
                DomainRule {
                    id: "ltr-retroelement".into(),
                    patterns: vec![
                        "integrase".into(),
                        "gag".into(),
                        "retroviral protease".into(),
                    ],
                    te_associated: true,
                    architecture_hint: Some(ArchitectureHint::LtrLike),
                    warning_code: None,
                    warning_message: None,
                },
                DomainRule {
                    id: "tir-transposase".into(),
                    patterns: vec!["transposase".into(), "dde".into(), "d d e".into()],
                    te_associated: true,
                    architecture_hint: Some(ArchitectureHint::TirLike),
                    warning_code: None,
                    warning_message: None,
                },
                DomainRule {
                    id: "non-ltr-retroelement".into(),
                    patterns: vec![
                        "reverse transcript".into(),
                        "rnase h".into(),
                        "ap endonuclease".into(),
                        "endonuclease".into(),
                    ],
                    te_associated: true,
                    architecture_hint: Some(ArchitectureHint::LineLike),
                    warning_code: None,
                    warning_message: None,
                },
                DomainRule {
                    id: "rolling-circle".into(),
                    patterns: vec![
                        "rolling circle".into(),
                        "rep protein".into(),
                        "helicase".into(),
                    ],
                    te_associated: true,
                    architecture_hint: Some(ArchitectureHint::Unknown),
                    warning_code: None,
                    warning_message: None,
                },
                DomainRule {
                    id: "tyrosine-recombinase".into(),
                    patterns: vec!["tyrosine recombinase".into()],
                    te_associated: true,
                    architecture_hint: Some(ArchitectureHint::Unknown),
                    warning_code: None,
                    warning_message: None,
                },
                DomainRule {
                    id: "protein-primed-polymerase".into(),
                    patterns: vec!["protein-primed".into(), "protein primed".into()],
                    te_associated: true,
                    architecture_hint: Some(ArchitectureHint::Unknown),
                    warning_code: None,
                    warning_message: None,
                },
            ],
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainHit {
    pub sequence_id: String,
    pub source_type: String,
    pub domain_name: String,
    pub accession: String,
    pub description: String,
    pub evalue: f64,
    pub conditional_evalue: f64,
    pub bitscore: f64,
    pub bias: f64,
    pub aa_start: usize,
    pub aa_end: usize,
    pub nucleotide_start: Option<usize>,
    pub nucleotide_end: Option<usize>,
    pub strand: Option<String>,
    pub frame: Option<i8>,
    pub coverage: f64,
    pub partial: bool,
}
pub fn parse_domtblout(path: &Path, source_type: &str) -> Result<Vec<DomainHit>, CurationError> {
    let file = File::open(path).map_err(|e| io_error(path, e))?;
    let mut out = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|e| io_error(path, e))?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<_> = line.split_whitespace().collect();
        if f.len() < 23 {
            return Err(CurationError::ExternalTool(format!(
                "malformed domtblout line: {line}"
            )));
        }
        let parse = |i: usize| {
            f[i].parse::<f64>().map_err(|_| {
                CurationError::ExternalTool(format!("invalid domtblout number '{}'", f[i]))
            })
        };
        let aa_start = f[17]
            .parse::<usize>()
            .map_err(|_| CurationError::ExternalTool("invalid domain start".into()))?;
        let aa_end = f[18]
            .parse::<usize>()
            .map_err(|_| CurationError::ExternalTool("invalid domain end".into()))?;
        let accession = f[1].to_string();
        let domain_name = f[0].to_string();
        let query = f[3].to_string();
        let model_len = f[2]
            .parse::<usize>()
            .map_err(|_| CurationError::ExternalTool("invalid HMM model length".into()))?;
        let query_len = f[5]
            .parse::<usize>()
            .map_err(|_| CurationError::ExternalTool("invalid protein query length".into()))?;
        let hmm_start = f[15]
            .parse::<usize>()
            .map_err(|_| CurationError::ExternalTool("invalid HMM alignment start".into()))?;
        let hmm_end = f[16]
            .parse::<usize>()
            .map_err(|_| CurationError::ExternalTool("invalid HMM alignment end".into()))?;
        if hmm_end < hmm_start {
            return Err(CurationError::ExternalTool(
                "HMM alignment end precedes start".into(),
            ));
        }
        let coverage = if model_len > 0 {
            (hmm_end - hmm_start + 1) as f64 / model_len as f64
        } else {
            0.0
        };
        out.push(DomainHit {
            sequence_id: query,
            source_type: source_type.into(),
            domain_name,
            accession,
            description: f[22..].join(" "),
            evalue: parse(6)?,
            conditional_evalue: parse(11)?,
            bitscore: parse(13)?,
            bias: parse(14)?,
            aa_start,
            aa_end,
            nucleotide_start: None,
            nucleotide_end: None,
            strand: None,
            frame: None,
            coverage,
            partial: coverage < 0.8 || aa_start <= 2 || aa_end + 1 >= query_len,
        })
    }
    Ok(out)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TsdEvidence {
    pub copy_id: String,
    pub candidate_length: usize,
    pub left_sequence: String,
    pub right_sequence: String,
    pub matches: usize,
    pub mismatches: usize,
    pub supported: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CurationWarning {
    pub code: String,
    pub message: String,
    pub severity: String,
}
impl CurationWarning {
    fn new(code: &str, message: &str, severity: &str) -> Self {
        Self {
            code: code.into(),
            message: message.into(),
            severity: severity.into(),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ConfidenceTier {
    High,
    Moderate,
    Low,
    Insufficient,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputFile {
    pub path: String,
    pub sha256: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentSummary {
    pub copies: usize,
    pub columns: usize,
    pub usable_copies: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocusDiagnostic {
    pub code: String,
    pub copy_ids: Vec<String>,
    pub subject_id: Option<String>,
    pub overlap_fraction: Option<f64>,
    pub message: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CurationParameters {
    pub threads: usize,
    pub minimum_copies: usize,
    pub redundancy_control: bool,
    pub near_identity: f64,
    pub minimum_pairwise_overlap: usize,
    pub scan_selected_copies: bool,
    pub max_copy_domain_targets: usize,
    pub domain_rule_source: String,
    pub boundary: BoundaryConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryOverride {
    pub schema_version: String,
    pub source_schema_version: String,
    pub left_1based: usize,
    pub right_1based_inclusive: usize,
    #[serde(default)]
    pub notes: String,
    pub created_at: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CurationResult {
    pub schema_version: String,
    pub reputils_version: String,
    pub inputs: BTreeMap<String, InputFile>,
    pub parameters: CurationParameters,
    pub alignment_summary: AlignmentSummary,
    pub copy_metrics: Vec<CopyMetrics>,
    pub selected_copies: Vec<String>,
    pub boundary_result: BoundaryResult,
    pub curator_override: Option<BoundaryOverride>,
    pub refined_consensus: RefinedConsensus,
    pub structural_evidence: Vec<StructuralEvidence>,
    pub tsd_evidence: Vec<TsdEvidence>,
    pub tsd_support_by_length: BTreeMap<usize, usize>,
    pub domain_evidence: Vec<DomainHit>,
    pub domain_rule_matches: Vec<DomainRuleMatch>,
    pub locus_diagnostics: Vec<LocusDiagnostic>,
    pub subfamily_evidence: SubfamilyEvidence,
    pub architecture_hints: Vec<ArchitectureHint>,
    pub warnings: Vec<CurationWarning>,
    pub confidence_tier: ConfidenceTier,
    pub output_files: BTreeMap<String, String>,
    pub external_commands: Vec<String>,
    pub external_tool_versions: BTreeMap<String, String>,
}

#[derive(Debug, Clone)]
pub struct CurateOptions {
    pub alignment: PathBuf,
    pub out: PathBuf,
    pub manifest: Option<PathBuf>,
    pub query: Option<PathBuf>,
    pub genome: Option<PathBuf>,
    pub coordinates: Option<PathBuf>,
    pub hmm_db: Option<PathBuf>,
    pub repeatpeps: Option<PathBuf>,
    pub threads: usize,
    pub hmmscan: Option<PathBuf>,
    pub hmmsearch: Option<PathBuf>,
    pub getorf: Option<PathBuf>,
    pub domain_rules: Option<PathBuf>,
    pub boundary_override: Option<PathBuf>,
    pub scan_selected_copies: bool,
    pub max_copy_domain_targets: usize,
    pub min_copies: usize,
    pub force: bool,
}

pub fn run(options: &CurateOptions) -> Result<CurationResult, CurationError> {
    if options.threads == 0 || options.min_copies == 0 || options.max_copy_domain_targets == 0 {
        return Err(CurationError::InvalidParameter(
            "--threads, --min-copies, and --max-copy-domain-targets must be positive".into(),
        ));
    }
    if options.getorf.is_some() {
        return Err(CurationError::InvalidParameter("--getorf is reserved for a later ORF-target mode; omit it to use the built-in six-frame translation".into()));
    }
    prepare_output(&options.out, options.force)?;
    let alignment = Alignment::from_fasta(&options.alignment, options.min_copies)?;
    let mut metadata = if let Some(p) = &options.manifest {
        parse_manifest(p)?
    } else {
        HashMap::new()
    };
    if let Some(p) = &options.coordinates {
        metadata.extend(parse_coordinates_tsv(p)?);
    }
    let provisional = if let Some(p) = &options.query {
        Some(read_first_fasta(p)?)
    } else {
        None
    };
    let mut stats = column_statistics(&alignment);
    let mut metrics = copy_metrics(&alignment, provisional.as_deref(), &metadata);
    let cfg = BoundaryConfig {
        minimum_body_length: 50.min(alignment.width()).max(1),
        minimum_copies: options.min_copies,
        ..BoundaryConfig::default()
    };
    let mut boundaries = infer_boundaries(&stats, &cfg, alignment.len())?;
    let curator_override = options
        .boundary_override
        .as_deref()
        .map(|path| read_boundary_override(path, alignment.width()))
        .transpose()?;
    if let Some(value) = &curator_override {
        boundaries.left = value.left_1based - 1;
        boundaries.right = value.right_1based_inclusive;
    }
    let selected = select_representatives(
        &alignment,
        &mut metrics,
        0.25,
        0.98,
        50.min(alignment.width()),
        true,
    );
    let mut consensus = refined_consensus(
        &alignment,
        &selected,
        &metrics,
        boundaries.left,
        boundaries.right,
    );
    for p in &consensus.profile {
        stats[p.alignment_position].consensus_base = Some(p.base);
        stats[p.alignment_position].consensus_confidence = Some(p.confidence)
    }
    for (m, copy) in metrics.iter_mut().zip(alignment.copies()) {
        m.identity_to_refined_consensus = sequence_identity(
            &copy.sequence[boundaries.left..boundaries.right],
            consensus.sequence.as_bytes(),
            1,
        )
        .map(|x| x.0)
    }
    let mut structural = terminal_repeat_evidence(&consensus.sequence, 4, 50);
    if let Some(e) = structural.first() {
        if e.identity >= 0.75 {
            boundaries.architecture_hint = Some(if e.repeat_type == "TIR" {
                ArchitectureHint::TirLike
            } else {
                ArchitectureHint::LtrLike
            });
            boundaries.winning_score.structural_support = e.identity;
            boundaries.winning_score.total += 0.10 * e.identity
        }
    }
    let subfamilies = subfamily_analysis(&alignment, 50.min(alignment.width()));
    let mut warnings = boundaries.warnings.clone();
    if curator_override.is_some() {
        warnings.push(CurationWarning::new("curator_boundary_override_applied", "An explicit curator boundary override was applied and automatic joint reranking was disabled", "info"));
    }
    if subfamilies.candidate_clusters > 1 && subfamilies.support_for_splitting > 0.35 {
        warnings.push(CurationWarning::new(
            "possible_mixed_subfamilies",
            "Pairwise called-site distances support more than one candidate cluster",
            "warning",
        ))
    }
    if subfamilies.insufficient_pairs > 0 {
        warnings.push(CurationWarning::new(
            "insufficient_overlap_for_subfamily_analysis",
            "Some copy pairs lacked enough mutually called sites",
            "info",
        ))
    }
    if metadata.is_empty() && (options.genome.is_some() || options.coordinates.is_some()) {
        warnings.push(CurationWarning::new("coordinate_evidence_unavailable","Genome/coordinates were supplied without resolvable manifest metadata; coordinate-aware TSD analysis was not run","warning"))
    }
    let tsd_evidence = if let Some(genome) = &options.genome {
        coordinate_tsd_evidence(genome, &metadata, 2, 20)?
    } else {
        Vec::new()
    };
    let mut tsd_support_by_length = BTreeMap::new();
    for hit in tsd_evidence.iter().filter(|x| x.supported) {
        *tsd_support_by_length
            .entry(hit.candidate_length)
            .or_insert(0) += 1;
    }
    let (domain_rules, domain_rule_source) = load_domain_rules(options.domain_rules.as_deref())?;
    let domain_targets =
        build_domain_targets(options, &alignment, &selected, &boundaries, &consensus);
    let (mut domains, mut commands, mut versions) = run_domains(options, &domain_targets)?;
    if curator_override.is_none() {
        let preliminary_boundaries = (boundaries.left, boundaries.right);
        rescore_boundaries(
            JointBoundaryEvidence {
                alignment: &alignment,
                stats: &stats,
                metrics: &metrics,
                selected: &selected,
                tsds: &tsd_evidence,
                domains: &domains,
                domain_reference_left: preliminary_boundaries.0,
            },
            &mut boundaries,
        );
        if (boundaries.left, boundaries.right) != preliminary_boundaries {
            warnings.push(CurationWarning::new(
            "boundary_changed_by_joint_scoring",
            "Structural, TSD, domain, and conserved-edge evidence changed the alignment-only boundary proposal",
            "info",
        ));
            consensus = refined_consensus(
                &alignment,
                &selected,
                &metrics,
                boundaries.left,
                boundaries.right,
            );
            for column in &mut stats {
                column.consensus_base = None;
                column.consensus_confidence = None;
            }
            for position in &consensus.profile {
                stats[position.alignment_position].consensus_base = Some(position.base);
                stats[position.alignment_position].consensus_confidence = Some(position.confidence);
            }
            for (metric, copy) in metrics.iter_mut().zip(alignment.copies()) {
                metric.identity_to_refined_consensus = sequence_identity(
                    &copy.sequence[boundaries.left..boundaries.right],
                    consensus.sequence.as_bytes(),
                    1,
                )
                .map(|value| value.0);
            }
            structural = terminal_repeat_evidence(&consensus.sequence, 4, 50);
            let final_targets =
                build_domain_targets(options, &alignment, &selected, &boundaries, &consensus);
            let (final_domains, final_commands, final_versions) =
                run_domains(options, &final_targets)?;
            domains = final_domains;
            commands.extend(final_commands);
            versions.extend(final_versions);
        }
    } else {
        rescore_preserving_explicit_boundary(
            JointBoundaryEvidence {
                alignment: &alignment,
                stats: &stats,
                metrics: &metrics,
                selected: &selected,
                tsds: &tsd_evidence,
                domains: &domains,
                domain_reference_left: boundaries.left,
            },
            &mut boundaries,
        );
    }
    if domains.iter().any(|d| d.partial) {
        warnings.push(CurationWarning::new(
            "partial_domain_evidence",
            "One or more HMM domains are partial",
            "info",
        ));
    }
    let (domain_hints, domain_warnings, domain_rule_matches) =
        interpret_domains(&domains, &domain_rules);
    warnings.extend(domain_warnings);
    let (locus_diagnostics, locus_warnings) =
        locus_and_flank_diagnostics(&metadata, &stats, &boundaries);
    warnings.extend(locus_warnings);
    let tier = confidence_tier(
        selected.len(),
        &boundaries,
        &subfamilies,
        &stats,
        &structural,
        &domains,
        &domain_hints,
    );
    let mut outputs = BTreeMap::new();
    for name in [
        "curated.consensus.fa",
        "curated.alignment.fa",
        "boundaries.tsv",
        "copies.tsv",
        "columns.tsv",
        "structural-evidence.tsv",
        "domains.tsv",
        "warnings.tsv",
        "result.json",
        "report.html",
    ] {
        outputs.insert(name.into(), options.out.join(name).display().to_string());
    }
    let mut inputs = BTreeMap::new();
    add_input(&mut inputs, "alignment", &options.alignment);
    if let Some(p) = &options.manifest {
        add_input(&mut inputs, "manifest", p)
    }
    if let Some(p) = &options.query {
        add_input(&mut inputs, "query", p)
    }
    if let Some(p) = &options.genome {
        add_input(&mut inputs, "genome", p)
    }
    if let Some(p) = &options.coordinates {
        add_input(&mut inputs, "coordinates", p)
    }
    if let Some(p) = &options.boundary_override {
        add_input(&mut inputs, "boundary_override", p)
    }
    let mut hints: Vec<ArchitectureHint> =
        boundaries.architecture_hint.clone().into_iter().collect();
    if domains.is_empty()
        && structural
            .first()
            .is_some_and(|hit| hit.identity >= 0.80 && hit.candidate_length >= 8)
        && !hints.contains(&ArchitectureHint::NonAutonomous)
    {
        hints.push(ArchitectureHint::NonAutonomous);
    }
    for hint in domain_hints {
        if !hints.contains(&hint) {
            hints.push(hint);
        }
    }
    let result = CurationResult {
        schema_version: SCHEMA_VERSION.into(),
        reputils_version: env!("CARGO_PKG_VERSION").into(),
        inputs,
        parameters: CurationParameters {
            threads: options.threads,
            minimum_copies: options.min_copies,
            redundancy_control: true,
            near_identity: 0.98,
            minimum_pairwise_overlap: 50.min(alignment.width()),
            scan_selected_copies: options.scan_selected_copies,
            max_copy_domain_targets: options.max_copy_domain_targets,
            domain_rule_source,
            boundary: cfg,
        },
        alignment_summary: AlignmentSummary {
            copies: alignment.len(),
            columns: alignment.width(),
            usable_copies: selected.len(),
        },
        copy_metrics: metrics,
        selected_copies: selected
            .iter()
            .map(|&i| alignment.copies[i].id.clone())
            .collect(),
        boundary_result: boundaries,
        curator_override: curator_override.clone(),
        refined_consensus: consensus,
        structural_evidence: structural,
        tsd_evidence,
        tsd_support_by_length,
        domain_evidence: domains,
        domain_rule_matches,
        locus_diagnostics,
        subfamily_evidence: subfamilies,
        architecture_hints: hints,
        warnings,
        confidence_tier: tier,
        output_files: outputs,
        external_commands: commands,
        external_tool_versions: versions,
    };
    write_outputs(&options.out, &alignment, &stats, &result)?;
    Ok(result)
}

fn prepare_output(out: &Path, force: bool) -> Result<(), CurationError> {
    if out.exists() && !force {
        return Err(CurationError::OutputExists(out.display().to_string()));
    }
    fs::create_dir_all(out).map_err(|e| io_error(out, e))?;
    Ok(())
}
fn read_boundary_override(
    path: &Path,
    alignment_width: usize,
) -> Result<BoundaryOverride, CurationError> {
    let file = File::open(path).map_err(|e| io_error(path, e))?;
    let value: BoundaryOverride = serde_json::from_reader(file)?;
    if value.schema_version != "reputils.curation-override.v1"
        || value.source_schema_version != SCHEMA_VERSION
    {
        return Err(CurationError::InvalidParameter(
            "boundary override has an unsupported schema/source version".into(),
        ));
    }
    if value.left_1based == 0
        || value.right_1based_inclusive < value.left_1based
        || value.right_1based_inclusive > alignment_width
    {
        return Err(CurationError::InvalidParameter(format!(
            "boundary override must satisfy 1 <= left <= right <= {alignment_width}"
        )));
    }
    Ok(value)
}

fn read_first_fasta(path: &Path) -> Result<Vec<u8>, CurationError> {
    let reader = fasta::Reader::from_file(path).map_err(|e| CurationError::Fasta(e.to_string()))?;
    let record = reader
        .records()
        .next()
        .ok_or_else(|| CurationError::Fasta("query FASTA is empty".into()))?
        .map_err(|e| CurationError::Fasta(e.to_string()))?;
    Ok(record
        .seq()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect())
}
fn add_input(map: &mut BTreeMap<String, InputFile>, key: &str, path: &Path) {
    map.insert(
        key.into(),
        InputFile {
            path: path.display().to_string(),
            sha256: checksum(path).ok(),
        },
    );
}
fn checksum(path: &Path) -> Result<String, CurationError> {
    let data = fs::read(path).map_err(|e| io_error(path, e))?;
    Ok(format!("{:x}", Sha256::digest(data)))
}

fn parse_manifest(path: &Path) -> Result<HashMap<String, CopyMetadata>, CurationError> {
    let f = File::open(path).map_err(|e| io_error(path, e))?;
    let value: serde_json::Value = serde_json::from_reader(f)?;
    if value.get("schema_version").and_then(|x| x.as_str()) != Some("rep.copy-recovery.v1") {
        return Err(CurationError::Manifest(
            "expected schema_version rep.copy-recovery.v1".into(),
        ));
    }
    let Some(rows) = value.get("copies").and_then(|x| x.as_array()) else {
        let copies = path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join("copies.tsv");
        if copies.is_file() {
            return parse_coordinates_tsv(&copies);
        }
        return Err(CurationError::Manifest(format!(
            "family manifest has no inline copies array and sibling {} is missing",
            copies.display()
        )));
    };
    let mut out = HashMap::new();
    for row in rows {
        let id = row
            .get("copy_id")
            .or_else(|| row.get("id"))
            .and_then(|x| x.as_str())
            .ok_or_else(|| CurationError::Manifest("copy lacks copy_id/id".into()))?;
        let get_u = |k: &str| row.get(k).and_then(|x| x.as_u64());
        let get_s = |k: &str| row.get(k).and_then(|x| x.as_str()).map(str::to_string);
        out.insert(
            id.into(),
            CopyMetadata {
                subject_id: get_s("subject_id"),
                genomic_start: get_u("genomic_start"),
                genomic_end: get_u("genomic_end"),
                strand: get_s("strand"),
                locus_id: get_s("locus_id"),
                unextended_start: get_u("unextended_start").or_else(|| get_u("genomic_start")),
                unextended_end: get_u("unextended_end").or_else(|| get_u("genomic_end")),
                extended_start: get_u("extended_start").or_else(|| get_u("actual_extended_start")),
                extended_end: get_u("extended_end").or_else(|| get_u("actual_extended_end")),
            },
        );
    }
    Ok(out)
}

fn parse_coordinates_tsv(path: &Path) -> Result<HashMap<String, CopyMetadata>, CurationError> {
    let file = File::open(path).map_err(|e| io_error(path, e))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .ok_or_else(|| CurationError::Manifest("empty coordinates table".into()))?
        .map_err(|e| io_error(path, e))?;
    let names: Vec<&str> = header.split('\t').collect();
    let index = |name: &str| names.iter().position(|x| *x == name);
    let (id_i, locus_style) = if let Some(i) = index("copy_id") {
        (i, false)
    } else if let Some(i) = index("locus_id") {
        (i, true)
    } else {
        return Err(CurationError::Manifest(
            "coordinates table requires copy_id or locus_id".into(),
        ));
    };
    let mut out = HashMap::new();
    for line in lines {
        let line = line.map_err(|e| io_error(path, e))?;
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let get_s = |name: &str| {
            index(name)
                .and_then(|i| f.get(i))
                .filter(|x| !x.is_empty())
                .map(|x| x.to_string())
        };
        let get_u = |name: &str| get_s(name).and_then(|x| x.parse().ok());
        let raw_id = f
            .get(id_i)
            .ok_or_else(|| CurationError::Manifest("coordinate row lacks its ID field".into()))?;
        let id = if locus_style {
            format!("locus_id={raw_id}")
        } else {
            raw_id.to_string()
        };
        out.insert(
            id,
            CopyMetadata {
                subject_id: get_s("subject_id"),
                genomic_start: get_u("genomic_start"),
                genomic_end: get_u("genomic_end"),
                strand: get_s("strand"),
                locus_id: get_s("locus_id"),
                unextended_start: get_u("unextended_start").or_else(|| get_u("genomic_start")),
                unextended_end: get_u("unextended_end").or_else(|| get_u("genomic_end")),
                extended_start: get_u("extended_start").or_else(|| get_u("actual_extended_start")),
                extended_end: get_u("extended_end").or_else(|| get_u("actual_extended_end")),
            },
        );
    }
    Ok(out)
}
fn coordinate_tsd_evidence(
    genome_path: &Path,
    metadata: &HashMap<String, CopyMetadata>,
    min_len: usize,
    max_len: usize,
) -> Result<Vec<TsdEvidence>, CurationError> {
    let reader =
        fasta::Reader::from_file(genome_path).map_err(|e| CurationError::Fasta(e.to_string()))?;
    let mut genome = HashMap::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| CurationError::Fasta(e.to_string()))?;
        genome.insert(
            rec.id().to_string(),
            rec.seq()
                .iter()
                .map(|x| x.to_ascii_uppercase())
                .collect::<Vec<_>>(),
        );
    }
    let mut out = Vec::new();
    for (id, m) in metadata {
        let (subject, start, end) =
            match (m.subject_id.as_ref(), m.unextended_start, m.unextended_end) {
                (Some(s), Some(a), Some(b)) => (s, a as usize, b as usize),
                _ => continue,
            };
        let seq = match genome.get(subject) {
            Some(x) => x,
            None => continue,
        };
        if start > end || end > seq.len() {
            continue;
        }
        for len in min_len..=max_len {
            if start < len || end + len > seq.len() {
                continue;
            }
            let mut left = seq[start - len..start].to_vec();
            let mut right = seq[end..end + len].to_vec();
            if m.strand.as_deref() == Some("-") {
                let old = left;
                left = revcomp(&right);
                right = revcomp(&old)
            }
            let matches = left.iter().zip(&right).filter(|(a, b)| a == b).count();
            let mismatches = len - matches;
            out.push(TsdEvidence {
                copy_id: id.clone(),
                candidate_length: len,
                left_sequence: String::from_utf8_lossy(&left).into_owned(),
                right_sequence: String::from_utf8_lossy(&right).into_owned(),
                matches,
                mismatches,
                supported: mismatches <= if len >= 8 { 1 } else { 0 },
            })
        }
    }
    Ok(out)
}

fn translate_codon(c: &[u8]) -> u8 {
    match c {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'X',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}
#[derive(Debug, Clone)]
struct DomainNucleotideTarget {
    id: String,
    source_type: String,
    sequence: Vec<u8>,
}

fn six_frame_translations(target_key: &str, seq: &[u8]) -> Vec<(String, Vec<u8>)> {
    let strands = [('+', seq.to_vec()), ('-', revcomp(seq))];
    let mut out = Vec::new();
    for (strand, sequence) in strands {
        for frame in 0..3 {
            let aa = sequence[frame..]
                .chunks_exact(3)
                .map(translate_codon)
                .collect();
            out.push((
                format!("{target_key}|strand={strand}|frame={}", frame + 1),
                aa,
            ));
        }
    }
    out
}
fn map_domain_coordinates(hit: &mut DomainHit, nucleotide_len: usize) {
    let mut strand = None;
    let mut frame = None;
    for part in hit.sequence_id.split('|') {
        if let Some(v) = part.strip_prefix("strand=") {
            strand = Some(v.to_string())
        }
        if let Some(v) = part.strip_prefix("frame=") {
            frame = v.parse::<i8>().ok()
        }
    }
    if let (Some(s), Some(fr)) = (strand.as_ref(), frame) {
        let offset = (fr as usize) - 1;
        let (start, end) = if s == "+" {
            (
                offset + (hit.aa_start - 1) * 3,
                (offset + hit.aa_end * 3).min(nucleotide_len),
            )
        } else {
            (
                nucleotide_len.saturating_sub(offset + hit.aa_end * 3),
                nucleotide_len.saturating_sub(offset + (hit.aa_start - 1) * 3),
            )
        };
        hit.nucleotide_start = Some(start);
        hit.nucleotide_end = Some(end);
        hit.frame = Some(if s == "+" { fr } else { -fr });
        hit.strand = strand.clone();
    }
}

fn build_domain_targets(
    options: &CurateOptions,
    alignment: &Alignment,
    selected: &[usize],
    boundaries: &BoundaryResult,
    consensus: &RefinedConsensus,
) -> Vec<DomainNucleotideTarget> {
    let mut targets = vec![DomainNucleotideTarget {
        id: "refined_consensus".into(),
        source_type: "refined_consensus".into(),
        sequence: consensus.sequence.as_bytes().to_vec(),
    }];
    if options.scan_selected_copies {
        for &index in selected.iter().take(options.max_copy_domain_targets) {
            let copy = &alignment.copies()[index];
            let sequence = copy.sequence[boundaries.left..boundaries.right]
                .iter()
                .copied()
                .filter(|base| !is_gap(*base))
                .collect();
            targets.push(DomainNucleotideTarget {
                id: copy.id.clone(),
                source_type: "selected_copy".into(),
                sequence,
            });
        }
    }
    targets
}

type DomainRun = (Vec<DomainHit>, Vec<String>, BTreeMap<String, String>);

fn run_domains(
    options: &CurateOptions,
    targets: &[DomainNucleotideTarget],
) -> Result<DomainRun, CurationError> {
    let Some(db) = options.hmm_db.as_ref().or(options.repeatpeps.as_ref()) else {
        return Ok((Vec::new(), Vec::new(), BTreeMap::new()));
    };
    let (exe, tool_name) = if let Some(path) = options
        .hmmscan
        .clone()
        .or_else(|| env::var_os("REPUTILS_HMMSCAN").map(PathBuf::from))
    {
        (path, "hmmscan")
    } else if let Some(path) = options
        .hmmsearch
        .clone()
        .or_else(|| env::var_os("REPUTILS_HMMSEARCH").map(PathBuf::from))
    {
        (path, "hmmsearch")
    } else {
        (PathBuf::from("hmmscan"), "hmmscan")
    };
    let target = options.out.join("domain-target.fa");
    let mut f = File::create(&target).map_err(|e| io_error(&target, e))?;
    let mut target_metadata = HashMap::new();
    for (index, nucleotide) in targets.iter().enumerate() {
        let key = format!("target{index}");
        for (id, protein) in six_frame_translations(&key, &nucleotide.sequence) {
            target_metadata.insert(
                id.clone(),
                (
                    nucleotide.id.clone(),
                    nucleotide.source_type.clone(),
                    nucleotide.sequence.len(),
                ),
            );
            writeln!(f, ">{}\n{}", id, String::from_utf8_lossy(&protein))
                .map_err(|e| io_error(&target, e))?;
        }
    }
    let dom = options.out.join("domains.domtblout");
    let command = format!(
        "{} --domtblout {} --cpu {} {} {}",
        exe.display(),
        dom.display(),
        options.threads,
        db.display(),
        target.display()
    );
    let status = Command::new(&exe)
        .arg("--domtblout")
        .arg(&dom)
        .arg("--cpu")
        .arg(options.threads.to_string())
        .arg(db)
        .arg(&target)
        .status()
        .map_err(|e| {
            CurationError::ExternalTool(format!("could not execute {}: {e}", exe.display()))
        })?;
    if !status.success() {
        return Err(CurationError::ExternalTool(format!(
            "{tool_name} exited with {status}"
        )));
    }
    let mut hits = parse_domtblout(&dom, "unknown")?;
    for hit in &mut hits {
        if let Some((original_id, source_type, nucleotide_len)) =
            target_metadata.get(&hit.sequence_id)
        {
            map_domain_coordinates(hit, *nucleotide_len);
            hit.sequence_id.clone_from(original_id);
            hit.source_type.clone_from(source_type);
        }
    }
    let mut versions = BTreeMap::new();
    if let Ok(o) = Command::new(&exe).arg("-h").output() {
        let v = String::from_utf8_lossy(&o.stdout)
            .lines()
            .next()
            .unwrap_or("unknown")
            .to_string();
        versions.insert(tool_name.into(), v);
    }
    Ok((hits, vec![command], versions))
}

fn load_domain_rules(path: Option<&Path>) -> Result<(DomainRuleSet, String), CurationError> {
    match path {
        Some(path) => {
            let file = File::open(path).map_err(|e| io_error(path, e))?;
            let rules: DomainRuleSet = serde_json::from_reader(file)?;
            if rules.rules.is_empty() {
                return Err(CurationError::InvalidParameter(
                    "domain rule file contains no rules".into(),
                ));
            }
            for rule in &rules.rules {
                if rule.id.trim().is_empty()
                    || rule.patterns.is_empty()
                    || rule.patterns.iter().any(|p| p.trim().is_empty())
                {
                    return Err(CurationError::InvalidParameter(
                        "each domain rule needs a non-empty id and patterns".into(),
                    ));
                }
            }
            Ok((rules, path.display().to_string()))
        }
        None => Ok((DomainRuleSet::default(), "internal-default".into())),
    }
}

fn interpret_domains(
    hits: &[DomainHit],
    rules: &DomainRuleSet,
) -> (
    Vec<ArchitectureHint>,
    Vec<CurationWarning>,
    Vec<DomainRuleMatch>,
) {
    let mut hints = Vec::new();
    let mut warnings = Vec::new();
    let mut matches = Vec::new();
    let mut strands = HashSet::new();
    let mut te_like = false;
    for hit in hits {
        if let Some(strand) = &hit.strand {
            strands.insert(strand.clone());
        }
        let text = format!("{} {} {}", hit.domain_name, hit.accession, hit.description)
            .to_ascii_lowercase();
        for rule in &rules.rules {
            if rule
                .patterns
                .iter()
                .any(|pattern| text.contains(&pattern.to_ascii_lowercase()))
            {
                matches.push(DomainRuleMatch {
                    rule_id: rule.id.clone(),
                    sequence_id: hit.sequence_id.clone(),
                    domain_name: hit.domain_name.clone(),
                });
                te_like |= rule.te_associated;
                if let Some(hint) = &rule.architecture_hint {
                    if !hints.contains(hint) {
                        hints.push(hint.clone());
                    }
                }
                if let Some(code) = &rule.warning_code {
                    warnings.push(CurationWarning::new(
                        code,
                        rule.warning_message
                            .as_deref()
                            .unwrap_or("Configured domain rule emitted a warning"),
                        "warning",
                    ));
                }
            }
        }
    }
    matches.sort_by(|a, b| {
        (&a.rule_id, &a.sequence_id, &a.domain_name).cmp(&(
            &b.rule_id,
            &b.sequence_id,
            &b.domain_name,
        ))
    });
    matches.dedup_by(|a, b| {
        a.rule_id == b.rule_id && a.sequence_id == b.sequence_id && a.domain_name == b.domain_name
    });
    if strands.len() > 1 {
        warnings.push(CurationWarning::new(
            "domains_on_both_strands",
            "Domain evidence occurs on both translation strands",
            "warning",
        ));
    }
    if !hits.is_empty() && !te_like {
        warnings.push(CurationWarning::new("possible_host_gene", "Detected domains match no configured TE-associated rule; ordinary host sequence remains possible", "warning"));
    }
    let conflicting_architectures = hints
        .iter()
        .filter(|hint| !matches!(hint, ArchitectureHint::Unknown))
        .count()
        > 1;
    if conflicting_architectures {
        warnings.push(CurationWarning::new(
            "possible_nested_element",
            "Configured domain rules support more than one TE architecture; nested or chimeric sequence is possible",
            "warning",
        ));
        warnings.push(CurationWarning::new(
            "possible_chimeric_consensus",
            "Conflicting domain architectures occur in the curated sequence",
            "warning",
        ));
    }
    let consensus_hits: Vec<_> = hits
        .iter()
        .filter(|hit| hit.source_type == "refined_consensus")
        .collect();
    let reverse_transcriptase = consensus_hits.iter().find(|hit| {
        format!("{} {}", hit.domain_name, hit.description)
            .to_ascii_lowercase()
            .contains("reverse transcript")
    });
    let integrase = consensus_hits.iter().find(|hit| {
        format!("{} {}", hit.domain_name, hit.description)
            .to_ascii_lowercase()
            .contains("integrase")
    });
    if let (Some(rt), Some(int)) = (reverse_transcriptase, integrase) {
        if let (Some(rt_start), Some(int_start), Some(strand)) = (
            rt.nucleotide_start,
            int.nucleotide_start,
            rt.strand.as_deref(),
        ) {
            if (strand == "+" && int_start < rt_start) || (strand == "-" && int_start > rt_start) {
                warnings.push(CurationWarning::new("domain_order_inconsistent", "Integrase and reverse-transcriptase order is inconsistent with a simple LTR-like pol arrangement", "warning"));
            }
        }
    }
    if hints.contains(&ArchitectureHint::LtrLike) && reverse_transcriptase.is_none() {
        warnings.push(CurationWarning::new(
            "expected_domain_missing",
            "LTR-like domain evidence lacks detectable reverse transcriptase support",
            "info",
        ));
    }
    if hits.iter().any(|hit| hit.partial) {
        warnings.push(CurationWarning::new(
            "domain_truncated_by_boundary",
            "At least one domain is partial or reaches a translated target end",
            "warning",
        ));
    }
    (hints, warnings, matches)
}
/// Jointly rescore the alignment-derived boundary hypotheses.
///
/// `total = alignment score + 0.12 terminal support + 0.08 TSD support
///        + 0.10 domain containment - 0.20 domain cut - 0.10 conserved-edge penalty`.
/// Domain coordinates are mapped back to alignment space using `domain_reference_left`.
#[derive(Clone, Copy)]
struct JointBoundaryEvidence<'a> {
    alignment: &'a Alignment,
    stats: &'a [AlignmentColumnStats],
    metrics: &'a [CopyMetrics],
    selected: &'a [usize],
    tsds: &'a [TsdEvidence],
    domains: &'a [DomainHit],
    domain_reference_left: usize,
}

fn rescore_boundaries(evidence: JointBoundaryEvidence<'_>, boundary: &mut BoundaryResult) {
    let JointBoundaryEvidence {
        alignment,
        stats,
        metrics,
        selected,
        tsds,
        domains,
        domain_reference_left,
    } = evidence;
    let mut candidates = vec![BoundaryCandidate {
        left: boundary.left,
        right: boundary.right,
        score: boundary.winning_score.total,
        components: boundary.winning_score.clone(),
    }];
    candidates.extend(boundary.alternatives.clone());
    let tsd_support = tsds
        .iter()
        .filter(|hit| hit.supported)
        .map(|hit| hit.copy_id.as_str())
        .collect::<HashSet<_>>()
        .len() as f64
        / alignment.len().max(1) as f64;
    for candidate in &mut candidates {
        let consensus = refined_consensus(
            alignment,
            selected,
            metrics,
            candidate.left,
            candidate.right,
        );
        let terminal = terminal_repeat_evidence(&consensus.sequence, 4, 50);
        let structural_support = terminal
            .iter()
            .map(|hit| hit.identity * (hit.candidate_length as f64 / 20.0).min(1.0))
            .fold(0.0, f64::max);
        let mut contained = 0usize;
        let mut cut = 0usize;
        let relevant: Vec<_> = domains
            .iter()
            .filter(|hit| hit.source_type == "refined_consensus")
            .filter_map(|hit| {
                Some((
                    domain_reference_left + hit.nucleotide_start?,
                    domain_reference_left + hit.nucleotide_end?,
                ))
            })
            .collect();
        for (start, end) in &relevant {
            if *start >= candidate.left && *end <= candidate.right {
                contained += 1;
            } else if *start < candidate.right && *end > candidate.left {
                cut += 1;
            }
        }
        let domain_denominator = relevant.len().max(1) as f64;
        let domain_containment = contained as f64 / domain_denominator;
        let domain_cut_penalty = cut as f64 / domain_denominator;
        let edge_width = 25.min((candidate.right - candidate.left).max(1) / 4).max(1);
        let edge_support = stats[candidate.left..candidate.left + edge_width]
            .iter()
            .chain(&stats[candidate.right - edge_width..candidate.right])
            .map(|column| column.occupancy * column.major_base_frequency_among_called_bases)
            .sum::<f64>()
            / (2 * edge_width) as f64;
        let concentration = candidate.components.boundary_concentration.clamp(0.0, 1.0);
        let conserved_flank_penalty = edge_support * (1.0 - concentration);
        candidate.components.structural_support = structural_support;
        candidate.components.tsd_support = tsd_support;
        candidate.components.domain_containment = if relevant.is_empty() {
            0.0
        } else {
            domain_containment
        };
        candidate.components.domain_cut_penalty = domain_cut_penalty;
        candidate.components.conserved_flank_penalty = conserved_flank_penalty;
        let alignment_score = 0.65 * candidate.components.alignment_body_support
            + 0.20 * candidate.components.boundary_concentration
            + 0.15 * candidate.components.occupancy_transition
            - 0.20 * candidate.components.insufficient_copy_penalty;
        candidate.score = alignment_score
            + 0.12 * structural_support
            + 0.08 * tsd_support
            + 0.10 * candidate.components.domain_containment
            - 0.20 * domain_cut_penalty
            - 0.10 * conserved_flank_penalty;
        candidate.components.total = candidate.score;
    }
    candidates.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| a.left.cmp(&b.left))
            .then_with(|| a.right.cmp(&b.right))
    });
    candidates.dedup_by(|a, b| a.left == b.left && a.right == b.right);
    if let Some(winner) = candidates.first().cloned() {
        boundary.left = winner.left;
        boundary.right = winner.right;
        boundary.winning_score = winner.components;
        boundary.alternatives = candidates.into_iter().skip(1).take(9).collect();
        let raw: Vec<f64> = stats
            .iter()
            .map(|s| {
                0.50 * s.occupancy
                    + 0.35 * s.major_base_frequency_among_called_bases
                    + 0.15 * (1.0 - (s.entropy / 2.0).min(1.0))
            })
            .collect();
        let support = smooth(&raw, 11.min(stats.len()).max(1));
        boundary.left_confidence =
            boundary_left_score(boundary.left, stats, &support, alignment.len()).clamp(0.0, 1.0);
        boundary.right_confidence =
            boundary_right_score(boundary.right, stats, &support, alignment.len()).clamp(0.0, 1.0);
    }
}

fn rescore_preserving_explicit_boundary(
    evidence: JointBoundaryEvidence<'_>,
    boundary: &mut BoundaryResult,
) {
    let chosen = (boundary.left, boundary.right);
    let stats = evidence.stats;
    let alignment = evidence.alignment;
    rescore_boundaries(evidence, boundary);
    let automatic = BoundaryCandidate {
        left: boundary.left,
        right: boundary.right,
        score: boundary.winning_score.total,
        components: boundary.winning_score.clone(),
    };
    let mut all = vec![automatic];
    all.extend(boundary.alternatives.clone());
    if let Some(index) = all
        .iter()
        .position(|candidate| (candidate.left, candidate.right) == chosen)
    {
        let selected = all.remove(index);
        boundary.left = selected.left;
        boundary.right = selected.right;
        boundary.winning_score = selected.components;
        boundary.alternatives = all.into_iter().take(9).collect();
    }
    let raw: Vec<f64> = stats
        .iter()
        .map(|column| {
            0.50 * column.occupancy
                + 0.35 * column.major_base_frequency_among_called_bases
                + 0.15 * (1.0 - (column.entropy / 2.0).min(1.0))
        })
        .collect();
    let support = smooth(&raw, 11.min(stats.len()).max(1));
    boundary.left_confidence =
        boundary_left_score(boundary.left, stats, &support, alignment.len()).clamp(0.0, 1.0);
    boundary.right_confidence =
        boundary_right_score(boundary.right, stats, &support, alignment.len()).clamp(0.0, 1.0);
}

fn locus_and_flank_diagnostics(
    metadata: &HashMap<String, CopyMetadata>,
    stats: &[AlignmentColumnStats],
    boundary: &BoundaryResult,
) -> (Vec<LocusDiagnostic>, Vec<CurationWarning>) {
    let mut diagnostics = Vec::new();
    let rows: Vec<_> = metadata.iter().collect();
    for i in 0..rows.len() {
        for j in i + 1..rows.len() {
            let (id_a, a) = rows[i];
            let (id_b, b) = rows[j];
            if a.subject_id != b.subject_id || a.subject_id.is_none() {
                continue;
            }
            if let (Some(start_a), Some(end_a), Some(start_b), Some(end_b)) = (
                a.genomic_start.or(a.unextended_start),
                a.genomic_end.or(a.unextended_end),
                b.genomic_start.or(b.unextended_start),
                b.genomic_end.or(b.unextended_end),
            ) {
                let overlap = end_a.min(end_b).saturating_sub(start_a.max(start_b));
                let shorter = end_a
                    .saturating_sub(start_a)
                    .min(end_b.saturating_sub(start_b))
                    .max(1);
                let fraction = overlap as f64 / shorter as f64;
                if fraction >= 0.25 {
                    diagnostics.push(LocusDiagnostic {
                        code: "overlapping_source_loci".into(), copy_ids: vec![id_a.clone(), id_b.clone()],
                        subject_id: a.subject_id.clone(), overlap_fraction: Some(fraction),
                        message: "Copies derive from substantially overlapping source intervals and may not be independent".into(),
                    });
                }
            }
            if a.locus_id.is_some() && a.locus_id == b.locus_id {
                diagnostics.push(LocusDiagnostic {
                    code: "duplicate_locus_id".into(),
                    copy_ids: vec![id_a.clone(), id_b.clone()],
                    subject_id: a.subject_id.clone(),
                    overlap_fraction: None,
                    message: "Multiple alignment records resolve to the same locus ID".into(),
                });
            }
        }
    }
    let track = |slice: &[AlignmentColumnStats]| -> f64 {
        if slice.is_empty() {
            0.0
        } else {
            slice
                .iter()
                .map(|column| column.occupancy * column.major_base_frequency_among_called_bases)
                .sum::<f64>()
                / slice.len() as f64
        }
    };
    let left_outside = track(&stats[..boundary.left]);
    let right_outside = track(&stats[boundary.right..]);
    if boundary.left >= 25 && left_outside >= 0.55 {
        diagnostics.push(LocusDiagnostic {
            code: "conserved_sequence_outside_left_boundary".into(),
            copy_ids: Vec::new(),
            subject_id: None,
            overlap_fraction: None,
            message: "Substantial conserved alignment support remains left of the proposed element"
                .into(),
        });
    }
    if stats.len().saturating_sub(boundary.right) >= 25 && right_outside >= 0.55 {
        diagnostics.push(LocusDiagnostic {
            code: "conserved_sequence_outside_right_boundary".into(),
            copy_ids: Vec::new(),
            subject_id: None,
            overlap_fraction: None,
            message:
                "Substantial conserved alignment support remains right of the proposed element"
                    .into(),
        });
    }
    let span = boundary.right - boundary.left;
    if span >= 100 {
        let edge = 50.min(span / 4).max(1);
        if track(&stats[boundary.left..boundary.left + edge]) < 0.30 {
            diagnostics.push(LocusDiagnostic {
                code: "long_unstructured_left_extension".into(),
                copy_ids: Vec::new(),
                subject_id: None,
                overlap_fraction: None,
                message: "The retained left edge contains a long weakly supported extension".into(),
            });
        }
        if track(&stats[boundary.right - edge..boundary.right]) < 0.30 {
            diagnostics.push(LocusDiagnostic {
                code: "long_unstructured_right_extension".into(),
                copy_ids: Vec::new(),
                subject_id: None,
                overlap_fraction: None,
                message: "The retained right edge contains a long weakly supported extension"
                    .into(),
            });
        }
    }
    diagnostics.sort_by(|a, b| (&a.code, &a.copy_ids).cmp(&(&b.code, &b.copy_ids)));
    diagnostics.dedup_by(|a, b| a.code == b.code && a.copy_ids == b.copy_ids);
    let warnings = diagnostics
        .iter()
        .map(|diagnostic| {
            CurationWarning::new(
                &diagnostic.code,
                &diagnostic.message,
                if diagnostic.code == "overlapping_source_loci"
                    || diagnostic.code == "duplicate_locus_id"
                {
                    "warning"
                } else {
                    "info"
                },
            )
        })
        .collect();
    (diagnostics, warnings)
}

fn confidence_tier(
    usable: usize,
    boundary: &BoundaryResult,
    subfamilies: &SubfamilyEvidence,
    stats: &[AlignmentColumnStats],
    structural: &[StructuralEvidence],
    domains: &[DomainHit],
    domain_hints: &[ArchitectureHint],
) -> ConfidenceTier {
    if usable < 2 {
        return ConfidenceTier::Insufficient;
    }
    let missing = stats.iter().map(|column| column.gap_fraction).sum::<f64>() / stats.len() as f64;
    let structural_coherence = structural.first().map_or(0.4, |hit| {
        hit.identity * (hit.candidate_length as f64 / 20.0).min(1.0)
    });
    let domain_strands: HashSet<_> = domains
        .iter()
        .filter_map(|hit| hit.strand.as_deref())
        .collect();
    let conflicting_domain_hints = domain_hints
        .iter()
        .filter(|hint| !matches!(hint, ArchitectureHint::Unknown))
        .count()
        > 1;
    let domain_coherence = if domains.is_empty() {
        0.5
    } else if domain_strands.len() > 1 || conflicting_domain_hints {
        0.2
    } else {
        0.85
    };
    let boundary_agreement = (boundary.left_confidence + boundary.right_confidence) / 2.0;
    let score = (usable.min(6) as f64 / 6.0) * 0.25
        + boundary_agreement * 0.25
        + (1.0 - subfamilies.support_for_splitting) * 0.15
        + (1.0 - missing) * 0.15
        + structural_coherence * 0.10
        + domain_coherence * 0.10;
    if score >= 0.72 {
        ConfidenceTier::High
    } else if score >= 0.50 {
        ConfidenceTier::Moderate
    } else {
        ConfidenceTier::Low
    }
}
fn create(path: &Path) -> Result<File, CurationError> {
    File::create(path).map_err(|e| io_error(path, e))
}
fn write_outputs(
    out: &Path,
    alignment: &Alignment,
    stats: &[AlignmentColumnStats],
    r: &CurationResult,
) -> Result<(), CurationError> {
    let mut f = create(&out.join("curated.consensus.fa"))?;
    writeln!(
        f,
        ">curated_consensus schema={} confidence={:?}\n{}",
        SCHEMA_VERSION, r.confidence_tier, r.refined_consensus.sequence
    )
    .map_err(|e| io_error(&out.join("curated.consensus.fa"), e))?;
    let trimmed = alignment.trim(r.boundary_result.left, r.boundary_result.right);
    let mut f = create(&out.join("curated.alignment.fa"))?;
    for c in trimmed.copies {
        writeln!(f, ">{}\n{}", c.id, String::from_utf8_lossy(&c.sequence))
            .map_err(|e| io_error(&out.join("curated.alignment.fa"), e))?
    }
    write_boundaries(&out.join("boundaries.tsv"), &r.boundary_result)?;
    write_copies(&out.join("copies.tsv"), &r.copy_metrics)?;
    write_columns(&out.join("columns.tsv"), stats)?;
    write_structural(
        &out.join("structural-evidence.tsv"),
        &r.structural_evidence,
        &r.tsd_evidence,
    )?;
    write_domains(&out.join("domains.tsv"), &r.domain_evidence)?;
    write_warnings(&out.join("warnings.tsv"), &r.warnings)?;
    let jf = create(&out.join("result.json"))?;
    serde_json::to_writer_pretty(jf, r)?;
    let html = render_report(r, alignment, stats);
    fs::write(out.join("report.html"), html).map_err(|e| io_error(&out.join("report.html"), e))?;
    Ok(())
}
fn write_boundaries(p: &Path, b: &BoundaryResult) -> Result<(), CurationError> {
    let mut f = create(p)?;
    writeln!(f,"rank\tleft_1based\tright_1based_inclusive\tleft_confidence\tright_confidence\ttotal_score\talignment_body_support\tboundary_concentration\toccupancy_transition\tstructural_support\ttsd_support\tdomain_containment\tdomain_cut_penalty\tconserved_flank_penalty\tinsufficient_copy_penalty").map_err(|e|io_error(p,e))?;
    let mut all = vec![BoundaryCandidate {
        left: b.left,
        right: b.right,
        score: b.winning_score.total,
        components: b.winning_score.clone(),
    }];
    all.extend(b.alternatives.clone());
    for (i, x) in all.iter().enumerate() {
        writeln!(f,"{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",i+1,x.left+1,x.right,b.left_confidence,b.right_confidence,x.score,x.components.alignment_body_support,x.components.boundary_concentration,x.components.occupancy_transition,x.components.structural_support,x.components.tsd_support,x.components.domain_containment,x.components.domain_cut_penalty,x.components.conserved_flank_penalty,x.components.insufficient_copy_penalty).map_err(|e|io_error(p,e))?
    }
    Ok(())
}
fn opt<T: ToString>(x: &Option<T>) -> String {
    x.as_ref().map(ToString::to_string).unwrap_or_default()
}
fn write_copies(p: &Path, v: &[CopyMetrics]) -> Result<(), CurationError> {
    let mut f = create(p)?;
    writeln!(f,"copy_id\tfirst_non_gap_column\tlast_non_gap_column\taligned_span\tcalled_bases\tgap_bases\tunknown_bases\tambiguous_bases\tinternal_gap_runs\tterminal_gap_left\tterminal_gap_right\tidentity_to_provisional_consensus\tidentity_to_refined_consensus\tcovered_fraction\tquality_score\tselected_for_consensus\texclusion_reason\tsubject_id\tgenomic_start\tgenomic_end\tstrand\tlocus_id\tunextended_start\tunextended_end\textended_start\textended_end").map_err(|e|io_error(p,e))?;
    for x in v {
        writeln!(f,"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.5}\t{:.5}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",x.copy_id,opt(&x.first_non_gap_column.map(|z|z+1)),opt(&x.last_non_gap_column.map(|z|z+1)),x.aligned_span,x.called_bases,x.gap_bases,x.unknown_bases,x.ambiguous_bases,x.internal_gap_runs,x.terminal_gap_left,x.terminal_gap_right,opt(&x.identity_to_provisional_consensus),opt(&x.identity_to_refined_consensus),x.covered_fraction,x.quality_score,x.selected_for_consensus,opt(&x.exclusion_reason),opt(&x.metadata.subject_id),opt(&x.metadata.genomic_start),opt(&x.metadata.genomic_end),opt(&x.metadata.strand),opt(&x.metadata.locus_id),opt(&x.metadata.unextended_start),opt(&x.metadata.unextended_end),opt(&x.metadata.extended_start),opt(&x.metadata.extended_end)).map_err(|e|io_error(p,e))?
    }
    Ok(())
}
fn write_columns(p: &Path, v: &[AlignmentColumnStats]) -> Result<(), CurationError> {
    let mut f = create(p)?;
    writeln!(f,"position\tnon_gap_count\tcalled_base_count\tambiguous_count\tunknown_count\toccupancy\tgap_fraction\tn_fraction\tmajor_base\tmajor_base_frequency_among_called_bases\tmajor_base_frequency_among_all_sequences\tentropy\tcopy_start_count\tcopy_end_count\tgap_open_count\tgap_close_count\tconsensus_base\tconsensus_confidence").map_err(|e|io_error(p,e))?;
    for x in v {
        writeln!(f,"{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",x.position+1,x.non_gap_count,x.called_base_count,x.ambiguous_count,x.unknown_count,x.occupancy,x.gap_fraction,x.n_fraction,opt(&x.major_base),x.major_base_frequency_among_called_bases,x.major_base_frequency_among_all_sequences,x.entropy,x.copy_start_count,x.copy_end_count,x.gap_open_count,x.gap_close_count,opt(&x.consensus_base),opt(&x.consensus_confidence)).map_err(|e|io_error(p,e))?
    }
    Ok(())
}
fn write_structural(
    p: &Path,
    terminal: &[StructuralEvidence],
    tsds: &[TsdEvidence],
) -> Result<(), CurationError> {
    let mut f = create(p)?;
    writeln!(f,"evidence_type\tcopy_id\trepeat_type\tcandidate_length\talignment_score\tidentity\tcoverage\tleft_start\tleft_end\tright_start\tright_end\tstrand_relation\tsupporting_copies\tleft_sequence\tright_sequence\tmatches\tmismatches\tsupported").map_err(|e|io_error(p,e))?;
    for x in terminal {
        writeln!(
            f,
            "terminal_repeat\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t\t\t\t",
            x.repeat_type,
            x.candidate_length,
            x.alignment_score,
            x.identity,
            x.coverage,
            x.left_start + 1,
            x.left_end,
            x.right_start + 1,
            x.right_end,
            x.strand_relation,
            x.supporting_copies
        )
        .map_err(|e| io_error(p, e))?;
    }
    for x in tsds {
        writeln!(
            f,
            "coordinate_tsd\t{}\tTSD\t{}\t\t{}\t1.0\t\t\t\t\tdirect\t1\t{}\t{}\t{}\t{}\t{}",
            x.copy_id,
            x.candidate_length,
            x.matches as f64 / x.candidate_length as f64,
            x.left_sequence,
            x.right_sequence,
            x.matches,
            x.mismatches,
            x.supported
        )
        .map_err(|e| io_error(p, e))?;
    }
    Ok(())
}
fn write_domains(p: &Path, v: &[DomainHit]) -> Result<(), CurationError> {
    let mut f = create(p)?;
    writeln!(f,"sequence_id\tsource_type\tdomain_name\taccession\tdescription\tevalue\tconditional_evalue\tbitscore\tbias\taa_start\taa_end\tnucleotide_start\tnucleotide_end\tstrand\tframe\tcoverage\tpartial").map_err(|e|io_error(p,e))?;
    for x in v {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            x.sequence_id,
            x.source_type,
            x.domain_name,
            x.accession,
            x.description,
            x.evalue,
            x.conditional_evalue,
            x.bitscore,
            x.bias,
            x.aa_start,
            x.aa_end,
            opt(&x.nucleotide_start),
            opt(&x.nucleotide_end),
            opt(&x.strand),
            opt(&x.frame),
            x.coverage,
            x.partial
        )
        .map_err(|e| io_error(p, e))?
    }
    Ok(())
}
fn write_warnings(p: &Path, v: &[CurationWarning]) -> Result<(), CurationError> {
    let mut f = create(p)?;
    writeln!(f, "code\tseverity\tmessage").map_err(|e| io_error(p, e))?;
    for x in v {
        writeln!(f, "{}\t{}\t{}", x.code, x.severity, x.message).map_err(|e| io_error(p, e))?
    }
    Ok(())
}
fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
fn render_report(r: &CurationResult, a: &Alignment, stats: &[AlignmentColumnStats]) -> String {
    let tracks = stats
        .iter()
        .map(|s| {
            format!(
                "{{p:{},o:{:.4},i:{:.4},e:{:.4},s:{},z:{},go:{},gc:{}}}",
                s.position + 1,
                s.occupancy,
                s.major_base_frequency_among_called_bases,
                s.entropy,
                s.copy_start_count,
                s.copy_end_count,
                s.gap_open_count,
                s.gap_close_count
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let warnings = r
        .warnings
        .iter()
        .map(|w| {
            format!(
                "<li><code>{}</code> — {}</li>",
                esc(&w.code),
                esc(&w.message)
            )
        })
        .collect::<String>();
    let copies = r
        .copy_metrics
        .iter()
        .map(|c| {
            format!(
                "<tr><td>{}</td><td>{:.3}</td><td>{:.3}</td><td>{}</td><td>{}</td></tr>",
                esc(&c.copy_id),
                c.covered_fraction,
                c.quality_score,
                c.selected_for_consensus,
                esc(c.exclusion_reason.as_deref().unwrap_or(""))
            )
        })
        .collect::<String>();
    let alignment = a
        .copies
        .iter()
        .map(|c| {
            format!(
                "<div><b>{}</b> <code>{}</code></div>",
                esc(&c.id),
                esc(&String::from_utf8_lossy(&c.sequence))
            )
        })
        .collect::<String>();
    format!(
        r#"<!doctype html><html><head><meta charset="utf-8"><title>reputils curation report</title><style>body{{font:14px system-ui;max-width:1200px;margin:auto;padding:2rem}}code{{overflow-wrap:anywhere}}table{{border-collapse:collapse}}td,th{{border:1px solid #bbb;padding:.3rem}}.track{{height:130px;width:100%;border:1px solid #bbb}}.boundary{{background:#fff3cd;padding:1rem}}.alignment{{overflow:auto;white-space:nowrap;font-family:monospace}}.muted{{color:#555}}</style></head><body><h1>TE family curation report</h1><p class="muted">Semi-automatic proposal; boundaries and architecture require curator review.</p><section class="boundary"><h2>Assessment: {:?}</h2><p>Proposed alignment interval: <b>{}–{}</b> (1-based inclusive); confidence left {:.2}, right {:.2}; architecture hint {:?}.</p><details><summary>Alternative boundary hypotheses</summary><pre>{}</pre></details></section><h2>Input summary</h2><p>{} copies × {} columns; {} selected representatives.</p><h2>Consensus</h2><code>{}</code><h2>Evidence tracks</h2><canvas id="tracks" class="track"></canvas><p>Blue occupancy; green called-base identity; orange normalized entropy. Raw start/end and gap-transition counts are embedded below.</p><details><summary>Track data</summary><pre id="track-data">[{}]</pre></details><h2>Copies</h2><table><tr><th>ID</th><th>covered</th><th>quality</th><th>selected</th><th>reason</th></tr>{}</table><h2>Alignment around proposed boundaries</h2><div class="alignment">{}</div><h2>Terminal evidence</h2><pre>{}</pre><h2>Protein domains</h2><pre>{}</pre><h2>Subfamily analysis</h2><pre>{}</pre><h2>Locus and host-flank diagnostics</h2><pre>{}</pre><h2>Warnings and limitations</h2><ul>{}</ul><h2>Curator override</h2><p>Edit the one-based inclusive boundaries and download a small override record. This does not silently rewrite the curation outputs.</p><label>Left <input id="override-left" type="number" min="1"></label> <label>Right <input id="override-right" type="number" min="1"></label><br><label>Notes<br><textarea id="override-notes" rows="4" cols="80"></textarea></label><br><button type="button" onclick="downloadOverride()">Download boundary override JSON</button><p id="override-status" class="muted"></p><h2>Parameters and provenance</h2><pre>{}</pre><script>const d={};const c=document.getElementById('tracks'),x=c.getContext('2d');c.width=c.clientWidth*devicePixelRatio;c.height=130*devicePixelRatio;function line(k,color,scale=1){{x.strokeStyle=color;x.beginPath();d.forEach((q,i)=>{{let xx=i/(d.length-1||1)*c.width,yy=c.height-q[k]/scale*c.height;(i?x.lineTo(xx,yy):x.moveTo(xx,yy))}});x.stroke()}}line('o','#1769aa');line('i','#2e7d32');line('e','#ef6c00',2);const boundaryText=document.querySelector('.boundary p').textContent;const boundaryMatch=boundaryText.match(/interval:\s*(\d+).*(\d+)\s*\(/);if(boundaryMatch){{document.getElementById('override-left').value=boundaryMatch[1];document.getElementById('override-right').value=boundaryMatch[2];}}function downloadOverride(){{const left=Number(document.getElementById('override-left').value),right=Number(document.getElementById('override-right').value),status=document.getElementById('override-status');if(!Number.isInteger(left)||!Number.isInteger(right)||left<1||right<left){{status.textContent='Invalid boundaries: require 1 <= left <= right.';return;}}const payload={{schema_version:'reputils.curation-override.v1',source_schema_version:'reputils.curation.v1',left_1based:left,right_1based_inclusive:right,notes:document.getElementById('override-notes').value,created_at:new Date().toISOString()}};const blob=new Blob([JSON.stringify(payload,null,2)+'\n'],{{type:'application/json'}}),url=URL.createObjectURL(blob),link=document.createElement('a');link.href=url;link.download='boundary-override.json';link.click();URL.revokeObjectURL(url);status.textContent='Override downloaded; apply it explicitly in a later curation step.';}}</script></body></html>"#,
        r.confidence_tier,
        r.boundary_result.left + 1,
        r.boundary_result.right,
        r.boundary_result.left_confidence,
        r.boundary_result.right_confidence,
        r.boundary_result.architecture_hint,
        esc(&serde_json::to_string_pretty(&r.boundary_result.alternatives).unwrap_or_default()),
        r.alignment_summary.copies,
        r.alignment_summary.columns,
        r.alignment_summary.usable_copies,
        esc(&r.refined_consensus.sequence),
        tracks,
        copies,
        alignment,
        esc(&serde_json::to_string_pretty(&serde_json::json!({"terminal_repeats": r.structural_evidence, "coordinate_tsds": r.tsd_evidence, "tsd_support_by_length": r.tsd_support_by_length})).unwrap_or_default()),
        esc(&serde_json::to_string_pretty(&serde_json::json!({"hits": r.domain_evidence, "matched_rules": r.domain_rule_matches})).unwrap_or_default()),
        esc(&serde_json::to_string_pretty(&r.subfamily_evidence).unwrap_or_default()),
        esc(&serde_json::to_string_pretty(&r.locus_diagnostics).unwrap_or_default()),
        warnings,
        esc(&serde_json::to_string_pretty(&r.parameters).unwrap_or_default()),
        tracks
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ali(rows: &[(&str, &str)]) -> Alignment {
        Alignment::from_records(
            rows.iter()
                .map(|(id, s)| AlignedCopy {
                    id: (*id).into(),
                    sequence: s.as_bytes().to_vec(),
                })
                .collect(),
            1,
        )
        .unwrap()
    }
    #[test]
    fn validates_alignment() {
        assert!(matches!(
            Alignment::from_records(vec![], 1),
            Err(CurationError::EmptyAlignment)
        ));
        assert!(matches!(
            Alignment::from_records(
                vec![
                    AlignedCopy {
                        id: "a".into(),
                        sequence: b"AC".to_vec()
                    },
                    AlignedCopy {
                        id: "b".into(),
                        sequence: b"A".to_vec()
                    }
                ],
                1
            ),
            Err(CurationError::UnequalLength { .. })
        ));
    }
    #[test]
    fn column_denominators_and_entropy() {
        let a = ali(&[("a", "A-NR"), ("b", "ACN-")]);
        let s = column_statistics(&a);
        assert_eq!(s[1].non_gap_count, 1);
        assert_eq!(s[2].unknown_count, 2);
        assert_eq!(s[3].ambiguous_count, 1);
        assert_eq!(s[0].major_base_frequency_among_called_bases, 1.0);
        assert_eq!(s[0].entropy, 0.0)
    }
    #[test]
    fn starts_ends_and_gap_metrics() {
        let a = ali(&[("a", "--AC-G--"), ("b", "---CTG--")]);
        let s = column_statistics(&a);
        assert_eq!(s[2].copy_start_count, 1);
        assert_eq!(s[3].copy_start_count, 1);
        let m = copy_metrics(&a, None, &HashMap::new());
        assert_eq!(m[0].terminal_gap_left, 2);
        assert_eq!(m[0].terminal_gap_right, 2);
        assert_eq!(m[0].internal_gap_runs, 1)
    }
    #[test]
    fn consensus_iupac_and_deletion() {
        let a = ali(&[("a", "A-C"), ("b", "AGC"), ("c", "CGC"), ("d", "CGC")]);
        let mut m = copy_metrics(&a, None, &HashMap::new());
        for x in &mut m {
            x.quality_score = 1.0
        }
        let c = refined_consensus(&a, &[0, 1, 2, 3], &m, 0, 3);
        assert_eq!(c.sequence, "MGC");
        assert!(c.profile[1].deletion_polymorphism)
    }
    #[test]
    fn missing_not_identity() {
        let a = ali(&[("a", "A---C"), ("b", "AT--C")]);
        assert_eq!(
            sequence_identity(&a.copies[0].sequence, &a.copies[1].sequence, 2),
            Some((1.0, 2))
        );
        assert!(sequence_identity(&a.copies[0].sequence, &a.copies[1].sequence, 3).is_none())
    }
    #[test]
    fn symmetric_boundary_inference() {
        let a = ali(&[
            ("a", "----ACGTACGT----"),
            ("b", "----ACGTACGT----"),
            ("c", "----ACGTACGT----"),
        ]);
        let s = column_statistics(&a);
        let b = infer_boundaries(
            &s,
            &BoundaryConfig {
                smoothing_window: 3,
                minimum_body_length: 6,
                minimum_occupancy: 0.5,
                minimum_copies: 3,
                extension: 0,
            },
            3,
        )
        .unwrap();
        assert!(b.left <= 5);
        assert!(b.right >= 11)
    }
    #[test]
    fn terminal_repeat_types() {
        let x = terminal_repeat_evidence("ACGTAAAACGT", 4, 4);
        assert!(x.iter().any(|e| e.repeat_type == "LTR"));
        let x = terminal_repeat_evidence("ACGAAAAATCGT", 4, 4);
        assert!(x.iter().any(|e| e.repeat_type == "TIR"))
    }
    #[test]
    fn manifest_v1() {
        let p = std::env::temp_dir().join(format!("reputils-manifest-{}.json", std::process::id()));
        fs::write(&p,r#"{"schema_version":"rep.copy-recovery.v1","copies":[{"copy_id":"x","subject_id":"chr1","unextended_start":10}]}"#).unwrap();
        let m = parse_manifest(&p).unwrap();
        assert_eq!(m["x"].subject_id.as_deref(), Some("chr1"));
        let _ = fs::remove_file(p);
    }
    #[test]
    fn result_serializes() {
        let a = ali(&[("a", "ACGT"), ("b", "ACGT"), ("c", "ACGT")]);
        let s = column_statistics(&a);
        let b = infer_boundaries(
            &s,
            &BoundaryConfig {
                smoothing_window: 1,
                minimum_body_length: 2,
                minimum_occupancy: 0.1,
                minimum_copies: 3,
                extension: 0,
            },
            3,
        )
        .unwrap();
        assert_eq!(b.right, 4);
        assert!(serde_json::to_string(&b)
            .unwrap()
            .contains("alignment_body_support"))
    }
    #[test]
    fn parses_domtblout() {
        let p = std::env::temp_dir().join(format!("reputils-dom-{}.txt", std::process::id()));
        fs::write(&p,"RT PF00001.1 100 seq1 - 200 1e-20 50.0 0.0 1 1 1e-10 1e-12 45.0 0.0 1 90 5 95 3 98 0.95 reverse transcriptase\n").unwrap();
        let h = parse_domtblout(&p, "refined").unwrap();
        assert_eq!(h[0].aa_start, 5);
        assert_eq!(h[0].domain_name, "RT");
        let _ = fs::remove_file(p);
    }
    #[test]
    fn redundancy_selection_is_deterministic() {
        let a = ali(&[("a", "ACGTACGT"), ("b", "ACGTACGT"), ("c", "TCGTACGA")]);
        let mut m = copy_metrics(&a, None, &HashMap::new());
        for x in &mut m {
            x.quality_score = 1.0
        }
        let selected = select_representatives(&a, &mut m, 0.0, 0.98, 4, true);
        assert_eq!(selected, vec![0, 2]);
        assert_eq!(
            m[1].exclusion_reason.as_deref(),
            Some("near_identical_redundancy")
        );
    }
    #[test]
    fn subfamily_warning_distances_use_called_overlap() {
        let a = ali(&[
            ("a", "AAAA----"),
            ("b", "AAAT----"),
            ("c", "----CCCC"),
            ("d", "----CCCT"),
        ]);
        let s = subfamily_analysis(&a, 3);
        assert!(s.insufficient_pairs >= 4);
        assert_eq!(s.candidate_clusters, 2);
    }
    #[test]
    fn coordinate_tsd_is_adjacent_not_internal() {
        let p = std::env::temp_dir().join(format!("reputils-genome-{}.fa", std::process::id()));
        fs::write(&p, ">chr1\nGGACGTAAAAAAACGTCC\n").unwrap();
        let mut m = HashMap::new();
        m.insert(
            "copy1".into(),
            CopyMetadata {
                subject_id: Some("chr1".into()),
                unextended_start: Some(6),
                unextended_end: Some(12),
                strand: Some("+".into()),
                ..CopyMetadata::default()
            },
        );
        let hits = coordinate_tsd_evidence(&p, &m, 4, 4).unwrap();
        assert!(hits[0].supported);
        assert_eq!(hits[0].left_sequence, "ACGT");
        m.get_mut("copy1").unwrap().unextended_start = Some(2);
        let false_hit = coordinate_tsd_evidence(&p, &m, 4, 4).unwrap();
        assert!(false_hit.is_empty());
        let _ = fs::remove_file(p);
    }
    #[test]
    fn maps_plus_and_minus_domain_coordinates() {
        let base = DomainHit {
            sequence_id: "refined_consensus|strand=+|frame=2".into(),
            source_type: "x".into(),
            domain_name: "x".into(),
            accession: "x".into(),
            description: String::new(),
            evalue: 0.0,
            conditional_evalue: 0.0,
            bitscore: 0.0,
            bias: 0.0,
            aa_start: 2,
            aa_end: 4,
            nucleotide_start: None,
            nucleotide_end: None,
            strand: None,
            frame: None,
            coverage: 0.0,
            partial: false,
        };
        let mut plus = base.clone();
        map_domain_coordinates(&mut plus, 30);
        assert_eq!(
            (plus.nucleotide_start, plus.nucleotide_end, plus.frame),
            (Some(4), Some(13), Some(2))
        );
        let mut minus = base;
        minus.sequence_id = "refined_consensus|strand=-|frame=2".into();
        map_domain_coordinates(&mut minus, 30);
        assert_eq!(
            (minus.nucleotide_start, minus.nucleotide_end, minus.frame),
            (Some(17), Some(26), Some(-2))
        );
    }
    #[test]
    fn six_frame_translation_has_six_targets() {
        let x = six_frame_translations("target", b"ATGAAATAG");
        assert_eq!(x.len(), 6);
        assert_eq!(x[0].1, b"MKX");
    }
    #[test]
    fn asymmetric_and_low_conservation_boundaries() {
        let a = ali(&[
            ("a", "----AAAACGTACCCC----"),
            ("b", "------TGCATGCCCC----"),
            ("c", "--------ACGTACCC----"),
            ("d", "----CGTATGCACCCC----"),
        ]);
        let s = column_statistics(&a);
        let b = infer_boundaries(
            &s,
            &BoundaryConfig {
                smoothing_window: 3,
                minimum_body_length: 6,
                minimum_occupancy: 0.4,
                minimum_copies: 3,
                extension: 0,
            },
            4,
        )
        .unwrap();
        assert!(b.right >= 16);
        assert!(b.right_confidence > b.left_confidence);
        assert!(b.right - b.left >= 6);
    }

    #[test]
    fn isolated_conserved_column_does_not_define_body() {
        let a = ali(&[
            ("a", "A---ACGTACGT---"),
            ("b", "A---TGCATGCA---"),
            ("c", "A---CGTACGTA---"),
        ]);
        let s = column_statistics(&a);
        let b = infer_boundaries(
            &s,
            &BoundaryConfig {
                smoothing_window: 3,
                minimum_body_length: 6,
                minimum_occupancy: 0.4,
                minimum_copies: 3,
                extension: 0,
            },
            3,
        )
        .unwrap();
        assert!(b.left >= 2);
        assert!(b.right >= 11);
    }

    #[cfg(unix)]
    #[test]
    fn mocked_hmmer_execution() {
        use std::os::unix::fs::PermissionsExt;
        let dir = std::env::temp_dir().join(format!("reputils-hmmer-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let exe = dir.join("mock-hmmscan");
        fs::write(&exe, r#"#!/bin/sh
if [ "$1" = "-h" ]; then echo mock-hmmscan-1.0; exit 0; fi
while [ "$1" != "--domtblout" ]; do shift; done
out=$2
printf 'RT PF00001.1 100 target0|strand=+|frame=1 - 3 1e-20 50.0 0.0 1 1 1e-10 1e-12 45.0 0.0 1 90 1 2 1 2 0.95 reverse transcriptase\n' > "$out"
"#).unwrap();
        let mut mode = fs::metadata(&exe).unwrap().permissions();
        mode.set_mode(0o755);
        fs::set_permissions(&exe, mode).unwrap();
        let db = dir.join("db.hmm");
        fs::write(&db, "").unwrap();
        let options = CurateOptions {
            alignment: PathBuf::new(),
            out: dir.clone(),
            manifest: None,
            query: None,
            genome: None,
            coordinates: None,
            hmm_db: Some(db),
            repeatpeps: None,
            threads: 1,
            hmmscan: Some(exe),
            hmmsearch: None,
            getorf: None,
            domain_rules: None,
            boundary_override: None,
            scan_selected_copies: false,
            max_copy_domain_targets: 10,
            min_copies: 1,
            force: true,
        };
        let (hits, commands, versions) = run_domains(
            &options,
            &[DomainNucleotideTarget {
                id: "refined_consensus".into(),
                source_type: "refined_consensus".into(),
                sequence: b"ATGAAATAG".to_vec(),
            }],
        )
        .unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].nucleotide_start, Some(0));
        assert_eq!(commands.len(), 1);
        assert!(versions.contains_key("hmmscan"));
        let _ = fs::remove_dir_all(dir);
    }
    #[test]
    fn parses_real_rep_family_manifest_layout() {
        let dir = std::env::temp_dir().join(format!("reputils-rep-layout-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let manifest = dir.join("manifest.json");
        fs::write(&manifest, r#"{"schema_version":"rep.copy-recovery.v1","family_id":"fam","status":"success","raw_hits":1,"merged_loci":1,"selected":1,"warnings":[],"stages":{}}"#).unwrap();
        fs::write(dir.join("copies.tsv"), "family_id\tlocus_id\tsubject_id\tgenomic_start\tgenomic_end\tstrand\trequested_extended_start\trequested_extended_end\tactual_extended_start\tactual_extended_end\nfamily\tlocus-000001\tchr1\t10\t20\t+\t0\t30\t2\t28\n").unwrap();
        let parsed = parse_manifest(&manifest).unwrap();
        let copy = &parsed["locus_id=locus-000001"];
        assert_eq!(copy.unextended_start, Some(10));
        assert_eq!(copy.extended_start, Some(2));
        let _ = fs::remove_dir_all(dir);
    }
    fn empty_components(score: f64) -> BoundaryScoreComponents {
        BoundaryScoreComponents {
            alignment_body_support: score,
            boundary_concentration: 0.5,
            occupancy_transition: 0.2,
            structural_support: 0.0,
            tsd_support: 0.0,
            domain_containment: 0.0,
            domain_cut_penalty: 0.0,
            conserved_flank_penalty: 0.0,
            insufficient_copy_penalty: 0.0,
            total: score,
        }
    }

    #[test]
    fn configurable_domain_rules_record_matches() {
        let rules = DomainRuleSet {
            rules: vec![DomainRule {
                id: "custom-kinase".into(),
                patterns: vec!["kinase".into()],
                te_associated: false,
                architecture_hint: None,
                warning_code: Some("custom_domain_warning".into()),
                warning_message: Some("Configured host-like domain".into()),
            }],
        };
        let hit = DomainHit {
            sequence_id: "consensus".into(),
            source_type: "refined_consensus".into(),
            domain_name: "Protein_kinase".into(),
            accession: "PF00069".into(),
            description: "kinase domain".into(),
            evalue: 1e-20,
            conditional_evalue: 1e-20,
            bitscore: 50.0,
            bias: 0.0,
            aa_start: 1,
            aa_end: 10,
            nucleotide_start: Some(0),
            nucleotide_end: Some(30),
            strand: Some("+".into()),
            frame: Some(1),
            coverage: 0.9,
            partial: false,
        };
        let (_, warnings, matches) = interpret_domains(&[hit], &rules);
        assert_eq!(matches[0].rule_id, "custom-kinase");
        assert!(warnings
            .iter()
            .any(|warning| warning.code == "custom_domain_warning"));
        assert!(warnings
            .iter()
            .any(|warning| warning.code == "possible_host_gene"));
    }

    #[test]
    fn selected_copy_domain_targets_are_bounded() {
        let alignment = ali(&[("a", "ACGTACGT"), ("b", "ACGTACGT"), ("c", "ACGTACGT")]);
        let options = CurateOptions {
            alignment: PathBuf::new(),
            out: PathBuf::new(),
            manifest: None,
            query: None,
            genome: None,
            coordinates: None,
            hmm_db: None,
            repeatpeps: None,
            threads: 1,
            hmmscan: None,
            hmmsearch: None,
            getorf: None,
            domain_rules: None,
            boundary_override: None,
            scan_selected_copies: true,
            max_copy_domain_targets: 2,
            min_copies: 1,
            force: false,
        };
        let boundary = BoundaryResult {
            left: 0,
            right: 8,
            left_confidence: 1.0,
            right_confidence: 1.0,
            architecture_hint: None,
            alternatives: Vec::new(),
            winning_score: empty_components(1.0),
            warnings: Vec::new(),
        };
        let consensus = RefinedConsensus {
            sequence: "ACGTACGT".into(),
            profile: Vec::new(),
        };
        let targets = build_domain_targets(&options, &alignment, &[0, 1, 2], &boundary, &consensus);
        assert_eq!(targets.len(), 3);
        assert_eq!(targets[2].id, "b");
        assert_eq!(targets[2].source_type, "selected_copy");
    }

    #[test]
    fn joint_scoring_can_choose_terminally_supported_alternative() {
        let alignment = ali(&[
            ("a", "ACGTAAAACGTGGGG"),
            ("b", "ACGTAAAACGTGGGG"),
            ("c", "ACGTAAAACGTGGGG"),
        ]);
        let stats = column_statistics(&alignment);
        let mut metrics = copy_metrics(&alignment, None, &HashMap::new());
        for metric in &mut metrics {
            metric.quality_score = 1.0;
        }
        let mut boundary = BoundaryResult {
            left: 0,
            right: 15,
            left_confidence: 0.5,
            right_confidence: 0.5,
            architecture_hint: None,
            alternatives: vec![BoundaryCandidate {
                left: 0,
                right: 11,
                score: 0.8,
                components: empty_components(0.8),
            }],
            winning_score: empty_components(0.8),
            warnings: Vec::new(),
        };
        rescore_boundaries(
            JointBoundaryEvidence {
                alignment: &alignment,
                stats: &stats,
                metrics: &metrics,
                selected: &[0, 1, 2],
                tsds: &[],
                domains: &[],
                domain_reference_left: 0,
            },
            &mut boundary,
        );
        assert_eq!(boundary.right, 11);
        assert!(boundary.winning_score.structural_support > 0.0);
    }

    #[test]
    fn locus_overlap_and_conserved_flank_are_diagnosed() {
        let alignment = ali(&[("a", "ACGTACGT"), ("b", "ACGTACGT")]);
        let stats = column_statistics(&alignment);
        let mut metadata = HashMap::new();
        metadata.insert(
            "a".into(),
            CopyMetadata {
                subject_id: Some("chr1".into()),
                genomic_start: Some(10),
                genomic_end: Some(30),
                ..CopyMetadata::default()
            },
        );
        metadata.insert(
            "b".into(),
            CopyMetadata {
                subject_id: Some("chr1".into()),
                genomic_start: Some(20),
                genomic_end: Some(40),
                ..CopyMetadata::default()
            },
        );
        let boundary = BoundaryResult {
            left: 0,
            right: 8,
            left_confidence: 1.0,
            right_confidence: 1.0,
            architecture_hint: None,
            alternatives: Vec::new(),
            winning_score: empty_components(1.0),
            warnings: Vec::new(),
        };
        let (diagnostics, warnings) = locus_and_flank_diagnostics(&metadata, &stats, &boundary);
        assert!(diagnostics
            .iter()
            .any(|value| value.code == "overlapping_source_loci"));
        assert!(warnings
            .iter()
            .any(|value| value.code == "overlapping_source_loci"));
    }

    #[test]
    fn boundary_override_schema_and_range_are_validated() {
        let path =
            std::env::temp_dir().join(format!("reputils-override-{}.json", std::process::id()));
        fs::write(&path, r#"{"schema_version":"reputils.curation-override.v1","source_schema_version":"reputils.curation.v1","left_1based":2,"right_1based_inclusive":7,"notes":"reviewed","created_at":null}"#).unwrap();
        let value = read_boundary_override(&path, 10).unwrap();
        assert_eq!((value.left_1based, value.right_1based_inclusive), (2, 7));
        assert!(read_boundary_override(&path, 6).is_err());
        let _ = fs::remove_file(path);
    }
}
