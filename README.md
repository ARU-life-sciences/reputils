# reputils

`reputils` is a semi-automatic curation tool for unknown transposable-element (TE) families. It combines evidence from an alignment of extended genomic copies to propose—not assert—family boundaries, select representative copies, build an uncertainty-aware consensus, and produce a reviewable HTML report plus versioned machine-readable results.

The original lightweight `con`, `div`, `dot`, `html`, `tir`, `tsd`, and `ttc` commands remain available. `ttc` is the legacy fixed-threshold trimmer; use `curate` when you need combined evidence, alternatives, warnings, and explicit confidence.

## Where reputils fits

Three related tasks should not be conflated:

1. **Copy recovery and alignment generation** finds genomic instances and creates an alignment. The companion `rep` program can provide extended copies, coordinates, and a `rep.copy-recovery.v1` manifest. `reputils` does not run BLAST, merge loci, or make MAFFT its primary pipeline.
2. **Curation** decides which aligned region is supported as the repeat, assesses copies, and refines the consensus. This is the role of `reputils curate`.
3. **Classification** interprets terminal structures and protein domains. `reputils` reports conservative architecture hints; it is not a fully automatic TE classifier or a replacement for RepeatMasker/RepeatModeler.

The minimum standalone input is an equal-length FASTA alignment containing at least three copies by default. A provisional consensus, manifest, coordinates, genome, and protein-domain HMMs are optional.

## Installation

A Rust toolchain is required.

```bash
git clone https://github.com/ARU-life-sciences/reputils
cd reputils
cargo build --release
./target/release/reputils --help
```

## Complete curation example

```bash
reputils curate \
  --alignment copies.aligned.fa \
  --manifest family.manifest.json \
  --query provisional-consensus.fa \
  --genome genome.fa \
  --hmm-db te-domains.hmm \
  --threads 8 \
  --out curated-family/
```

Only `--alignment` and `--out` are mandatory:

```bash
reputils curate --alignment copies.aligned.fa --out curated-family/
```

Options:

- `--manifest`: a `rep.copy-recovery.v1` JSON manifest. For a real `rep` family manifest, the sibling `copies.tsv` is discovered automatically and `locus_id=…` FASTA IDs are resolved; standalone manifests may instead contain an inline `copies` array.
- `--coordinates`: a tab-separated table whose header includes `copy_id`; recognised metadata columns are listed below. Values override matching manifest rows.
- `--query`: provisional consensus FASTA used for per-copy identity.
- `--genome`: source genome FASTA. With `subject_id`, `unextended_start`, and `unextended_end`, immediately adjacent flanks are tested for TSDs. Coordinates are interpreted as zero-based, half-open.
- `--hmm-db` or `--repeatpeps`: a local HMM database searched against protein targets derived from the refined consensus.
- `--hmmscan`, `--hmmsearch`: executable paths. `REPUTILS_HMMSCAN` and `REPUTILS_HMMSEARCH` are also recognised. `hmmscan` is preferred; `hmmsearch` is used when no hmmscan path/environment override is supplied.
- `--domain-rules`: optional JSON rule table; see `examples/domain-rules.json`. Each rule has an ID, case-insensitive patterns, an association category (`te-associated`, `host-enriched`, `ubiquitous`, or `unknown`), an optional architecture hint, and an optional warning. Legacy `te_associated` booleans remain accepted.
- `--scan-selected-copies` and `--max-copy-domain-targets`: opt in to HMM searches of a deterministic, bounded selected-copy subset in addition to the refined consensus.
- `--boundary-override`: apply a validated `reputils.curation-override.v1` file exported by the HTML report. Explicit curator boundaries are preserved and not automatically reranked.
- `--getorf`: optional EMBOSS `getorf` executable path; `REPUTILS_GETORF` is also recognised. When configured, ORFs from the refined consensus (and the bounded selected-copy targets, if requested) replace built-in six-frame translations. The invoked command and version are recorded. Without it, the deterministic six-frame fallback remains available.
- `--split-subfamilies`: conservatively emit cluster-specific alignments and consensuses under `subfamilies/`. Output requires agreement across linkage methods, threshold perturbation, copy resampling, distance separation, and cluster size; otherwise the one-family result is retained.
- `--min-subfamily-copies` and `--min-split-support`: configure the emitted-cluster minimum (default 3) and separation support (default 0.50). Split support below 0.50 is rejected.
- `--confidence-calibration`: optional labeled `reputils.confidence-calibration.v1` JSON table. A probability is reported only when the heuristic score falls in a validated bin derived from an externally labeled data set.
- `--threads`, `--min-copies`, `--force`: external-tool threads, validation minimum, and permission to overwrite known output files in an existing directory.

The command never silently discards a copy. Every copy appears in `copies.tsv` and `result.json` with selection state and exclusion reason.

## Input validation

The curation core returns structured errors for:

- an empty alignment;
- unequal aligned lengths;
- duplicate record IDs;
- unsupported nucleotide characters;
- fewer records than `--min-copies`;
- malformed manifests or coordinate tables;
- invalid parameters;
- requested domain analysis when an external executable fails or is unavailable.

Accepted symbols are A/C/G/T/U, IUPAC ambiguity symbols, N, `?`, `-`, and `.` (case-insensitive). U is treated as T for nucleotide calls; `-` and `.` are gaps. Gaps, unknowns, ambiguous bases, and called A/C/G/T bases remain distinct evidence categories.

## Evidence and algorithms

### Alignment columns

`columns.tsv` contains one row per input alignment column. Denominators are explicit:

- `occupancy = non_gap_count / all_sequences`; ambiguous and unknown symbols are occupied.
- `gap_fraction = gap_count / all_sequences`.
- `n_fraction = unknown_count / all_sequences`; this counts N and `?`, not other IUPAC ambiguity.
- `major_base_frequency_among_called_bases = major A/C/G/T count / called A/C/G/T count`.
- `major_base_frequency_among_all_sequences = major A/C/G/T count / all_sequences`.
- entropy is Shannon entropy in bits over called A/C/G/T only; gaps, N, and ambiguity do not become nucleotides.
- copy starts/ends are the first/last non-gap aligned columns. Gap opens/closes count transitions between adjacent columns within each row.

Consensus call and confidence fields are populated only for retained columns. Confidence reflects weighted support and copy count; it is not a calibrated probability.

### Copy metrics and representative selection

Completeness, called/gap/unknown counts, terminal gaps, internal gap runs, coverage, provisional/refined identity, and a transparent quality score are calculated for every row. Terminal missingness is not treated as an internal deletion. Pairwise identity uses mutually called A/C/G/T sites only.

Copies are sorted deterministically by quality. Clearly unusable fragments are marked excluded. With redundancy control, a copy at or above 98% identity to a selected copy across the minimum mutual overlap is retained in metadata but does not add another full vote to consensus calling.

### Boundary inference

The initial, deliberately understandable model smooths a per-column body track:

```text
body support = 0.50 × occupancy
             + 0.35 × called-base identity
             + 0.15 × (1 − entropy / 2)
```

It locates sustained supported intervals, searches neighbouring recurrent copy starts/ends, and ranks boundary pairs using:

```text
pair score = 0.65 × mean internal body support
           + 0.20 × boundary start/end concentration
           + 0.15 × occupancy-transition support
           − 0.20 × insufficient-copy penalty
```

All retained alternatives are then jointly rescored with `+0.12 × terminal support + 0.08 × coordinate-TSD support + 0.10 × domain containment − 0.20 × domain-cut penalty − 0.10 × conserved-edge penalty`. If the winner changes, the final consensus, copy identities, terminal evidence, and optional domain search are regenerated. Every winning and alternative component is written to `boundaries.tsv` and JSON. Left and right confidence describe agreement among signals, not a statistical posterior. Diffuse 5′ starts may produce a LINE-like/asymmetric hint without forcing both ends to behave symmetrically.

The `confidence_assessment.boundary_robustness` object reruns alignment-only inference under three sensitivity analyses:

- deterministic leave-one-copy-out inference (up to 20 copies), including endpoint shifts and influential copy IDs;
- 20 deterministic copy-bootstrap samples with replacement;
- 27 combinations of smoothing window, minimum occupancy, and minimum body length.

It reports endpoint recovery fractions and combined minimum, lower-decile, median, upper-decile, and maximum positions. These ranges measure algorithm/data sensitivity; they are not statistical confidence or credible intervals.

### Refined consensus

Only selected representatives vote, weighted by copy quality. Terminal missing sequence does not vote as a deletion. Internal gaps are tracked separately and columns with substantial gap/called support are marked as deletion polymorphisms. Strong single-base support produces A/C/G/T; genuine mixed support produces an IUPAC code; sparse or weak support produces N. The FASTA remains standard nucleotide sequence.

### Structural evidence

Terminal windows across several lengths are compared directly for LTR-like evidence and against the reverse complement for TIR-like evidence. Short weak matches cannot dominate boundary inference. `structural-evidence.tsv` contains terminal-repeat rows and, when coordinates and a genome are available, per-copy TSD rows.

Coordinate-aware TSD analysis compares only the immediately adjacent left and right flank for lengths 2–20 bp. Exact matches and near-exact matches (one mismatch for lengths of at least 8 bp) are distinguished from unsupported candidates. Family support is aggregated by length in `result.json`. Arbitrary internal matching substrings are not accepted as coordinate TSDs.

### Subfamilies

Distances use mutually called overlapping sites and a minimum overlap. The primary deterministic single-link membership is compared with average- and complete-linkage results, thresholds ±0.03, and leave-one-copy-out clustering. JSON reports method cluster counts, pairwise co-assignment agreement, threshold stability, copy-resampling stability, within/between distance, and split support.

By default this evidence only generates warnings. With `--split-subfamilies`, outputs require at least 0.90 method agreement, 0.80 threshold stability, 0.80 copy-resampling stability, equal cluster counts across all three linkage methods, sufficient separation, and the configured minimum cluster size/support. Each accepted cluster receives its own boundary inference and trimmed consensus; the primary family consensus is always retained.

### Locus and host-flank diagnostics

When coordinates are available, substantially overlapping source intervals and duplicate locus IDs are reported because they reduce independence. Conserved sequence outside the proposal and long weakly supported retained edges are emitted as structured `locus_diagnostics` plus warnings. These signals do not automatically delete sequence.

### Protein domains

`reputils` does not implement profile-HMM search. It either produces six-frame translations or invokes configured EMBOSS `getorf`, then executes HMMER, parses domain-table output, maps amino-acid hits back to zero-based half-open nucleotide coordinates on both strands, and retains partial hits. A successful `getorf` run that finds no ORFs yields an empty domain result rather than failing curation. Executed commands and tool versions are recorded.

Compatible partial hits to the same model can be chained across translation frames when they occur on one strand, in compatible nucleotide order, with no more than 150 nt between fragments and limited overlap. Chains require at least two frames and combined model coverage of 0.50. `domain-chains.tsv` distinguishes consensus-only chains from chains recurring at two or more independent selected-copy loci. Chaining is evidence of a possible frameshift; it never edits the nucleotide sequence.

Domain names, accessions, and descriptions are interpreted with a validated JSON-configurable rule table (or the extensible internal defaults) covering reverse transcriptase, RNase H, integrase, retroviral protease/gag, endonuclease, DDE/transposase, rolling-circle Rep, helicase, tyrosine recombinase, and protein-primed polymerase concepts. These yield hints and warnings, not automatic rejection.

The host assessment treats `host-enriched` separately from ubiquitous and unknown domains. Its component score records consensus host-domain coverage, host-hit fraction, recurrence across independent loci, boundary uncertainty, lack of terminal-structure counterevidence, and absence of TE-domain support. It distinguishes absent, mixed, review-level, recurrent, and strong host-like signals. No category automatically rejects or trims a family.

Supply databases yourself. Compatible sources may include Dfam profile HMMs, RepeatPeps-related local resources, and curated Pfam-derived TE models. Check each source's licence and conversion/indexing requirements; the repository does not redistribute or download third-party databases and does not claim that every source is freely redistributable.

## Outputs and schema

```text
curated-family/
├── curated.consensus.fa
├── curated.alignment.fa
├── boundaries.tsv
├── copies.tsv
├── columns.tsv
├── structural-evidence.tsv
├── domains.tsv
├── domain-chains.tsv
├── warnings.tsv
├── result.json
└── report.html
```

When a requested subfamily split passes the conservative checks, the directory additionally contains:

```text
curated-family/subfamilies/
├── cluster-1.consensus.fa
├── cluster-1.alignment.fa
├── cluster-2.consensus.fa
└── cluster-2.alignment.fa
```

`result.json` uses `schema_version: "reputils.curation.v1"`. Its top-level fields are:

- `schema_version`, `reputils_version`;
- `inputs` with paths and SHA-256 checksums where readable;
- `parameters` and alignment summary;
- all copy metrics and selected copy IDs;
- winning and alternative boundaries with component scores;
- refined consensus and per-position confidence/profile;
- terminal-repeat, TSD, TSD-by-length, domain, cross-frame chain, and subfamily evidence;
- architecture hints, warnings, overall confidence tier, factor scores, and boundary robustness;
- optional emitted subfamily consensuses and a host-domain assessment;
- output file paths, external commands, and external tool versions.

Rust/JSON coordinates are zero-based and half-open unless a field says otherwise. Curator-facing TSV boundary and sequence positions are one-based inclusive. Empty optional TSV fields mean unavailable, not zero.

The self-contained `report.html` includes a curator boundary editor that downloads a versioned override JSON, and the proposal, alternatives via the embedded result evidence, occupancy/identity/entropy tracks, raw start/end and gap-transition data, copy completeness and selection, the alignment, structural evidence, domains, subfamily analysis, warnings, parameters, and provenance. It uses no browser-side network resources. A curator should review and may override the proposed interval before downstream use. The downloaded `boundary-override.json` can be applied with `--boundary-override`; its one-based inclusive interval is range checked and recorded in `result.json`.

## Confidence and warnings

Overall tiers are `high`, `moderate`, `low`, or `insufficient`. A documented weighted heuristic combines usable-copy count, boundary agreement, completeness, mixed-subfamily evidence, terminal-structure coherence, and protein-domain coherence. Leave-one-out, bootstrap, and parameter-perturbation analyses report left/right stability separately. Neither the tier nor stability is a statistical probability.

An optional calibrated probability is emitted only with a validated external table such as:

```json
{
  "schema_version": "reputils.confidence-calibration.v1",
  "source": "description of the held-out, curator-labeled data set",
  "bins": [
    {
      "minimum_score": 0.0,
      "maximum_score": 0.5,
      "estimated_probability_correct": 0.35,
      "sample_size": 120
    },
    {
      "minimum_score": 0.5,
      "maximum_score": 1.0,
      "estimated_probability_correct": 0.82,
      "sample_size": 180
    }
  ]
}
```

The example numbers are illustrative and must not be used as a real calibration. Bins must be ordered and non-overlapping, probabilities must lie in `[0,1]`, and sample sizes must be positive. The result records the matching bin's source and sample size. Without a table, the JSON explicitly marks confidence as heuristic/internal. Warning codes make uncertainty actionable, including uncertain or sensitivity-unstable boundaries, influential copies, insufficient overlap, possible mixed subfamilies, cross-frame domains, missing coordinate evidence, domains on both strands, and possible host sequence.

## Benchmarking and calibration

`reputils benchmark` evaluates existing curation JSON files against curator-provided acceptable endpoint intervals:

```bash
reputils benchmark \
  --labels benchmark-labels.json \
  --bins 5 \
  --minimum-bin-size 10 \
  --out benchmark-results/
```

The input schema is `reputils.curation-benchmark.v1`:

```json
{
  "schema_version": "reputils.curation-benchmark.v1",
  "source": "description of the independent curator review",
  "families": [
    {
      "family_id": "family-1",
      "result_json": "curated-family/result.json",
      "accepted_left": {
        "minimum_1based": 91,
        "maximum_1based": 96
      },
      "accepted_right": {
        "minimum_1based": 1240,
        "maximum_1based": 1245
      },
      "accepted_architectures": ["tir-like"],
      "warning_labels": {
        "possible_mixed_subfamilies": false,
        "possible_host_gene": true
      }
    }
  ]
}
```

Relative result paths are resolved against the label file. Warning precision/recall considers only warning codes explicitly labeled true or false; omitted codes are unassessed. `benchmark.json` uses `reputils.curation-benchmark-result.v1` and records per-family endpoint errors, interval coverage, architecture accuracy, warning counts, result checksums, aggregate metrics, and calibration status. `families.tsv` provides the main metrics in tabular form.

When a score bin contains at least `--minimum-bin-size` labeled families, `calibration.json` is emitted using `reputils.confidence-calibration.v1`. Each generated bin records correct/total counts and a 95% Wilson interval as well as the empirical fraction. The table includes the label-file SHA-256 checksum. Generation does not prove transportability: validate the table on families independent of those used to fit it before passing it to `curate --confidence-calibration`.

## Known limitations

- Old, highly diverged elements may have too little homology for stable boundaries.
- Very few independent copies produce low or insufficient confidence.
- Mixed families and nested repeats can distort one-consensus analysis.
- Assembly errors and frameshifts can fragment structural/domain evidence.
- LINE 5′ truncation makes left boundaries intrinsically less certain.
- Non-autonomous elements may have strong termini but no coding domains.
- Coordinate-aware TSD calls assume correct zero-based, half-open unextended intervals and do not fetch/merge loci.
- Subfamily splitting is intentionally conservative; it does not recursively classify each emitted cluster and remains sensitive to sparse pairwise overlap.
- Cross-frame domain chaining identifies compatible evidence but does not infer an exact mutation, repair frameshifts, or replace specialist gene prediction.
- Host-domain evidence is a review signal, not full host-gene classification or an automatic rejection rule.
- Calibration-table quality is the user's responsibility; reputils does not ship a calibration data set or claim calibrated probabilities without one.
- Selected-copy protein searches are opt-in and bounded.

## Legacy commands

```text
con   consensus from an aligned FASTA
 div  sliding-window diversity
 dot  short-sequence self dotplots
html  legacy combined stdout HTML
 tir  quick consensus/reverse-complement alignment
 tsd  legacy matching substrings near alignment ends
 ttc  legacy fixed missingness/identity threshold trimmer
```

Run `reputils <command> --help` for its established options. The legacy `html` spelling bug for `--tsd_max_window` and unsafe `ttc` extension arithmetic were corrected; output formats otherwise remain unchanged.

## Development

Standard checks require no HMMER, MAFFT, or genome database:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all
```

External command tests should use mock executables and tiny fixture databases.
