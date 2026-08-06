# Generation of this data

## Dotplot

```bash
# trim the alignment; a prediction from a tardigrade genome.
reputils ttc -f ramVar1-256#Unknown.gaps95.fa > ttc.fa
# make a dotplot
reputils dot -f ttc.fa --wsize 15 --wstep 1 --nmatches 5 --dir ../dot
```

## HTML

```bash
# an example html file
reputils html -f examples/ramVar1-256#Unknown.gaps95.fa > examples/example.html
```

## Domain rules

`domain-rules.json` demonstrates explicit `te-associated`, `host-enriched`, and `ubiquitous` rule categories. Patterns are illustrative; curate the mapping for the HMM database in use.

## Curation benchmark

`benchmark-labels.json` documents the `reputils.curation-benchmark.v1` structure. Its result path is a placeholder and must be replaced with a real curator-reviewed `result.json` before running:

```bash
reputils benchmark \
  --labels examples/benchmark-labels.json \
  --out benchmark-results/
```
