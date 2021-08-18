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