# Repeat sequence utilities

Inspired by the TE workshop recently (run by Alex Suh). Designed to be used on files downstream of RepeatModeler2/RepeatMasker2. Might have some utility beyond this.

## Usage

Everything is pretty much printed to stdout.

```
reputils 0.1.0
Max Brown <mb39@sanger.ac.uk>
reputils - some functions to aid TE discovery.

USAGE:
    reputils [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    con     Make a consensus out of a multiple alignment fasta. Optimised for TE's.
    div     Calculate diversity along sliding windows of an alignment.
    help    Prints this message or the help of the given subcommand(s)
    tir     Take a consensus and quickly check for terminal inverted repeats (TIR)
```

### Make a consensus

```
reputils-con 
Make a consensus out of a multiple alignment fasta. Optimised for TE's.

USAGE:
    reputils con [OPTIONS] --fasta <fasta>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -a, --append <append>    Append the consensus to the input fasta? Otherwise print consensus alone. [default: true]
                             [possible values: true, false]
    -f, --fasta <fasta>      The multiple alignment file in fasta format.
    -n, --name <name>        Name of the consensus sequence header. [default: CONS]
```

### Diversity in windows over a TE alignment

It's more clear when the alignment is trimmed to include only the putative TE.

```
reputils-div 
Calculate diversity along sliding windows of an alignment.

USAGE:
    reputils div --fasta <fasta> --plot <plot> --step <step> --window <window>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>      The consensus sequence file in fasta format.
    -p, --plot <plot>        Plot the diversity across windows of the alignment? Output is a PNG. [default: false]
                             [possible values: true, false]
    -s, --step <step>        The step size of the window to iterate over. If equal to window, then windows are non-
                             overlapping. [default: 25]
    -w, --window <window>    The size of the window to iterate over. [default: 25]
```

### Presence of TIR's

Quickly check whether your consensus sequence has TIR's.

```
reputils-tir 
Take a consensus and quickly check for terminal inverted repeats (TIR)

USAGE:
    reputils tir --fasta <fasta> --show <show>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>    The consensus sequence file in fasta format.
    -s, --show <show>      Pretty print the alignment? [default: false]  [possible values: true, false]
```

TODO list:
- Check beginning & end of sequences for motif (probably hard to do and not worth it.)
- Windows of diversity along alignment. Done.
- Revcomp sequence & self align? Done
- simple fasta stats of consensus sequences, sequence length, and length distribution (if multiple fastas).
- Make consensus of a fasta. Done
- Can we detect TSD's????