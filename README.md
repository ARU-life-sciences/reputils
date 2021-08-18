# Repeat sequence utilities (reputils)

Help identify transposable elements (TE's) in manual curation by bringing together several lines of evidence. The input is mainly an alignment of a consensus TE sequence which has been blasted against the original genome, and aligned. These alignments are often fragmented and gappy.

The program of greatest use is probably `reputils html`. It brings together all of the functionality so far implemented into a single html document.

## Usage

Building <a href="https://www.rust-lang.org/tools/install">requires Rust</a>. 

```bash
git clone https://github.com/tolkit/reputils
cd reputils
cargo build --release
# ./target/release/reputils is the executable
# show help
./target/release/reputils --help
```

Everything is pretty much printed to stdout or to a PNG.

```
reputils 0.2.0
Max Brown <mb39@sanger.ac.uk>
reputils - some functions to aid TE identification.

USAGE:
    reputils [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    con     Make a consensus out of a multiple alignment fasta. Optimised for TE's.
    div     Calculate diversity along sliding windows of an alignment.
    dot     Make (self) dotplots from a fasta file. Suitable really only for short(ish) sequences.
    help    Prints this message or the help of the given subcommand(s)
    html    Render an HTML to gather several lines of identification evidence for a TE.
    tir     Take a consensus and quickly check for terminal inverted repeats (TIR)
    tsd     Try to find the Target Site Duplication of a TE. Prints a table.
    ttc     Trim an alignment to the core TE sequence.
```

### HTML overview

Combines all of the other functions so far into a single HTML document. See the `/examples` folder for an example. To see it in action, go <a href="https://htmlpreview.github.io/?https://github.com/tolkit/reputils/blob/main/examples/example.html"><b>here</b></a>.

```
reputils-html 
Render an HTML to gather several lines of identification evidence for a TE.

USAGE:
    reputils html --div_window_size <div_window_size> --div_window_step <div_window_step> --dot_nmatch <dot_nmatch> --dot_wsize <dot_wsize> --dot_wstep <dot_wstep> --fasta <fasta> --trim_extend <trim_extend> --trim_iden <trim_iden> --trim_miss <trim_miss> --trim_next_hit <trim_next_hit> --tsd_len <tsd_len> --tsd_min_window <tsd_max_window> --tsd_min_window <tsd_min_window>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --div_window_size <div_window_size>    The size of the window to iterate over. [default: 10]
        --div_window_step <div_window_step>    The step size of the window to iterate over. If equal to window, then
                                               windows are non-overlapping. [default: 3]
        --dot_nmatch <dot_nmatch>              Number of matches to tolerate a positive match. [default: 1]
        --dot_wsize <dot_wsize>                Window size to iterate over sequence. [default: 10]
        --dot_wstep <dot_wstep>                Window step size for window iterator. [default: 3]
    -f, --fasta <fasta>                        The multiple alignment file in fasta format.
        --trim_extend <trim_extend>            Extend alingment either end by number of bases specified. [default: 30]
        --trim_iden <trim_iden>                % identity in a column for the column to be considered a hit. [default:
                                               0.85]
        --trim_miss <trim_miss>                % missing data tolerated in a column. [default: 0.1]
        --trim_next_hit <trim_next_hit>        Isolated hits of well conserved columns leads to bad trimming. Play with
                                               this number? [default: 1]
        --tsd_len <tsd_len>                    Number of bases from beginning or end of alignment to query. [default:
                                               30]
        --tsd_min_window <tsd_max_window>      TSD's are searched for <= to this length. [default: 12]
        --tsd_min_window <tsd_min_window>      TSD's are searched for >= to this length. [default: 2]
```

### Make a consensus

```
reputils-con 
Make a consensus out of a multiple alignment fasta. Optimised for TE's.

USAGE:
    reputils con [FLAGS] [OPTIONS] --fasta <fasta>

FLAGS:
    -a, --append     Append the consensus to the input fasta.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>    The multiple alignment file in fasta format.
    -n, --name <name>      Name of the consensus sequence header. [default: CONS]
```

### Diversity in windows over a TE alignment

It's more clear when the alignment is trimmed to include only the putative TE.

```
reputils-div 
Calculate diversity along sliding windows of an alignment.

USAGE:
    reputils div [FLAGS] [OPTIONS] --fasta <fasta> --step <step> --window <window>

FLAGS:
    -h, --help       Prints help information
    -p, --plot       Plot the diversity across windows of the alignment. Output is a PNG.
    -V, --version    Prints version information

OPTIONS:
    -d, --dir <dir>          Directory to put plot in. [default: .]
    -f, --fasta <fasta>      The consensus sequence file in fasta format.
    -n, --name <name>        Name of the plot/PNG. [default: div_plot]
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
    reputils tir [FLAGS] --fasta <fasta>

FLAGS:
    -h, --help       Prints help information
    -s, --show       Pretty print the alignment.
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>    The consensus sequence file in fasta format.
```

### Trim alignment to core TE sequence

This script will take an alignment and trim it to the TE, plus any TSD's (hopefully). It needs a bit of testing, but worked on the Mariners I was looking at. TE's with 5' truncation may not work with this.

```
reputils-ttc 
Trim an alignment to the core TE sequence.

USAGE:
    reputils ttc --extend <extend> --fasta <fasta> --identity <identity> --missing <missing> --next_hit <next_hit>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -e, --extend <extend>        Extend the extracted alignment by `e` many bases either side of the alignment.
                                 [default: 15]
    -f, --fasta <fasta>          The multiple alignment sequence file in fasta format.
    -i, --identity <identity>    % identity in a column for the column to be considered a hit. [default: 0.8]
    -m, --missing <missing>      % missing data tolerated in a column. [default: 0.1]
    -n, --next_hit <next_hit>    Isolated hits of well conserved columns leads to bad trimming. Play with this number?
                                 [default: 1]
```

### Help to identify TSD's

Looks at the either end of a *trimmed* alignment (must be trimmed). I don't know how useful this actually is (it might confuse things more). But here it is:

```
reputils-tsd 
Try to find the Target Site Duplication of a TE. Prints a table.

USAGE:
    reputils tsd --fasta <fasta> --length <length> --maximum <maximum> --minimum <minimum>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>        The multiple alignment sequence file in fasta format.
    -l, --length <length>      Number of bases from beginning and end of alignment to query. [default: 20]
    -x, --maximum <maximum>    TSD's are searched for <= to this length. [default: 12]
    -m, --minimum <minimum>    TSD's are searched for >= to this length. [default: 2]

```

### Dotplot of sequences

Takes a fasta file and self compares each sequence.

```
reputils-dot 
Make (self) dotplots from a fasta file. Suitable really only for short(ish) sequences.

USAGE:
    reputils dot --dir <dir> --fasta <fasta> --nmatches <nmatches> --wsize <wsize> --wstep <wstep>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -d, --dir <dir>              Dirname where output plots should go. [default: dot]
    -f, --fasta <fasta>          The multiple alignment sequence file in fasta format.
    -n, --nmatches <nmatches>    Number of matches to tolerate a positive match. [default: 1]
    -i, --wsize <wsize>          Window size to iterate over sequence. [default: 10]
    -t, --wstep <wstep>          Window step size for window iterator. [default: 4]
```

<img src="examples/BDGG01000017.1_186586-190792.png">

### Performance

Performance will take a dip with large sequences (>10/100kb), especially `reputils dot`. If you want massive dotplots, there are many other more efficient programs out there!

### TODO's

- simple fasta stats of consensus sequences, sequence length, and length distribution (if multiple fastas)?
- Detect 5' truncation?