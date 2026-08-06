use clap::{App, Arg};
use std::process;

use reputils::con::make_consensus;
use reputils::curate::{run as run_curate, CurateOptions};
use reputils::div::diversity_windows;
use reputils::dot::dot;
use reputils::html::render_html;
use reputils::tir::revcomp_alignment;
use reputils::tsd::find_tsds;
use reputils::ttc::ttc;
use std::path::PathBuf;

fn main() {
    let matches = App::new("reputils")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("reputils - some functions to aid TE identification.")
        .subcommand(
            clap::SubCommand::with_name("html")
                .about("Render an HTML to gather several lines of identification evidence for a TE.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The multiple alignment file in fasta format."),
                )
                .arg(
                    Arg::with_name("trim_extend")
                        .long("trim_extend")
                        .takes_value(true)
                        .required(true)
                        .default_value("30")
                        .help("Extend alingment either end by number of bases specified."),
                )
                .arg(
                    Arg::with_name("trim_next_hit")
                        .long("trim_next_hit")
                        .takes_value(true)
                        .required(true)
                        .default_value("1")
                        .help("Isolated hits of well conserved columns leads to bad trimming. Play with this number?"),
                )
                .arg(
                    Arg::with_name("trim_miss")
                        .long("trim_miss")
                        .takes_value(true)
                        .required(true)
                        .default_value("0.1")
                        .help("% missing data tolerated in a column."),
                )
                .arg(
                    Arg::with_name("trim_iden")
                        .long("trim_iden")
                        .takes_value(true)
                        .required(true)
                        .default_value("0.85")
                        .help("% identity in a column for the column to be considered a hit."),
                )
                .arg(
                    Arg::with_name("dot_wsize")
                        .long("dot_wsize")
                        .takes_value(true)
                        .required(true)
                        .default_value("10")
                        .help("Window size to iterate over sequence."),
                )
                .arg(
                    Arg::with_name("dot_wstep")
                        .long("dot_wstep")
                        .takes_value(true)
                        .required(true)
                        .default_value("3")
                        .help("Window step size for window iterator."),
                )
                .arg(
                    Arg::with_name("dot_nmatch")
                        .long("dot_nmatch")
                        .takes_value(true)
                        .required(true)
                        .default_value("1")
                        .help("Number of matches to tolerate a positive match."),
                )
                .arg(
                    Arg::with_name("tsd_len")
                        .long("tsd_len")
                        .takes_value(true)
                        .required(true)
                        .default_value("30")
                        .help("Number of bases from beginning or end of alignment to query."),
                )
                .arg(
                    Arg::with_name("tsd_min_window")
                        .long("tsd_min_window")
                        .takes_value(true)
                        .required(true)
                        .default_value("2")
                        .help("TSD's are searched for >= to this length."),
                )
                .arg(
                    Arg::with_name("tsd_max_window")
                        .long("tsd_max_window")
                        .takes_value(true)
                        .required(true)
                        .default_value("12")
                        .help("TSD's are searched for <= to this length."),
                )
                .arg(
                    Arg::with_name("div_window_size")
                        .long("div_window_size")
                        .takes_value(true)
                        .required(true)
                        .default_value("10")
                        .help("The size of the window to iterate over."),
                )
                .arg(
                    Arg::with_name("div_window_step")
                        .long("div_window_step")
                        .takes_value(true)
                        .required(true)
                        .default_value("3")
                        .help("The step size of the window to iterate over. If equal to window, then windows are non-overlapping."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("con")
                .about("Make a consensus out of a multiple alignment fasta. Optimised for TE's.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The multiple alignment file in fasta format."),
                )
                .arg(
                    Arg::with_name("name")
                        .short("n")
                        .long("name")
                        .takes_value(true)
                        .required(false)
                        .default_value("CONS")
                        .help("Name of the consensus sequence header."),
                )
                .arg(
                    Arg::with_name("append")
                        .short("a")
                        .long("append")
                        .help("Append the consensus to the input fasta."),
                ),
        )
        .subcommand(
            clap::SubCommand::with_name("tir")
                .about("Take a consensus and quickly check for terminal inverted repeats (TIR)")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The consensus sequence file in fasta format."),
                )
                .arg(
                    Arg::with_name("show")
                        .short("s")
                        .long("show")
                        .help("Pretty print the alignment."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("div")
                .about("Calculate diversity along sliding windows of an alignment.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The consensus sequence file in fasta format."),
                )
                .arg(
                    Arg::with_name("window")
                        .short("w")
                        .long("window")
                        .takes_value(true)
                        .required(true)
                        .default_value("25")
                        .help("The size of the window to iterate over."),
                )
                .arg(
                    Arg::with_name("step")
                        .short("s")
                        .long("step")
                        .takes_value(true)
                        .required(true)
                        .default_value("25")
                        .help("The step size of the window to iterate over. If equal to window, then windows are non-overlapping."),
                )
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .takes_value(true)
                        .default_value(".")
                        .help("Directory to put plot in."),
                )
                .arg(
                    Arg::with_name("name")
                        .short("n")
                        .long("name")
                        .takes_value(true)
                        .default_value("div_plot")
                        .help("Name of the plot/PNG."),
                )
                .arg(
                    Arg::with_name("plot")
                        .short("p")
                        .long("plot")
                        .help("Plot the diversity across windows of the alignment. Output is a PNG."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("ttc")
                .about("Trim an alignment to the core TE sequence.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The multiple alignment sequence file in fasta format."),
                )
                .arg(
                    Arg::with_name("extend")
                        .short("e")
                        .long("extend")
                        .takes_value(true)
                        .required(true)
                        .default_value("15")
                        .help("Extend the extracted alignment by `e` many bases either side of the alignment."),
                )
                .arg(
                    Arg::with_name("next_hit")
                        .short("n")
                        .long("next_hit")
                        .takes_value(true)
                        .required(true)
                        .default_value("1")
                        .help("Isolated hits of well conserved columns leads to bad trimming. Play with this number?"),
                )
                .arg(
                    Arg::with_name("missing")
                        .short("m")
                        .long("missing")
                        .takes_value(true)
                        .required(true)
                        .default_value("0.1")
                        .help("% missing data tolerated in a column."),
                )
                .arg(
                    Arg::with_name("identity")
                        .short("i")
                        .long("identity")
                        .takes_value(true)
                        .required(true)
                        .default_value("0.8")
                        .help("% identity in a column for the column to be considered a hit."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("dot")
                .about("Make (self) dotplots from a fasta file. Suitable really only for short(ish) sequences.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The multiple alignment sequence file in fasta format."),
                )
                .arg(
                    Arg::with_name("wsize")
                        .short("i")
                        .long("wsize")
                        .takes_value(true)
                        .required(true)
                        .default_value("10")
                        .help("Window size to iterate over sequence."),
                )
                .arg(
                    Arg::with_name("wstep")
                        .short("t")
                        .long("wstep")
                        .takes_value(true)
                        .required(true)
                        .default_value("4")
                        .help("Window step size for window iterator."),
                )
                .arg(
                    Arg::with_name("nmatches")
                        .short("n")
                        .long("nmatches")
                        .takes_value(true)
                        .required(true)
                        .default_value("1")
                        .help("Number of matches to tolerate a positive match."),
                )
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .takes_value(true)
                        .required(true)
                        .default_value("dot")
                        .help("Dirname where output plots should go."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("tsd")
                .about("Try to find the Target Site Duplication of a TE. Prints a table.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The multiple alignment sequence file in fasta format."),
                )
                .arg(
                    Arg::with_name("length")
                        .short("l")
                        .long("length")
                        .takes_value(true)
                        .required(true)
                        .default_value("20")
                        .help("Number of bases from beginning or end of alignment to query."),
                )
                .arg(
                    Arg::with_name("minimum")
                        .short("m")
                        .long("minimum")
                        .takes_value(true)
                        .required(true)
                        .default_value("2")
                        .help("TSD's are searched for >= to this length."),
                )
                .arg(
                    Arg::with_name("maximum")
                        .short("x")
                        .long("maximum")
                        .takes_value(true)
                        .required(true)
                        .default_value("12")
                        .help("TSD's are searched for <= to this length."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("curate")
                .about("Propose evidence-based TE family boundaries, consensus, and curation report.")
                .arg(Arg::with_name("alignment").long("alignment").takes_value(true).required(true).help("Aligned extended genomic copies in FASTA format."))
                .arg(Arg::with_name("out").long("out").takes_value(true).required(true).help("Output directory."))
                .arg(Arg::with_name("manifest").long("manifest").takes_value(true).help("Optional rep.copy-recovery.v1 JSON manifest."))
                .arg(Arg::with_name("query").long("query").takes_value(true).help("Optional provisional consensus FASTA."))
                .arg(Arg::with_name("genome").long("genome").takes_value(true).help("Optional source genome FASTA."))
                .arg(Arg::with_name("coordinates").long("coordinates").takes_value(true).help("Optional per-copy coordinate table."))
                .arg(Arg::with_name("hmm_db").long("hmm-db").takes_value(true).help("Profile-HMM database for hmmscan."))
                .arg(Arg::with_name("repeatpeps").long("repeatpeps").takes_value(true).help("Alternative protein-domain HMM database."))
                .arg(Arg::with_name("threads").long("threads").takes_value(true).default_value("1"))
                .arg(Arg::with_name("hmmscan").long("hmmscan").takes_value(true).help("hmmscan executable path; also REPUTILS_HMMSCAN."))
                .arg(Arg::with_name("hmmsearch").long("hmmsearch").takes_value(true).help("hmmsearch executable path; also REPUTILS_HMMSEARCH."))
                .arg(Arg::with_name("getorf").long("getorf").takes_value(true).help("getorf executable path; also REPUTILS_GETORF."))
                .arg(Arg::with_name("domain_rules").long("domain-rules").takes_value(true).help("Optional JSON domain interpretation rule file."))
                .arg(Arg::with_name("boundary_override").long("boundary-override").takes_value(true).help("Apply a reputils.curation-override.v1 JSON boundary override."))
                .arg(Arg::with_name("confidence_calibration").long("confidence-calibration").takes_value(true).help("Optional labeled reputils.confidence-calibration.v1 JSON table."))
                .arg(Arg::with_name("scan_selected_copies").long("scan-selected-copies").help("Also scan selected genomic copies for protein domains."))
                .arg(Arg::with_name("split_subfamilies").long("split-subfamilies").help("Emit additional consensuses only for strongly supported candidate subfamilies."))
                .arg(Arg::with_name("max_copy_domain_targets").long("max-copy-domain-targets").takes_value(true).default_value("10").help("Maximum selected copies searched for domains."))
                .arg(Arg::with_name("min_copies").long("min-copies").takes_value(true).default_value("3"))
                .arg(Arg::with_name("force").long("force").help("Replace known curation output files in an existing directory."))
        )
        .get_matches();

    let subcommand = matches.subcommand();
    match subcommand.0 {
        "curate" => {
            let matches = subcommand.1.unwrap();
            let path = |name: &str| matches.value_of(name).map(PathBuf::from);
            let options = CurateOptions {
                alignment: PathBuf::from(matches.value_of("alignment").unwrap()),
                out: PathBuf::from(matches.value_of("out").unwrap()),
                manifest: path("manifest"),
                query: path("query"),
                genome: path("genome"),
                coordinates: path("coordinates"),
                hmm_db: path("hmm_db"),
                repeatpeps: path("repeatpeps"),
                threads: matches
                    .value_of("threads")
                    .unwrap()
                    .parse()
                    .unwrap_or_else(|_| {
                        eprintln!("error: --threads must be a positive integer");
                        process::exit(2)
                    }),
                hmmscan: path("hmmscan"),
                hmmsearch: path("hmmsearch"),
                getorf: path("getorf"),
                domain_rules: path("domain_rules"),
                boundary_override: path("boundary_override"),
                confidence_calibration: path("confidence_calibration"),
                scan_selected_copies: matches.is_present("scan_selected_copies"),
                split_subfamilies: matches.is_present("split_subfamilies"),
                max_copy_domain_targets: matches
                    .value_of("max_copy_domain_targets")
                    .unwrap()
                    .parse()
                    .unwrap_or_else(|_| {
                        eprintln!("error: --max-copy-domain-targets must be a positive integer");
                        process::exit(2)
                    }),
                min_copies: matches
                    .value_of("min_copies")
                    .unwrap()
                    .parse()
                    .unwrap_or_else(|_| {
                        eprintln!("error: --min-copies must be a positive integer");
                        process::exit(2)
                    }),
                force: matches.is_present("force"),
            };
            match run_curate(&options) {
                Ok(result) => eprintln!(
                    "[+] Curation complete: {} ({:?} confidence)",
                    options.out.display(),
                    result.confidence_tier
                ),
                Err(error) => {
                    eprintln!("error: {}", error);
                    process::exit(1);
                }
            }
        }
        "con" => {
            let matches = subcommand.1.unwrap();
            make_consensus(matches);
        }
        "tir" => {
            let matches = subcommand.1.unwrap();
            revcomp_alignment(matches);
        }
        "div" => {
            let matches = subcommand.1.unwrap();
            diversity_windows(matches);
        }
        "tsd" => {
            let matches = subcommand.1.unwrap();
            find_tsds(matches);
        }
        "ttc" => {
            let matches = subcommand.1.unwrap();
            ttc(matches);
        }
        "dot" => {
            let matches = subcommand.1.unwrap();
            dot(matches).unwrap();
        }
        "html" => {
            let matches = subcommand.1.unwrap();
            render_html(matches);
        }
        _ => {
            eprintln!(
                "[-]\tRun with a subcommand, run with '--help' or '-h' for options. Exiting."
            );
            process::exit(1);
        }
    }
}
