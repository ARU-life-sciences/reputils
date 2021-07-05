use clap::{App, Arg};
use std::process;

// TODO list
//  : alignment stats? Is that useful?
use reputils::con::con::make_consensus;
use reputils::div::div::diversity_windows;
use reputils::tir::tir::revcomp_alignment;
// use reputils::tsd::tsd::find_tsds;
// use reputils::ttc::ttc::ttc;

fn main() {
    let matches = App::new("reputils")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("reputils - some functions to aid TE discovery.")
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
                        .takes_value(true)
                        .required(false)
                        .default_value("true")
                        .possible_values(&["true", "false"])
                        .help("Append the consensus to the input fasta? Otherwise print consensus alone."),
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
                        .takes_value(true)
                        .required(true)
                        .default_value("false")
                        .possible_values(&["true", "false"])
                        .help("Pretty print the alignment?"),
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
                    Arg::with_name("plot")
                        .short("p")
                        .long("plot")
                        .takes_value(true)
                        .required(true)
                        .default_value("false")
                        .possible_values(&["true", "false"])
                        .help("Plot the diversity across windows of the alignment? Output is a PNG."),
                )
        )
        // .subcommand(
        //     clap::SubCommand::with_name("tsd")
        //         .about("Determine possible target site duplications")
        //         .arg(
        //             Arg::with_name("fasta")
        //                 .short("f")
        //                 .long("fasta")
        //                 .takes_value(true)
        //                 .required(true)
        //                 .help("The consensus sequence file in fasta format."),
        //         )
        // )
        .get_matches();

    let subcommand = matches.subcommand();
    match subcommand.0 {
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
        // "tsd" => {
        //     let matches = subcommand.1.unwrap();
        //     find_tsds(matches);
        // }
        // "ttc" => {
        //     let matches = subcommand.1.unwrap();
        //     ttc(matches);
        // }
        _ => {
            eprintln!(
                "[-]\tSubcommand invalid, run with '--help' or '-h' for subcommand options. Exiting."
            );
            process::exit(1);
        }
    }
}