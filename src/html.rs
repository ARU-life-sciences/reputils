// module for making an html document

use bio::alignment::pairwise::*;
use bio::io::fasta;
use clap::value_t;
use itertools::izip;

use crate::con::get_consensus;
use crate::dot::match_case;
use crate::utils::alignment::{Alignment, Sequence};
use crate::utils::revcomp::reverse_complement;
use crate::utils::windows::SeqWindows;

pub fn render_html(matches: &clap::ArgMatches) {
    // parse command line args
    let fasta = matches.value_of("fasta").unwrap();
    let trim_extend = value_t!(matches.value_of("trim_extend"), usize).unwrap_or_else(|e| e.exit());
    let trim_next_hit =
        value_t!(matches.value_of("trim_next_hit"), usize).unwrap_or_else(|e| e.exit());
    let trim_miss = value_t!(matches.value_of("trim_miss"), f64).unwrap_or_else(|e| e.exit());
    let trim_iden = value_t!(matches.value_of("trim_iden"), f64).unwrap_or_else(|e| e.exit());
    let dot_wsize = value_t!(matches.value_of("dot_wsize"), usize).unwrap_or_else(|e| e.exit());
    let dot_wstep = value_t!(matches.value_of("dot_wstep"), usize).unwrap_or_else(|e| e.exit());
    let dot_nmatch = value_t!(matches.value_of("dot_nmatch"), usize).unwrap_or_else(|e| e.exit());
    let tsd_len = value_t!(matches.value_of("tsd_len"), usize).unwrap_or_else(|e| e.exit());
    let tsd_min_window =
        value_t!(matches.value_of("tsd_min_window"), usize).unwrap_or_else(|e| e.exit());
    let tsd_max_window =
        value_t!(matches.value_of("tsd_max_window"), usize).unwrap_or_else(|e| e.exit());
    let div_window_size =
        value_t!(matches.value_of("div_window_size"), usize).unwrap_or_else(|e| e.exit());
    let div_window_step =
        value_t!(matches.value_of("div_window_step"), usize).unwrap_or_else(|e| e.exit());

    let reader = fasta::Reader::from_file(fasta).expect("[-]\tPath invalid.");
    let mut read_number = 0i32;

    // read the fasta into Alignment
    // we need to do this three times
    // and make three copies of the alignment
    // large alignments will require large amounts of memory...
    let mut matrix = Alignment::new();
    let mut tsd_matrix = Alignment::new();
    let mut div_window_matrix = Alignment::new();

    for record in reader.records() {
        let fasta_record = record.expect("[-]\tError during fasta record parsing.");
        matrix.add_sequence(Sequence {
            name: fasta_record.id().to_string(),
            sequence: fasta_record.seq().to_vec(),
        });
        tsd_matrix.add_sequence(Sequence {
            name: fasta_record.id().to_string(),
            sequence: fasta_record.seq().to_vec(),
        });
        div_window_matrix.add_sequence(Sequence {
            name: fasta_record.id().to_string(),
            sequence: fasta_record.seq().to_vec(),
        });
        read_number += 1;
    }
    eprintln!("[+]\tAlignments in memory.");

    // find the blocks and trim
    let blocks = matrix.find_blocks(trim_miss, trim_iden);
    let trimmed = blocks.trim(matrix, trim_extend, trim_next_hit, true);

    // perhaps do something more sensible with this.
    let ok_trimmed = match trimmed {
        Some(t) => t,
        None => Alignment::new(),
    };
    eprintln!("[+]\tTrimmed alignment.");

    // for the header of the html
    let mut seq_names = String::new();
    seq_names += fasta;

    //
    // Making consensus:
    // Get the profile of an Alignment, then call get_consensus
    // Delete gaps.
    //
    let profile = ok_trimmed.get_profile();
    let mut consensus = get_consensus(profile, read_number);
    consensus.retain(|&e| e != 45);
    eprintln!("[+]\tConsensus sequence generated.");

    //
    // Make the dotplot:
    // SVG rendition of the dotplot from crate::dot::dot
    // Collection of <rect> 's made from iterating over the matrix
    //
    let seq_counter = SeqWindows::new(&consensus, dot_wsize, dot_wstep);
    let no_its = seq_counter.count();

    let seq_windows = SeqWindows::new(&consensus, dot_wsize, dot_wstep);

    let seq_windows_vec: Vec<&[u8]> = seq_windows.collect();

    // Base 1d array
    let mut matrix_raw = vec![false; no_its * no_its];
    // Vector of 'width' elements slices
    let mut matrix_base: Vec<_> = matrix_raw.as_mut_slice().chunks_mut(no_its).collect();
    // Final 2d array `&mut [&mut [_]]`
    let dot_matrix: &mut [&mut [bool]] = matrix_base.as_mut_slice();

    // compare windows to themselves
    for row in 0..no_its {
        for column in 0..no_its {
            // base or kmer 1
            let k1 = seq_windows_vec[row];
            // base or kmer 2
            let k2 = seq_windows_vec[column];
            dot_matrix[row][column] = match_case(k1, k2, dot_nmatch)
        }
    }
    eprintln!("[+]\tDotplot matrix generated.");

    // make svg of dotplot
    let mut dot_plot = String::new();

    let dot_svg_dim = dot_matrix.len();
    // add the initial tag
    dot_plot.push_str(&format!(
        r###"
                <svg viewBox="0 0 {} {}" width=80%
                    id="svg_dotplot"
                    xmlns="http://www.w3.org/2000/svg" >
                <rect width="100%" height="100%" fill="#F5F5DC" />
            "###,
        dot_svg_dim, dot_svg_dim
    ));

    // populate the svg matrix.
    for (row, _x) in dot_matrix.iter().enumerate() {
        for (col, val) in _x.iter().enumerate() {
            if *val {
                // do we want id's here for something later?
                dot_plot += &format!(
                    r###"
                            <rect x="{}" y ="{}" width="1" height="1" />
                        "###,
                    row, col
                )
            }
        }
    }
    dot_plot.push_str("</svg>");

    eprintln!("[+]\tSVG dotplot made.");

    //
    // Self align the consensus:
    // Consensus vs revcomp consensus
    // Print to an html table
    //
    eprintln!("[+]\tSelf aligning for TIR identification.");
    let reverse_consensus = reverse_complement(&consensus);

    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // gap open score: -5, gap extension score: -1
    let mut aligner =
        Aligner::with_capacity(consensus.len(), reverse_consensus.len(), -5, -1, &score);
    let alignment = aligner.semiglobal(&consensus, &reverse_consensus);
    let alignment_pretty = alignment.pretty(&consensus, &reverse_consensus);

    // merge every third entry in vec
    // kind of annoying really; probably be better to modify the actual rust bio api.
    let ali_print_parts: Vec<&str> = alignment_pretty.split('\n').collect();
    let ali_print_parts_f: Vec<&str> = ali_print_parts
        .iter()
        .filter(|e| !e.is_empty())
        .enumerate()
        .filter_map(|(i, e)| if i % 3 == 0 { Some(*e) } else { None })
        .collect();

    let ali_print_parts_m: Vec<&str> = ali_print_parts
        .iter()
        .filter(|e| !e.is_empty())
        .enumerate()
        .filter_map(|(i, e)| if i % 3 == 1 { Some(*e) } else { None })
        .collect();

    let ali_print_parts_r: Vec<&str> = ali_print_parts
        .iter()
        .filter(|e| !e.is_empty())
        .enumerate()
        .filter_map(|(i, e)| if i % 3 == 2 { Some(*e) } else { None })
        .collect();

    let ali_print_parts_f_str = ali_print_parts_f.join("");
    let ali_print_parts_m_str = ali_print_parts_m.join("");
    let ali_print_parts_r_str = ali_print_parts_r.join("");

    let mut fwd = String::new();
    let mut matches = String::new();
    let mut rev = String::new();

    // the three rows of the table
    fwd += "<tr>";
    matches += "<tr>";
    rev += "<tr>";

    for (f, _m, r) in izip!(
        ali_print_parts_f_str.chars(),
        ali_print_parts_m_str.chars(),
        ali_print_parts_r_str.chars()
    ) {
        // iterate over the chars and highlight matching
        // something like <p>AAGGC <span class="matching">TTCGGA</span>TAGGCATA...
        if f == r {
            fwd += r###"<td><span class="matching">"###;
            fwd += &f.to_string();
            fwd += r###"</span></td>"###;

            rev += r###"<td><span class="matching">"###;
            rev += &r.to_string();
            rev += r###"</span></td>"###;

            matches += "<td>|</td>";
        } else {
            fwd += r###"<td><span class="not-matching">"###;
            fwd += &f.to_string();
            fwd += r###"</span></td>"###;

            rev += r###"<td><span class="not-matching">"###;
            rev += &r.to_string();
            rev += r###"</span></td>"###;

            matches += "<td> </td>";
        }
    }

    fwd += "</tr>";
    matches += "</tr>";
    rev += "</tr>";

    let consensus_formatted = std::str::from_utf8(&consensus).unwrap().to_string();

    //
    // Target Site Duplication identification
    // Not entirely sure about this code, but it may help...
    // write a table of potential TSD's to an html table
    //
    let tsd = tsd_matrix
        .to_tsd_hash(tsd_len, tsd_min_window, tsd_max_window)
        .merge(true)
        .unwrap();

    // write this table to html
    let mut tsd_table = String::new();
    // start with <tr>'s as table can be put below
    let max_tsd_lens = tsd
        .iter()
        .max_by(|x, y| x.1.len().cmp(&y.1.len()))
        .unwrap()
        .1
        .len();

    for (id, tsds) in &tsd {
        // row will be id
        // column number dictated by tsds.len()
        // add a row

        tsd_table += "<tr>";
        tsd_table += &format!("<td>{}</td>", id);

        let mut tsd_string = String::new();

        for i in 0..max_tsd_lens {
            if i == 0 {
                match tsds.get(i) {
                    Some(x) => tsd_string += &x.to_ascii_uppercase().to_string(),
                    None => tsd_string += "<i>None detected</i>",
                }
            } else {
                match tsds.get(i) {
                    Some(x) => tsd_string += &format!(", {}", x.to_ascii_uppercase()),
                    None => tsd_string += "",
                }
            }
        }

        tsd_table += &format!("<td>{}</td>", tsd_string);

        tsd_table += "</tr>";
    }

    //
    // Windows of diversity across TE
    // Much functionality poached from tidk::plot (https://github.com/tolkit/telomeric-identifier/blob/main/src/plot.rs)
    // With some extra x & y axis labels
    //
    eprintln!("[+]\tMaking diversity windows.");

    let div_window_data = div_window_matrix.div_windows(div_window_size, div_window_step, true);

    // get the maximum y value
    let div_y_max = div_window_data
        .iter()
        .map(|(_x, _y, z)| z)
        .max_by(|x, y| x.partial_cmp(y).expect("Tried to compare a NaN"))
        .unwrap();

    // make svg of dotplot
    let mut div_plot = String::new();
    let div_height = 700;

    // width is the length of the sequence
    // but this is scaled by width = 90%
    let div_svg_width = div_window_data.len();
    // add the initial tag
    div_plot.push_str(&format!(
        r###"
                <svg viewBox="0 0 {} {}" width=90%
                    id="svg_divplot"
                    xmlns="http://www.w3.org/2000/svg" >
                <rect width="100%" height="100%" fill="#F5F5DC" />
            "###,
        div_svg_width, div_height
    ));

    // scale a range [min, max] to custom range [a, b]
    // our range will be [0, height of plot] == 300
    fn scale_y(y: f64, a: f64, b: f64, min: f64, max: f64) -> f64 {
        (((b - a) * (y - min)) / (max - min)) + a
    }

    // make the path element
    let mut path = String::new();
    const MARGIN: usize = 150;
    // margin at either side of the plot
    let div_width_incl_margins = div_svg_width - (MARGIN);

    // distance between points on the x axis
    let div_x_bin: f64 = div_width_incl_margins as f64 / div_svg_width as f64;

    // move to the initial point
    path += &format!(
        "M{},{}",
        MARGIN as f64 / 2f64,
        div_height as f64
            - MARGIN as f64 / 2.0
            - scale_y(
                div_window_data[0].2 as f64,
                0.0,
                div_height as f64 - MARGIN as f64 / 2.0,
                0.0,
                *div_y_max as f64
            )
    );

    // keep track of x bins
    let mut div_bin: f64 = div_x_bin;
    // hold the <text> nodes
    let mut div_x_axis_text = String::new();
    // 5 labels on the x axis.
    let div_x_axis_divisions = div_svg_width / 5;
    let mut division_counter = 0;

    // make divisions
    let mut test_vec = Vec::new();
    for _ in 0..5 {
        test_vec.push(division_counter);
        division_counter += div_x_axis_divisions;
    }

    // iterate over the data and populate the path element
    for (index, element) in div_window_data.iter().skip(1).enumerate() {
        // make the path element
        path += &format!(
            "L{},{}",
            div_bin + (MARGIN as f64 / 2.0),
            div_height as f64
                - MARGIN as f64 / 2.0
                - scale_y(
                    element.2 as f64,
                    0.0,
                    div_height as f64 - MARGIN as f64 / 2.0,
                    0.0,
                    *div_y_max as f64
                )
        );
        // create the axis text elements
        if test_vec.contains(&index) {
            div_x_axis_text += &format!(
                r###"<text x="{x}" y="{y}" class="x_axis_text">{text}</text>"###,
                x = div_bin + MARGIN as f64 / 2.0,
                y = div_height - (MARGIN / 4),
                text = format!("{} bp", element.0)
            );
        }

        // increment bins
        div_bin += div_x_bin;
    }

    div_plot += &format!(
        r###"<path d="{}" id="div_path" stroke="black" fill="none" stroke-width="1" />"###,
        path
    );
    // add x axis text
    div_plot += &div_x_axis_text;
    // add y axis text
    let div_y_axis_divisions = div_height / 5;
    let diversity_divisions = div_y_max / 5.0;
    let mut mut_diversity_divisions = div_height as f32;
    let mut diversity_text = 0f32;
    let mut diversity_text_string = String::new();

    for _ in 0..5 {
        if diversity_text != 0.0 {
            diversity_text_string += &format!(
                r###"<text x="{}" y="{}" class="y_axis_text">{:.2}</text>"###,
                0, mut_diversity_divisions, diversity_text
            );
        }

        diversity_text += diversity_divisions;
        mut_diversity_divisions -= div_y_axis_divisions as f32;
    }

    div_plot += &diversity_text_string;

    div_plot += "</svg>";

    //
    // The final html:
    // Add all the strings created above
    // Some minimal styling; print to stdout.
    // It's not properly formatted, but it works.
    //

    let html = format!(
        r###"
            <!DOCTYPE html>
            <html>
                <head>
                    <meta charset="utf#-8">
                    <title>reputils html</title>
                </head>
                <style>
                    #consensus_sequence {{
                        word-wrap: break-word; 
                        max-width: 80ch;
                    }}
                    .matching {{
                        color: red; 
                    }}
                    table {{
                        display: block;
                        overflow-x: auto;
                        white-space: nowrap;
                        border: 3px solid black;
                        border-collapse: collapse;
                    }}
                    .tsds td {{
                        border-left:1px solid #000;
                        padding: 0 15px;
                        padding-bottom: 1em;
                    }}
                    tr th {{
                        border-bottom:2px solid #000;
                    }}
                    .x_axis_text {{
                        font-size: 35px;
                        font-family: sans-serif;
                    }}
                    .y_axis_text {{
                        font-size: 35px;
                        font-family: sans-serif;
                    }}
                </style>
                <body>
                    <h2>{title}</h2>
                    <h3>Consensus sequence generated:</h3>
                    <p>>{fasta_header}</p>
                    <p id="consensus_sequence">{consensus}</p>
                    <h3>Dotplot of consensus sequence:</h3>
                    <p id="parameters">
                        Window size - {window_size}; step size - {step_size}; number of mismatches tolerated - {mismatches}
                    </p>
                    {dotplot}
                    <h3>TIR:</h3>
                    {terminal_inverted_repeat}
                    <h3>TSD's:</h3>
                    <p>A table of potential target site duplications. K-mers shown in the table are present at either ends of the sequence.</p>
                    {target_site_duplication_table}
                    <h3>Nuleotide diversity in windows:</h3>
                    <p>X-axis indicates number of base pairs into TE and the y-axis shows the nucleotide diversity (pi) at that point in the alignment.</p>
                    <p>Note, gaps present in the alignment will be depicted here.</p>
                    <p>Window size - {diversity_window_size}; step size - {diversity_step_size}</p>
                    {diversity_windows_plot}
                </body>
            </html>
            "###,
        title = seq_names,
        fasta_header = seq_names,
        consensus = consensus_formatted,
        window_size = dot_wsize,
        step_size = dot_wstep,
        mismatches = dot_nmatch,
        dotplot = dot_plot,
        terminal_inverted_repeat = format!(
            "<p>Forward consensus aligned to reverse complement of consensus:</p>
                <p>Scroll left/right to view more of alignment.</p>
                <table>
                {}{}{}
                </table>",
            fwd, matches, rev
        ),
        target_site_duplication_table = format!(
            r###"<table class="tsds">
                    <tr id="first_row">
                        <th><b>Alignment ID: <b></th>
                        <th><b>TSD's &rarr;<b></th>
                    </tr>
                    {}
                    </table>"###,
            tsd_table
        ),
        diversity_windows_plot = div_plot,
        diversity_window_size = div_window_size,
        diversity_step_size = div_window_step,
    );
    println!("{}", html);
}
