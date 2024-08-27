use clap::Parser;
use coyan_fta::{formula::CNFFormat, solver::get_solver_from_path};
use rand::Rng;
use random_fault_trees::{RFTConfig, RFaultTree};
use serde_json::json;
use std::{fmt::Debug, str::FromStr, time::Instant};
mod random_fault_trees;

/// CMD Arguments
#[derive(Parser, Debug)]
#[command(
    author = "Nazareno Garagiola",
    version = "0.2",
    about = "
        Module to create Random Fault Trees using discrete probabilities.
        - |Basic Events| = n_nodes * rate_be.
        - |Gates| = n_nodes - |Basic Events|
        - Requires: Sum(Gate rates) = 1. 
    "
)]
struct Args {
    /// Total number of nodes.
    #[arg(short, long)]
    n_nodes: usize,
    /// Output file, writes .dft files.
    #[arg(short, long)]
    output: String,
    /// Number in (0,1]. Rate of Basic Events from the total amount of Nodes.
    #[arg(long, default_value_t = 0.5, requires = "rate_and")]
    rate_be: f64,
    /// Number in (0,1]. Rate of AND gates from the total amount of gates.
    #[arg(long, default_value_t = 0.5, requires = "rate_or")]
    rate_and: f64,
    /// Number in (0,1]. Rate of OR gates from the total amount of gates.
    #[arg(long, default_value_t = 0.5)]
    rate_or: f64,
    /// Number in (0,1]. Rate of VOT gates from the total amount of gates.
    /// If the FT contains a VOT gate, cannot be processed and solved, Can be solved but after writing to dft.
    #[arg(long, default_value_t = 0.0, conflicts_with = "solver_path")]
    rate_vot: f64,
    ///Number in (0,1], multiplies the random generated float of the probability, so it can be a smaller probability.
    #[arg(long, default_value_t = 1e-4)]
    prob_multiplier: f64,
    /// In which percentage of the last gates start to put the Basic Events if they were not used before.
    #[arg(long, default_value_t = 0.6)]
    perc_last: f64,
    /// Specify the max number of children that a gate can have.
    #[arg(long, default_value_t = 5)]
    max_n_children: usize,
    /// Execution timeout for the WMC solver in seconds.
    #[arg(long, default_value_t = 100)]
    timeout_s: u64,
    /// In which percentage of the last gates start to put the Basic Events if they were not used before. [Default=random]
    #[arg(long, value_parser = clap::value_parser!(u64))]
    seed: Option<u64>,
    /// Solver path and arguments.
    /// First is the solvers path, then the prefix for the args and then the arguments
    #[arg(short, long, conflicts_with = "rate_vot")]
    solver_path: Option<String>,
    /// Output format for the CNF formual. The format gives the extension to the file. Currently supports MC21 and MCC.
    #[arg(long, default_value = "MC21")]
    format: Option<String>,
}

fn main() {
    let args = Args::parse();

    let n_nodes = args.n_nodes;
    let rates = vec![args.rate_be, args.rate_and, args.rate_or, args.rate_vot];
    let output_filename = args.output;
    let seed = args.seed;

    let seed = match seed {
        None => {
            let mut rng = rand::thread_rng();
            rng.gen_range(u16::MIN..u16::MAX) as u64
        }
        Some(u) => u,
    };

    let format = match args.format {
        Some(f) => CNFFormat::from_str(&f).expect("Unsupported format. Try MCC or MC21."),
        None => CNFFormat::MC21,
    };

    let output_filename = if !output_filename.ends_with(".dft") {
        format!("{}.dft", output_filename)
    } else {
        output_filename
    };

    let config = RFTConfig::from_vec(rates);
    let solver_cmd = &args.solver_path;

    let start = Instant::now();
    let rft = RFaultTree::new_random(
        n_nodes,
        config,
        args.prob_multiplier,
        args.perc_last,
        seed,
        args.max_n_children,
    );

    match solver_cmd {
        Option::None => {
            rft.save_to_dft(output_filename);
            let duration = start.elapsed();
            println!(
                "{}",
                json!({
                    "time_elapsed": format!("{:?}", duration),
                })
            );
        }
        Option::Some(cmd) => {
            let solver = get_solver_from_path(&cmd);
            rft.save_to_dft(output_filename);
            let ft = rft.extract_ft();
            let wmc = solver.compute_probabilty(&ft, format, 1.0, args.timeout_s, None, false);
            let duration = start.elapsed();

            println!(
                "{}",
                json!({
                    "solver": solver._name(),
                    "tep": wmc,
                    "time_elapsed": format!("{:?}", duration),
                })
            );
        }
    }
}
