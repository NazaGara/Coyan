use clap::Parser;
use fault_tree::FaultTree;
use formula::CNFFormat;
use itertools::Itertools;
use rayon::prelude::*;
use serde_json::json;
use solver::*;
use std::path::Path;
use std::time::Instant;
use std::{fmt::Debug, str::FromStr};

mod fault_tree;
mod fault_tree_normalizer;
mod formula;
mod modularizer;
mod nodes;
mod solver;

#[derive(Parser, Debug)]
#[command(
    author = "Nazareno Garagiola",
    version = "0.2",
    about = "
        Coyan is a Rust project that transforms Static Fault Trees into CNF equisatisfiables formulas. 
        Then, with the use of a Model Counter, computes the Top Event Probabilty (TEP) of failure at one or more given time points.
    "
)]
struct Arguments {
    #[clap(subcommand)]
    command: Command,
}

#[derive(Parser, Clone, Debug)]
enum Command {
    #[clap(
        about = "Outputs information about the FT: the amount of basic events, of gates and the number of clauses generated by the method."
    )]
    Info(InfoCommand),
    #[clap(
        about = "Translates the FT implicit formula to a CNF equisatisfiable formula. Outputs a DIMACS file with a wcnf file."
    )]
    Translate(TranslateCommand),
    #[clap(
        about = "Executes a Solver to obtain the TEP of the FT at a given timepoint or timebound."
    )]
    Solve(SolveCommand),
    #[clap(
        about = "Modularize the input FT into all his modules, then compute the TEP of each module and replace the gate with a Basic Event, where the probability is the obtained TEP of the module."
    )]
    Modularize(ModCommand),
}

#[derive(Parser, Debug, Clone)]
struct InfoCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    input: String,
    /// Simplify the FT by removing one children gates.
    #[arg(long, default_value_t = false)]
    simplify: bool,
}
#[derive(Parser, Debug, Clone)]
struct SolveCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    input: String,
    /// Solver path and arguments.
    #[arg(short, long)]
    solver_path: String,
    /// Time bounds, creates a range of values according to the command arguments: [start, end, step].
    #[arg(
        long,
        value_delimiter = ' ',
        num_args = 3,
        conflicts_with = "timepoint"
    )]
    timebounds: Option<Vec<f64>>,
    /// Compute TEP of the FT a given timepoint. Conflicts with `timebounds`.
    #[arg(long, default_value_t = 1.0, conflicts_with = "timebounds")]
    timepoint: f64,
    /// Execution timeout for the WMC solver in seconds.
    #[arg(long, default_value_t = 300)]
    timeout_s: u64,
    /// Output format for the CNF formula. The format gives the extension to the file. Support values `MC21` and `MCC` [default: `MC21`]
    #[arg(long, default_value = "MC21")]
    format: String,
    /// Number of threads to use when time bounds are used.
    #[arg(long, default_value_t = 4)]
    num_threads: usize,
    /// Verbosity, if true, prints details at finish.
    #[arg(long, default_value_t = false)]
    verb: bool,
    /// Simplify the FT by removing one children gates.
    #[arg(long, default_value_t = true)]
    simplify: bool,
}

#[derive(Parser, Debug, Clone)]
struct TranslateCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    input: String,
    /// Output file, writes a .wcnf file.
    #[arg(short, long)]
    output: String,
    /// Output format for the CNF formula. The format gives the extension to the file. Support values `MC21` and `MCC` [default: `MC21`]
    #[arg(long, default_value = "MC21")]
    format: String,
    /// Time bounds, creates a range of values according to the command arguments: [start, end, step]. Conflicts with `timepoint`.
    #[arg(
        long,
        value_delimiter = ' ',
        num_args = 3,
        conflicts_with = "timepoint"
    )]
    timebounds: Option<Vec<f64>>,
    /// If provided, the weights will be written to a separate file.
    #[arg(long)]
    w_file: Option<String>,
    /// Compute TEP of the FT a given timepoint. Conflicts with `timebounds`.
    #[arg(long, default_value_t = 1.0, conflicts_with = "timebounds")]
    timepoint: f64,
    /// Simplify the FT by removing one children gates.
    #[arg(long, default_value_t = true)]
    simplify: bool,
}

#[derive(Parser, Debug, Clone)]
struct ModCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    input: String,
    /// Solver path and arguments.
    #[arg(short, long)]
    solver_path: String,
    // /// Compute TEP of the FT a given timepoint. Conflicts with `timebounds`.
    #[arg(long, default_value_t = 1.0)] //, conflicts_with = "timebounds")]
    timepoint: f64,
    /// Execution timeout for the WMC solver in seconds.
    #[arg(long, default_value_t = 300)]
    timeout_s: u64,
    // /// Output format for the CNF formula. The format gives the extension to the file. Support values `MC21` and `MCC` [default: `MC21`]
    #[arg(long, default_value = "MC21")]
    format: String,
    /// Number of threads to run in parallel to repalce BEs.
    #[arg(long, default_value_t = 1)]
    num_threads: usize,
    /// Verbosity, if true, prints details at finish.
    #[arg(long, default_value_t = true)]
    display: bool,
    /// Simplify the FT by removing one children gates.
    #[arg(long, default_value_t = true)]
    simplify: bool,
}

/// Outputs relevant information about the FT.
fn ft_info(command: InfoCommand) {
    let dft_filename = command.input;
    let simplify = command.simplify;
    // let mut ft: FaultTree<String> = FaultTree::new();
    let ft = FaultTree::new_from_file(&dft_filename, simplify);
    // ft.read_from_file(&dft_filename, simplify);
    let path = Path::new(dft_filename.as_str());
    let model_name = path.file_name().unwrap();
    let (num_be, num_gates, num_clauses) = ft.get_info();
    println!(
        "{}",
        json!({
            "model": model_name.to_str(),
            "num_basic_events": num_be,
            "num_gates": num_gates,
            "num_clauses": num_clauses,
        })
    );
}

/// Translates the FT explicit formula to a CNF file.
fn translate(command: TranslateCommand) {
    let dft_filename = command.input;
    let cnf_filename = command.output;
    let w_file = command.w_file;
    let simplify = command.simplify;
    let path = Path::new(dft_filename.as_str());
    let model_name = path.file_name().unwrap();
    let format =
        CNFFormat::from_str(&command.format).expect("Unsupported format. Try MCC or MC21.");

    // let mut ft: FaultTree<String> = FaultTree::new();
    let time_start = Instant::now();
    let ft = FaultTree::new_from_file(&dft_filename, simplify);
    // ft.read_from_file(&dft_filename, simplify);

    match command.timebounds {
        None => {
            let cnf_path = format!("{}_t={}.wcnf", cnf_filename, command.timepoint);
            ft.dump_cnf_to_file(cnf_path, format, command.timepoint, w_file);
            let duration = time_start.elapsed();
            println!(
                "{}",
                json!({
                    "model": model_name.to_str(),
                    "duration": format!("{:?}", duration),
                })
            );
        }
        Some(ts) => {
            let (start, end, step) = (ts[0], ts[1], ts[2]);
            let n_steps = if step != 0.0 {
                (end / step) as i64 + 1
            } else {
                2
            };
            let timebounds = (start as i64..n_steps)
                .into_iter()
                .map(|v| start + ((100 * v) as f64 * step).round() / 100.0)
                .collect_vec();
            for t in timebounds {
                let cnf_path = format!("{}_t={}.wcnf", cnf_filename, t);
                ft.dump_cnf_to_file(cnf_path, format, t.to_owned(), w_file.clone());
                let duration = time_start.elapsed();
                println!(
                    "{}",
                    json!({
                        "model": model_name.to_str(),
                        "duration": format!("{:?}", duration),
                    })
                );
            }
        }
    }
}

/// Compute TEP of FT, given a solver and the configuration needed.
/// Can perform multiple timepoints if range was given and use multi-threading
/// to handle each run.
fn compute_tep(command: SolveCommand) {
    let dft_filename = command.input;
    let solver_path = command.solver_path;
    let format =
        CNFFormat::from_str(&command.format).expect("Unsupported format. Try MCC or MC21.");

    rayon::ThreadPoolBuilder::new()
        .num_threads(command.num_threads)
        .build_global()
        .unwrap();
    let time_start = Instant::now();
    let ft = FaultTree::new_from_file(&dft_filename, command.simplify);
    match command.timebounds {
        None => {
            let solver: Box<dyn Solver> = get_solver_from_path(&solver_path);
            let tep = solver.compute_probabilty(&ft, format, command.timepoint, command.timeout_s);
            let duration = time_start.elapsed();
            if !command.verb {
                println!(
                    "{}",
                    json!({"TEP": tep,
                "timepoint": command.timepoint})
                )
            } else {
                let path = Path::new(dft_filename.as_str());
                let model_name = path.file_name().unwrap();
                println!(
                    "{}",
                    json!({
                        "model": model_name.to_str(),
                        "timepoint": command.timepoint,
                        "TEP": tep,
                        "duration": format!("{:?}", duration),
                    })
                );
            };
        }
        Some(ts) => {
            let (start, end, step) = (ts[0], ts[1], ts[2]);
            let n_steps = if step != 0.0 {
                (end / step) as i64 + 1
            } else {
                2
            };
            let timebounds = (start as i64..n_steps)
                .into_iter()
                .map(|v| start + ((100 * v) as f64 * step).round() / 100.0)
                .collect_vec();
            let _probs: Vec<(f64, f64)> = timebounds
                .into_par_iter()
                .filter_map(move |t| {
                    let ft = &ft;
                    if t > end {
                        None
                    } else {
                        let solver = get_solver_from_path(&solver_path);
                        let tep = solver.compute_probabilty(ft, format, t, command.timeout_s);
                        let duration = time_start.elapsed();
                        if !command.verb {
                            println!(
                                "{}",
                                json!({"TEP": tep,
                            "timepoint": t})
                            )
                        } else {
                            let path = Path::new(dft_filename.as_str());
                            let model_name = path.file_name().unwrap();
                            println!(
                                "{}",
                                json!({
                                    "model": model_name.to_str(),
                                    "timepoint": command.timepoint,
                                    "TEP": tep,
                                    "duration": format!("{:?}", duration),
                                })
                            );
                        };
                        Some((t, tep))
                    }
                })
                .collect();
        }
    }
}

fn modularize_ft(command: ModCommand) {
    let dft_filename = command.input;
    let format =
        CNFFormat::from_str(&command.format).expect("Unsupported format. Try MCC or MC21.");
    let solver_path = command.solver_path;
    let solver: Box<dyn Solver + Sync> = get_solver_from_path(&solver_path);
    let path = Path::new(dft_filename.as_str());
    let model_name = path.file_name().unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(command.num_threads)
        .build_global()
        .unwrap();

    let mut ft = FaultTree::new_from_file(&dft_filename, command.simplify);
    let time_start = Instant::now();

    let mut module_ids = ft.modularize_ft();
    let n_modules = module_ids.len();

    module_ids.reverse();

    ft.replace_modules(
        &solver,
        module_ids,
        format,
        command.timepoint,
        command.timeout_s,
        command.num_threads,
        command.display,
    );

    let tep = solver.compute_probabilty(&ft, format, command.timepoint, command.timeout_s);
    let elapsed = time_start.elapsed();
    println!(
        "{}",
        json!({
            "#modules" : n_modules,
            "model": model_name.to_str(),
            "timepoint": command.timepoint,
            "TEP": tep,
            "duration": format!("{:?}", elapsed),
        })
    );
}

fn main() {
    let arguments = Arguments::parse();
    match arguments.command {
        Command::Info(command) => ft_info(command),
        Command::Solve(command) => compute_tep(command),
        Command::Translate(command) => translate(command),
        Command::Modularize(command) => modularize_ft(command),
    }
}
