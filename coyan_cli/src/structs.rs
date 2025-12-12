use clap::Parser;

#[derive(Parser, Debug, Clone)]
pub struct InfoCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    pub input: String,
    /// Simplify the FT by removing one children gates. [default: false]
    #[arg(short, long, default_value_t = false)]
    pub simplify: bool,
    /// Get the number of sub-modules of the FT. [default: false]
    #[arg(short, long, default_value_t = false)]
    pub modularize: bool,
    /// If provided, postprocess the CNF formula by passing a CNF preprocessor. [default: None]
    #[arg(short, long, default_value = None)]
    pub preprocess: Option<String>,
}
#[derive(Parser, Debug, Clone)]
pub struct SolveCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    pub input: String,
    /// Solver path and arguments.
    #[arg(short, long)]
    pub solver_path: String,
    /// Compute TEP at a specific timepoint
    #[arg(short, long, default_value_t = 1.0)]
    pub timepoint: f64,
    /// Computes the unavailability.
    /// See Table XI-2 of the Fault Tree Handbook. U.S. Nuclear Regulatory Commission (1981) for more information.
    #[arg(long, default_value_t = false)]
    pub unavailability: bool,
    /// Execution configuration parameters.
    #[command(flatten)]
    pub config: ExtraArgs,
}

#[derive(Parser, Debug, Clone)]
pub struct TranslateCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    pub input: String,
    /// Output file, writes a .wcnf file.
    #[arg(short, long)]
    pub output: String,
    /// If provided, the weights will be written to a separate file.
    #[arg(long)]
    pub w_file: Option<String>,
    /// Compute TEP of the FT a given timepoint.
    #[arg(short, long, default_value_t = 1.0)]
    pub timepoint: f64,
    /// Computes the unavailability.
    /// See Table XI-2 of the Fault Tree Handbook. U.S. Nuclear Regulatory Commission (1981) for more information.
    #[arg(long, default_value_t = false)]
    pub unavailability: bool,
    /// Execution configuration parameters.
    #[command(flatten)]
    pub config: ExtraArgs,
}

#[derive(Parser, Debug, Clone)]
pub struct ModCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    pub input: String,
    /// Solver path and arguments.
    #[arg(short, long)]
    pub solver_path: String,
    /// Compute TEP of the FT a given timepoint.
    #[arg(short, long, default_value_t = 1.0)]
    pub timepoint: f64,
    /// Execution configuration parameters.
    #[command(flatten)]
    pub config: ExtraArgs,
}

#[derive(Parser, Debug, Clone)]
pub struct ImportanceCommand {
    /// Input file containing the fault tree in GALILEO format.
    #[arg(short, long, required = true)]
    pub input: String,
    /// Solver path and arguments.
    #[arg(short, long)]
    pub solver_path: String,
    /// Timepoint to compute the true TEP and the measures for each basic event.
    #[arg(short, long, default_value_t = 1.0)]
    pub timepoint: f64,
    /// Execution configuration parameters.
    #[command(flatten)]
    pub config: ExtraArgs,
}

#[derive(Parser, Debug, Clone)]
pub struct ExtraArgs {
    /// Max cache size to distribute between the threads in KB. [default: 3500]
    #[arg(long, default_value_t = 3500)]
    pub max_cache_size: usize,
    /// Negate top gate if is an OR, to favor UnitPropagation. [default: false]
    #[arg(short, long, default_value_t = false)]
    pub negate_or: bool,
    /// Execution timeout for the WMC solver in seconds. [default: 300]
    #[arg(long, default_value_t = 300)]
    pub timeout_s: u64,
    /// Output format for the CNF formula. The format gives the extension to the file. Support values 'MC21' and 'MCC' [default: 'MC21']
    #[arg(long, default_value = "MC21")]
    pub format: String,
    /// Number of threads to use. [default: 1]
    #[arg(long, default_value_t = 1)]
    pub num_threads: usize,
    /// Verbosity in the output, if true prints more information. [default: false]
    #[arg(long, default_value_t = false)]
    pub verb: bool,
    /// Display progress messages and, if possible, progress bars. [default: false]
    #[arg(long, default_value_t = false)]
    pub display: bool,
    /// Simplify the FT by removing one children gates. [default: true]
    #[arg(long, default_value_t = true)]
    pub simplify: bool,
    /// If provided, postprocess the CNF formula by passing a CNF preprocessor. [default: None]
    #[arg(long, default_value = None)]
    pub preprocess: Option<String>,
}

/// CMD Arguments
#[derive(Parser, Debug, Clone)]
pub struct RandomGenerationCommand {
    /// Total number of nodes.
    #[arg(short, long)]
    pub n_nodes: usize,
    /// Output file, writes .dft files.
    #[arg(short, long)]
    pub output: String,
    /// Number in (0,1]. Rate of Basic Events from the total amount of Nodes.
    #[arg(long, default_value_t = 0.5)]
    pub rate_be: f64,
    /// Number in (0,1]. Rate of AND gates from the total amount of gates.
    #[arg(long, default_value_t = 0.5, requires = "rate_or")]
    pub rate_and: f64,
    /// Number in (0,1]. Rate of OR gates from the total amount of gates.
    #[arg(long, default_value_t = 0.5, requires = "rate_and")]
    pub rate_or: f64,
    /// Number in (0,1]. Rate of VOT gates from the total amount of gates.
    /// If the FT contains a VOT gate, cannot be processed and solved, Can be solved but after writing to dft.
    #[arg(long, default_value_t = 0.0, conflicts_with = "solver_path")]
    pub rate_vot: f64,
    ///Number in (0,1], multiplies the random generated float of the probability, so it can be a smaller probability.
    #[arg(long, default_value_t = 1e-4)]
    pub prob_multiplier: f64,
    /// In which percentage of the last gates start to put the Basic Events if they were not used before.
    #[arg(long, default_value_t = 0.6)]
    pub perc_last: f64,
    /// Specify the max number of children that a gate can have.
    #[arg(long, default_value_t = 5)]
    pub max_n_children: usize,
    /// Execution timeout for the WMC solver in seconds.
    #[arg(long, default_value_t = 300, requires = "solver_path")]
    pub timeout_s: u64,
    /// In which percentage of the last gates start to put the Basic Events if they were not used before.
    #[arg(long, value_parser = clap::value_parser!(u64))]
    pub seed: Option<u64>,
    /// Solver path.
    #[arg(short, long, conflicts_with = "rate_vot")]
    pub solver_path: Option<String>,
}
