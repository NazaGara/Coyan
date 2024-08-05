use crate::fault_tree::FaultTree;
use crate::formula::CNFFormat;
use itertools::Itertools;
use rand::distributions::Alphanumeric;
use rand::Rng;
use std::fs::{self, File};
use std::io::Write;
use std::process::{Command, Output, Stdio};
use std::time::Instant;

pub trait Solver {
    fn _name(&self) -> String;
    fn get_command(&self, timeout_s: u64) -> String;
    fn run_model(
        &self,
        ft: &FaultTree<String>,
        format: CNFFormat,
        timebound: f64,
        timeout_s: u64,
        preprocess: bool,
    ) -> Result<Output, &'static str>;
    fn get_tep(&self, result: Output) -> f64;
    fn compute_probabilty(
        &self,
        ft: &FaultTree<String>,
        format: CNFFormat,
        timebound: f64,
        timeout_s: u64,
        preprocess: bool,
    ) -> f64 {
        let top_is_or = ft.nodes[ft.root_id].kind.is_or();
        match self.run_model(ft, format, timebound, timeout_s, preprocess) {
            Ok(value) => {
                let wmc_res = self.get_tep(value);
                if top_is_or {
                    1.0 - wmc_res
                } else {
                    wmc_res
                }
            }
            Err(msg) => panic!("{:?}", msg),
        }
    }
    fn _set_cache_size(&mut self, new_cs: usize);
}

pub fn get_solver_from_path(path: &str) -> Box<dyn Solver + Sync> {
    match path.to_ascii_lowercase() {
        x if x.contains("sharpsat") => Box::new(SharpsatTDSolver::new(path)),
        x if x.contains("addmc") => Box::new(ADDMCSolver::new(path)),
        x if x.contains("gpmc") => Box::new(GPMCSolver::new(path)),
        _ => panic!("Solver not supported."),
    }
}

/// Struct to support the solver [SharpSAT-TD](https://github.com/Laakeri/sharpsat-td)
/// The description of the flags is taken from the repository.
pub struct SharpsatTDSolver {
    /// Path to the solver
    path: String,
    ///decot -> the number of seconds to run flowcutter to find a tree decomposition. Required. Recommended value 60-600 if running with a total time budjet of 1800-3600 seconds.
    decot: usize,
    ///decow -> the weight of the tree decomposition in the decision heuristic. Recommended value >1 if the heuristic should care about the tree decomposition.
    decow: usize,
    ///tpmdir -> the directory to store temporary files for running flowcutter. Required.
    tmpdir: String,
    ///prec -> the number of digits in output of weighted model counting. Does not affect the internal precision.
    precision: usize,
    ///cs -> limit of the cache size. If the memory upper bound is X megabytes, then the value here should be around x/2-500.
    cs: usize,
}

impl SharpsatTDSolver {
    pub fn new(path: &str) -> Self {
        SharpsatTDSolver {
            path: String::from(path),
            // we: true,
            decot: 1,
            decow: 1,
            tmpdir: String::from(".tmp"),
            precision: 20,
            cs: 3500,
        }
    }
}

impl Solver for SharpsatTDSolver {
    fn _name(&self) -> String {
        String::from("SharpSAT-TD")
    }

    fn get_command(&self, timeout_s: u64) -> String {
        format!(
            "timeout -s KILL {}s {} -WE -decot {} -decow {} -tmpdir {} -prec {} -cs {}",
            timeout_s, self.path, self.decot, self.decow, self.tmpdir, self.precision, self.cs
        )
    }

    fn _set_cache_size(&mut self, new_cs: usize) {
        self.cs = new_cs
    }

    fn run_model(
        &self,
        ft: &FaultTree<String>,
        format: CNFFormat,
        timebound: f64,
        timeout_s: u64,
        preprocess: bool,
    ) -> Result<Output, &'static str> {
        // Set unique tmp name for each thread. With 5 char the chance of taking a name in use is 26⁵.
        let rnd_ft_file: String = rand::thread_rng()
            .sample_iter(&Alphanumeric)
            .take(5)
            .map(char::from)
            .collect();
        let tmp_ft_file = format!("{}/{}", self.tmpdir, rnd_ft_file);

        ft.dump_cnf_to_file(tmp_ft_file.clone(), format, timebound, None, preprocess);
        let solver_cmd = format!("{} ./{}", self.get_command(timeout_s), tmp_ft_file);

        let child = Command::new("sh")
            .arg("-c")
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .arg(solver_cmd)
            .spawn()
            .expect("Failed to spawn child process");

        match child.wait_with_output() {
            Ok(out) => {
                let stderr = String::from_utf8(out.stderr.clone())
                    .unwrap()
                    .to_lowercase();
                if stderr.eq("") {
                    // empty means nothing went wrong. Cleand and go next.
                    let _ = fs::remove_file(tmp_ft_file);
                    Ok(out)
                } else if stderr.eq("killed\n") {
                    // If it has something, check if is the killed signal
                    Err("Execution timeout.")
                } else {
                    ft.dump_cnf_to_file(String::from("failed.dft"), format, timebound, None, false);
                    println!("{:?}", stderr);
                    Err("Something went wrong...")
                }
            }
            Err(_err) => Err("Solver Process had an error."),
        }
    }

    fn get_tep(&self, result: Output) -> f64 {
        let stdout =
            String::from_utf8(result.stdout).expect("failed to produce the stdout of the solver");
        let result_line = stdout
            .split("\n")
            .filter(|l| l.starts_with("c s exact arb float"))
            .join("");

        result_line
            .split(" ")
            .last()
            .expect("Something went wrong while reading solver output")
            .to_owned()
            .parse()
            .expect("Error while parsing value to float")
    }
}

/// Struct to support the solver [GPMC](https://git.trs.css.i.nagoya-u.ac.jp/k-hasimt/GPMC)
/// The description of the flags is taken from the repository. There are more flags, ommited here.
pub struct GPMCSolver {
    /// Path to the solver
    path: String,
    /// -mode=<0..3> -> 0: MC; 1: WMC; 2: PMC; 3: WPMC
    mode: usize,
    /// -prec=<1..intmax> -> set the precision of floating-point numbers (default: 15).
    prec: usize,
    /// -cs=<1..intmax> -> set maximum component cache size (MB) (default: 4000)
    cs: usize,
}

impl GPMCSolver {
    pub fn new(path: &str) -> Self {
        GPMCSolver {
            path: String::from(path),
            mode: 1,
            cs: 3500,
            prec: 15,
        }
    }
}

impl Solver for GPMCSolver {
    fn _name(&self) -> String {
        String::from("GPMC")
    }

    fn _set_cache_size(&mut self, new_cs: usize) {
        self.cs = new_cs
    }

    fn get_command(&self, timeout_s: u64) -> String {
        format!(
            "timeout -s KILL {}s {} -mode={} -cs={} -prec={}",
            timeout_s, self.path, self.mode, self.cs, self.prec
        )
    }

    fn run_model(
        &self,
        ft: &FaultTree<String>,
        format: CNFFormat,
        timebound: f64,
        timeout_s: u64,
        preprocess: bool,
    ) -> Result<Output, &'static str> {
        let solver_cmd = self.get_command(timeout_s);
        let model_text = ft.dump_cnf(format, timebound, preprocess);
        let mut child = Command::new("sh")
            .arg("-c")
            .arg(solver_cmd)
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("Failed to spawn child process");

        let stdin = child.stdin.as_mut().expect("Failed to open stdin");
        stdin
            .write_all(model_text.as_bytes())
            .expect("Failed to write to stdin");

        match child.wait_with_output() {
            Ok(out) => {
                let stderr = String::from_utf8(out.stderr.clone())
                    .unwrap()
                    .to_lowercase();
                if stderr.eq("") {
                    Ok(out)
                } else if stderr.eq("killed\n") {
                    // If it has something, check if is the killed signal
                    Err("Execution timeout.")
                } else {
                    println!("{:?}", stderr);
                    Err("Something went wrong.")
                }
            }
            Err(_err) => Err("Solver Process had an error."),
        }
    }
    fn get_tep(&self, result: Output) -> f64 {
        let stdout =
            String::from_utf8(result.stdout).expect("failed to produce the stdout of the solver");

        let result_line = stdout
            .split("\n")
            .filter(|l| l.starts_with("c s exact double"))
            .join("");

        result_line
            .split(" ")
            .last()
            .expect("Something went wrong while reading solver output")
            .to_owned()
            .parse()
            .expect("Error while parsing value to float")
    }
}

pub struct ADDMCSolver {
    /// Path to the solver
    path: String,
}

impl ADDMCSolver {
    pub fn new(path: &str) -> Self {
        ADDMCSolver {
            path: String::from(path),
        }
    }
}

impl Solver for ADDMCSolver {
    fn _name(&self) -> String {
        String::from("ADDMC")
    }

    fn _set_cache_size(&mut self, _new_cs: usize) {
        println!("WARNING!: ADDMC solver does not have any parameter to regulate the cache size or any memory consumption.")
    }

    fn get_command(&self, timeout_s: u64) -> String {
        format!("timeout -s KILL {}s {}", timeout_s, self.path)
    }

    fn run_model(
        &self,
        ft: &FaultTree<String>,
        format: CNFFormat,
        timebound: f64,
        timeout_s: u64,
        preprocess: bool,
    ) -> Result<Output, &'static str> {
        let solver_cmd = self.get_command(timeout_s);
        let model_text = ft.dump_cnf(format, timebound, preprocess);
        let mut child = Command::new("sh")
            .arg("-c")
            .arg(solver_cmd)
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("Failed to spawn child process");

        let stdin = child.stdin.as_mut().expect("Failed to open stdin");
        stdin
            .write_all(model_text.as_bytes())
            .expect("Failed to write to stdin");

        match child.wait_with_output() {
            Ok(out) => {
                let stderr = String::from_utf8(out.stderr.clone())
                    .unwrap()
                    .to_lowercase();
                if stderr.eq("") {
                    // empty means nothing went wrong
                    Ok(out)
                } else if stderr.eq("killed\n") {
                    // If it has something, check if is the killed signal
                    Err("Execution timeout.")
                } else {
                    println!("{:?}", stderr);
                    Err("Something went wrong.")
                }
            }
            Err(_err) => Err("Solver Process had an error."),
        }
    }
    fn get_tep(&self, result: Output) -> f64 {
        let stdout =
            String::from_utf8(result.stdout).expect("failed to produce the stdout of the solver");

        let result_line = stdout
            .split("\n")
            .filter(|l| l.starts_with("s wmc"))
            .join("");

        result_line
            .split(" ")
            .last()
            .expect("Something went wrong while reading solver output")
            .to_owned()
            .parse()
            .expect("Error while parsing value to float")
    }
}

#[derive(Debug, Clone)]
#[warn(dead_code)]
pub struct PMCOptions {
    /// Remove some variable by affine function detection. [default: off]
    affine: bool,
    /// Remove some variable by or gates detection. [default: off]
    or_gate: bool,
    /// Remove some variable by equivalence detection. [default: off]
    equiv: bool,
    /// Perform the vivification technique. [default: on]
    vivification: bool,
    /// Symplify the proble w.r.t the backbone. [default: on]
    lit_implied: bool,
    /// Try to remove a maximum number of literal. [default: off]
    eliminate_lit: bool,
    /// Generate a set of clause which one will be added to the result. [default: off]
    add_clause: bool,
}

impl PMCOptions {
    pub fn _new(
        affine: bool,
        or_gate: bool,
        equiv: bool,
        vivification: bool,
        lit_implied: bool,
        eliminate_lit: bool,
        add_clause: bool,
    ) -> Self {
        PMCOptions {
            affine,
            or_gate,
            equiv,
            vivification,
            lit_implied,
            eliminate_lit,
            add_clause,
        }
    }
    pub fn _eq_configuration() -> Self {
        PMCOptions {
            affine: false,
            or_gate: false,
            equiv: false,
            vivification: true,
            lit_implied: true,
            eliminate_lit: true,
            add_clause: false,
        }
    }
    pub fn _numeq_configuration() -> Self {
        PMCOptions {
            affine: true,
            or_gate: true,
            equiv: true,
            vivification: true,
            lit_implied: true,
            eliminate_lit: true,
            add_clause: false,
        }
    }
    fn to_cmd(&self) -> String {
        format!(
            "-{} -{} -{} -{} -{} -{} -{}",
            if self.affine { "affine" } else { "no-affine" },
            if self.or_gate { "orGate" } else { "no-orGate" },
            if self.equiv { "equiv" } else { "no-equiv" },
            if self.vivification {
                "vivification"
            } else {
                "no-vivification"
            },
            if self.lit_implied {
                "litImplied"
            } else {
                "no-litImplied"
            },
            if self.eliminate_lit {
                "eliminateLit"
            } else {
                "no-eliminateLit"
            },
            if self.add_clause {
                "addClause"
            } else {
                "no-addClause"
            },
        )
    }
}

// ./preproc_linux [options] <input-file> <result-output-file>
// It can also read stdin :-).
#[derive(Debug, Clone)]
pub struct PreProccessor {
    luby_restart: bool,
    rnd_init: bool,
    iterations: usize,
    options: PMCOptions,
}

// MAX Memory usage in megabytes: 2147483647.
// MAX CPU time allowed in seconds: 2147483647.
impl PreProccessor {
    pub fn new(opt: PMCOptions) -> Self {
        PreProccessor {
            luby_restart: true,
            rnd_init: true,
            iterations: 10,
            options: opt,
        }
    }

    // If something fails, it returns the normal CNF without preprocessing.
    pub fn execute(&self, problem_line: &String, formula_cnf: &String) -> String {
        let time_start = Instant::now();
        let tmp_file = "tmp.cnf";

        let command = format!(
            "./preproc_linux -iterate={} {} {} {} {}",
            self.iterations,
            if self.luby_restart {
                "-luby"
            } else {
                "-no-luby"
            },
            if self.rnd_init {
                "-rnd-init"
            } else {
                "-no-rnd-init"
            },
            self.options.to_cmd(),
            tmp_file
        );

        let mut f = File::create(String::from(tmp_file)).expect("unable to create file");
        f.write_all(problem_line.as_bytes())
            .expect("Error writing problem line to file");
        f.write_all(formula_cnf.as_bytes())
            .expect("Error writing the formula to file");

        println!("Previous: {:?}", problem_line);

        let child = Command::new("sh")
            .arg("-c")
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .arg(command)
            .spawn()
            .expect("Failed to spawn child process");

        match child.wait_with_output() {
            Ok(out) => {
                let _ = fs::remove_file(String::from(tmp_file));
                let processed = String::from_utf8(out.stdout)
                    .expect("failed to produce the stdout of the solver");
                let elapsed = time_start.elapsed();
                println!(
                    "Previous: {:?}. Elapsed: {:?}",
                    processed
                        .clone()
                        .split("\n")
                        .into_iter()
                        .filter(|l| l.starts_with("p cnf"))
                        .collect_vec(),
                    elapsed
                );
                processed
            }
            Err(_err) => format!("{}\n{}\n", problem_line, formula_cnf),
        }
    }
}
