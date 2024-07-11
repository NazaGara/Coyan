use crate::fault_trees::FaultTree;
use crate::formula::CNFFormat;
use itertools::Itertools;
use std::io::Write;
use std::process::{Command, Output, Stdio};

pub trait Solver {
    fn _name(&self) -> String;
    fn get_command(&self) -> String;
    fn run_model(&self, ft: &FaultTree<String>, format: CNFFormat, timebound: f64) -> Output;
    fn get_tep(&self, result: Output) -> f64;
    fn compute_probabilty(&self, ft: &FaultTree<String>, format: CNFFormat, timebound: f64) -> f64 {
        let out = self.run_model(ft, format, timebound);
        self.get_tep(out)
    }
}

pub fn get_solver_from_path(path: &str) -> Box<dyn Solver> {
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
            decow: 50,
            tmpdir: String::from(".tmp"),
            precision: 20,
            cs: 3500,
        }
    }

    fn _set_cache_size(&mut self, new_cs: usize) {
        self.cs = new_cs
    }
}

impl Solver for SharpsatTDSolver {
    fn _name(&self) -> String {
        String::from("SharpSAT-TD")
    }

    fn get_command(&self) -> String {
        format!(
            "{} -WE -decot {} -decow {} -tmpdir {} -prec {} -cs {}",
            self.path, self.decot, self.decow, self.tmpdir, self.precision, self.cs
        )
    }

    fn run_model(&self, ft: &FaultTree<String>, format: CNFFormat, timebound: f64) -> Output {
        let tmp_ft_file = format!("{}/tmp_ft.cnf", self.tmpdir);
        ft.dump_cnf_to_file(tmp_ft_file.clone(), format, timebound, None);

        let solver_cmd = format!("{} ./{}", self.get_command(), tmp_ft_file);
        let mut c: Command = Command::new("sh");
        c.arg("-c");
        c.arg(solver_cmd.clone());
        c.output().expect("failed to execute solver")
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
            cs: 4000,
            prec: 15,
        }
    }

    fn _set_cache_size(&mut self, new_cs: usize) {
        self.cs = new_cs
    }
}

impl Solver for GPMCSolver {
    fn _name(&self) -> String {
        String::from("GPMC")
    }

    fn get_command(&self) -> String {
        format!(
            "{} -mode={} -cs={} -prec={}",
            self.path, self.mode, self.cs, self.prec
        )
    }

    fn run_model(&self, ft: &FaultTree<String>, format: CNFFormat, timebound: f64) -> Output {
        let solver_cmd = self.get_command();
        let model_text = ft.dump_cnf(format, timebound);
        let mut child = Command::new("sh")
            .arg("-c")
            .arg(solver_cmd)
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("Failed to spawn child process");

        let mut stdin = child.stdin.take().expect("Failed to open stdin");
        std::thread::spawn(move || {
            stdin
                .write_all(model_text.as_bytes())
                .expect("Failed to write to stdin");
        });

        child.wait_with_output().expect("Failed to read stdout")
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

    fn get_command(&self) -> String {
        format!("{}", self.path)
    }

    fn run_model(&self, ft: &FaultTree<String>, format: CNFFormat, timebound: f64) -> Output {
        let solver_cmd = self.get_command();
        let model_text = ft.dump_cnf(format, timebound);
        let mut child = Command::new("sh")
            .arg("-c")
            .arg(solver_cmd)
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("Failed to spawn child process");

        let mut stdin = child.stdin.take().expect("Failed to open stdin");
        std::thread::spawn(move || {
            stdin
                .write_all(model_text.as_bytes())
                .expect("Failed to write to stdin");
        });

        child.wait_with_output().expect("Failed to read stdout")
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
