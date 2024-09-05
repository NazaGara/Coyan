use crate::fault_tree::FaultTree;
use crate::formula::CNFFormat;
use itertools::Itertools;
use rand::distributions::Alphanumeric;
use rand::Rng;
use std::fs;
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
        preprocess: Option<String>,
    ) -> Result<Output, &'static str>;
    fn get_tep(&self, result: Output) -> f64;
    fn compute_probabilty(
        &self,
        ft: &FaultTree<String>,
        format: CNFFormat,
        timepoint: f64,
        timeout_s: u64,
        preprocess: Option<String>,
        negate_top_or: bool,
    ) -> f64 {
        let top_is_or = ft.nodes[ft.root_id].kind.is_or();
        match self.run_model(ft, format, timepoint, timeout_s, preprocess) {
            Ok(value) => {
                let wmc_res = self.get_tep(value);
                if top_is_or && negate_top_or {
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
        x if x.contains("dmc") => Box::new(DMCSolver::new(path)),
        _ => panic!("Solver not supported. Supported solves: ADDMC - GPMC - SharpSAT-TD"),
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
        preprocess: Option<String>,
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
                // empty means nothing went wrong. Clean and go.
                if stderr.eq("") {
                    let _ = fs::remove_file(tmp_ft_file);
                    Ok(out)
                // If it has something, check if is the killed signal
                } else if stderr.eq("killed\n") {
                    Err("Execution timeout.")
                // Something else failed, print error.
                } else {
                    // ft.dump_cnf_to_file(String::from("failed.dft"), format, timebound, None, None);
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
        preprocess: Option<String>,
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

pub struct DMCSolver {
    /// Path to the solver. Required.
    dmc_path: String,
    /// Path to the htb tool. Required. It is assumed that htb and dmc are in the same directory.
    htb_path: String,
    ///tpmdir -> the directory to store temporary files for running the tree decomposition. Required.
    tmpdir: String,
}

impl DMCSolver {
    pub fn new(path: &str) -> Self {
        let dmc_path = String::from(path);
        let htb_path = format!(
            "{}{}",
            path.strip_suffix("dmc")
                .expect("DMC file does not end in dmc"),
            "htb"
        );
        DMCSolver {
            dmc_path,
            htb_path,
            tmpdir: String::from(".tmp"),
        }
    }

    fn compute_joint_tree(&self, timeout_s: u64, filepath: &str) -> (String, u64) {
        let mut c = Command::new("sh");
        let time_start = Instant::now();
        let cmd = format!(
            "timeout -s kill {} {} --cf ./{}",
            timeout_s,
            self.htb_path.to_owned(),
            filepath
        );
        c.arg("-c");
        c.arg(cmd);
        let out = c.output().expect("failed to execute tree decomposition");
        let duration = time_start.elapsed();
        (
            String::from_utf8(out.stdout).expect("failed to produce the stdout of the solver"),
            timeout_s - duration.as_secs(),
        )
    }
}

impl Solver for DMCSolver {
    fn _name(&self) -> String {
        String::from("DMC")
    }
    fn _set_cache_size(&mut self, _new_cs: usize) {
        println!("WARNING!: Memory consumption of the DMC solver is not implemented.")
    }

    fn get_command(&self, timeout_s: u64) -> String {
        format!("timeout -s KILL {}s {}", timeout_s, self.dmc_path)
    }

    fn run_model(
        &self,
        ft: &FaultTree<String>,
        format: CNFFormat,
        timebound: f64,
        timeout_s: u64,
        preprocess: Option<String>,
    ) -> Result<Output, &'static str> {
        let rnd_ft_file: String = rand::thread_rng()
            .sample_iter(&Alphanumeric)
            .take(5)
            .map(char::from)
            .collect();
        let tmp_ft_file = format!("{}/{}", self.tmpdir, rnd_ft_file);
        ft.dump_cnf_to_file(tmp_ft_file.clone(), format, timebound, None, preprocess);
        let (heuristic_tree, remaining_s) = self.compute_joint_tree(timeout_s, &tmp_ft_file);

        let solver_cmd: String = format!("{} --cf {}", self.get_command(remaining_s), tmp_ft_file);
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
                .write_all(heuristic_tree.as_bytes())
                .expect("Failed to write to stdin");
        });

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
        preprocess: Option<String>,
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
