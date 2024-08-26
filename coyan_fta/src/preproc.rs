use itertools::Itertools;
use std::fs::{self, File};
use std::io::Write;
use std::process::{Command, Stdio};
use std::time::Instant;

pub trait Preprocessor {
    fn execute(&self, problem_line: &String, formula_cnf: &String) -> String;
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
// MAX Memory usage in megabytes: 2147483647.
// MAX CPU time allowed in seconds: 2147483647.
#[derive(Debug, Clone)]
pub struct PMC {
    /// Use the Luby restart sequence
    luby_restart: bool,
    /// Randomize the initial activity
    rnd_init: bool,
    /// Number of time where the preprocessing technique is iterated. [default: 1]
    iterations: usize,
    options: PMCOptions,
}

impl Default for PMC {
    fn default() -> Self {
        PMC {
            luby_restart: true,
            rnd_init: false,
            iterations: 10,
            options: PMCOptions::_eq_configuration(),
        }
    }
}

impl Preprocessor for PMC {
    // If something fails, it returns the normal CNF without preprocessing.
    fn execute(&self, problem_line: &String, formula_cnf: &String) -> String {
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
                    "New: {:?}. Elapsed: {:?}",
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

pub struct BPlusE {
    /// Use the Luby restart sequence
    luby_restart: bool,
    /// Randomize the initial activity
    rnd_init: bool,
    /// Limit the solver for definability (0 means no limit). [default: 0]
    lim_solver: i32,
    /// Limit the maximal number of authorized resolution. [default: 500]
    max_num_res: i32,
}

impl Default for BPlusE {
    fn default() -> Self {
        BPlusE {
            luby_restart: true,
            rnd_init: false,
            lim_solver: 0,
            max_num_res: 500,
        }
    }
}

impl Preprocessor for BPlusE {
    fn execute(&self, problem_line: &String, formula_cnf: &String) -> String {
        let time_start = Instant::now();
        let tmp_file = "tmp.cnf";

        let command = format!(
            "./B+E_linux {} {} -limSolver={} -max#Res={} {}",
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
            self.lim_solver,
            self.max_num_res,
            tmp_file
        );

        let mut f = File::create(String::from(tmp_file)).expect("unable to create file");
        f.write_all(problem_line.as_bytes())
            .expect("Error writing problem line to file");
        f.write_all(formula_cnf.as_bytes())
            .expect("Error writing the formula to file");

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
                    "New: {:?}. Elapsed: {:?}",
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
