use crate::fault_trees::FaultTree;
use crate::formula::CNFFormat;
use itertools::Itertools;
use std::io::Write;
use std::process::{Command, Output, Stdio};

#[derive(Debug)]
pub enum Solver {
    GPMC(String),
    ADDMC(String),
}

impl Solver {
    pub fn _name(&self) -> String {
        match self {
            Solver::ADDMC(_) => "ADDMC".to_owned(),
            Solver::GPMC(_) => "GPMC".to_owned(),
        }
    }

    pub fn from_vec(cmd: Vec<String>) -> Solver {
        if cmd[0].to_lowercase().contains("gpmc") {
            let cmd = cmd.join(" -");
            Solver::GPMC(cmd)
        } else if cmd[0].to_lowercase().contains("addmc") {
            let cmd = cmd.join(" ");
            Solver::ADDMC(cmd)
        } else {
            panic!("Unsupported solver. Try GPMC or ADDMC.")
        }
    }

    pub fn get_command(&self) -> String {
        match self {
            Solver::ADDMC(cmd) => cmd.to_owned(),
            Solver::GPMC(cmd) => cmd.to_owned(),
        }
    }

    pub fn run_model(&self, ft: &FaultTree<String>, format: CNFFormat, timebound: f64) -> Output {
        let solver = self.get_command();
        let model_text = ft.dump_stdin(format, timebound);
        let mut child = Command::new("sh")
            .arg("-c")
            .arg(solver)
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

    pub fn get_prob_from_result(&self, result: Output) -> f64 {
        let stdout =
            String::from_utf8(result.stdout).expect("failed to produce the stdout of the solver");
        let result_line = match self {
            Solver::GPMC(_) => stdout
                .split("\n")
                .filter(|l| l.starts_with("c s exact double"))
                .join(""),
            Solver::ADDMC(_) => stdout
                .split("\n")
                .filter(|l| l.starts_with("s wmc"))
                .join(""),
        };
        result_line
            .split(" ")
            .last()
            .expect("Something went wrong while reading solver output")
            .to_owned()
            .parse()
            .expect("Error while parsing value to float")
    }

    pub fn compute_probabilty(
        &self,
        ft: &FaultTree<String>,
        format: CNFFormat,
        timebound: f64,
    ) -> f64 {
        let out = self.run_model(ft, format, timebound);
        self.get_prob_from_result(out)
    }
}
