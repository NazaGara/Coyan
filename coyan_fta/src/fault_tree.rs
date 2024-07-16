use index_vec::IndexVec;
use itertools::Itertools;
use nodes::{Node, NodeId, NodeType};
use std::fs::File;
use std::io::Write;
use std::sync::atomic::AtomicUsize;
use std::fs::read_to_string;
// use std::collections::HashMap;

use crate::formula::{CNFFormat, Formula};
use crate::nodes;
use crate::fault_tree_normalizer::FaultTreeNormalizer;

/// Helper reader function.
fn _read_lines(filename: &str) -> Vec<String> {
    read_to_string(filename)
        .unwrap()
        .lines()
        .map(String::from)
        .collect()
}

impl<T> From<FaultTreeNormalizer<T>> for FaultTree::<T>{
    fn from(ft_norm: FaultTreeNormalizer<T>) -> Self {
        FaultTree::<T> {
            // lookup_table: ft_norm.lookup_table,
            nodes: ft_norm.nodes,
            root_id: ft_norm.root_id,
            node_counter: ft_norm.node_counter,
        }
    }
}

/// A Fault Tree representation in Rust.
pub struct FaultTree<T> {
    // lookup_table: HashMap<String, NodeId>,
    nodes: IndexVec<NodeId, Node<T>>,
    pub root_id: NodeId,
    node_counter: AtomicUsize,
}


/// Internal representation of a Fault Tree.
/// Can be created empty or read from a file using the Normalizer.
/// They have some extra charactersistics:
/// - Handle the logic of the Tseitin Encoding
/// - Do not have information about names of nodes
/// - Do not have VOT gates
/// - DO not have negations in argujments, but in separeted gates.
impl FaultTree<String> {
    pub fn _empty() -> Self {
        FaultTree {
            nodes: IndexVec::new(),
            root_id: NodeId::new(0),
            node_counter: AtomicUsize::new(0),
        }
    }

    pub fn new_from_file(filename: &str, simplify: bool) -> Self {
        let mut ft_norm = FaultTreeNormalizer::new();
        ft_norm.read_from_file(filename, simplify);
        FaultTree::from(ft_norm)
    }

    pub fn get_count(&self) -> usize {
        self.node_counter.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Get useful information of the Fault Tree.
    pub fn get_info(&self) -> (usize, usize, usize) {
        let num_be = self
            .nodes
            .indices()
            .filter(|i| {
                let n = self.nodes.get(i.raw()).unwrap();
                match n.kind {
                    NodeType::BasicEvent(_, _, _) => true,
                    _ => false,
                }
            })
            .count();
        let num_gates = self
            .nodes
            .indices()
            .filter(|i| {
                let n = self.nodes.get(i.raw()).unwrap();
                match n.kind {
                    NodeType::And(_) => true,
                    NodeType::Or(_) => true,
                    NodeType::Vot(_, _) => true,
                    NodeType::Not(_) => true,
                    _ => false,
                }
            })
            .count();

        let f = self.apply_tseitin();

        (num_be, num_gates, f.num_clauses())
    }

    /// For each node on the tree, call to the tseitin transformation.
    pub fn apply_tseitin(&self) -> Formula<NodeId> {
        let mut args = vec![Formula::Atom(self.root_id)];
        for (nid, node) in self.nodes.iter_enumerated() {
            match node.tseitin_transformation(nid) {
                Formula::And(or_args) => args.extend(or_args),
                Formula::Or(literals) => args.push(Formula::Or(literals)),
                Formula::_True => {}
                _ => panic!("Something went wrong while translating the Tseitin transformation."),
            }
        }
        // args.sort_unstable();
        // args.dedup();
        Formula::And(args)
    }

    /// Save the fault tree CNF formula into a .wcnf o .cnf file depending on the format.
    pub fn dump_cnf_to_file(
        &self,
        filename: String,
        format: CNFFormat,
        timebound: f64,
        w_file: Option<String>,
    ) {
        let cnf_formula = self.apply_tseitin();
        let text_formula = cnf_formula.to_text();
        let n_vars = self.get_count();

        let n_clauses = cnf_formula.num_clauses();
        let (problem_line, weight_start) = match format {
            CNFFormat::MC21 => (
                format!("p cnf {} {}\n", n_vars, n_clauses),
                String::from("c p weight"),
            ),
            CNFFormat::MCC => (
                format!("p wcnf {} {}\n", n_vars, n_clauses),
                String::from("w"),
            ),
        };

        let (gate_weights, be_weights) = self.get_weights(weight_start, timebound);
        let be_str = be_weights.join("\n");
        let gate_str = gate_weights.join("\n");

        let mut formula_str = text_formula
            .replace(" ∧ ", " 0 \n")
            .replace(" V ", " ")
            .replace("(", "")
            .replace(")", "");
        formula_str.push_str(" 0 \n");

        let mut f = File::create(filename).expect("unable to create file");
        f.write_all(&problem_line.as_bytes())
            .expect("Error writing problem line to file");
        f.write_all(&formula_str.as_bytes())
            .expect("Error writing the formula to file");
        match w_file {
            None => {
                f.write_all(&be_str.as_bytes())
                    .expect("Error writing the BE weights to file");
                f.write_all(&"\n".as_bytes())
                    .expect("Error writing . to file");
                f.write_all(&gate_str.as_bytes())
                    .expect("Error writing the Gate weights to file");
            }
            Some(w_filename) => {
                let mut w_f =
                    File::create(format!("{}.w", w_filename)).expect("unable to create file");
                w_f.write_all(&be_str.as_bytes())
                    .expect("Error writing the BE weights to file");
                w_f.write_all(&"\n".as_bytes())
                    .expect("Error writing . to file");
                w_f.write_all(&gate_str.as_bytes())
                    .expect("Error writing the Gate weights to file");
            }
        }
    }

    pub fn dump_cnf(&self, format: CNFFormat, timebound: f64) -> String {
        let cnf_formula = self.apply_tseitin();
        let text_formula = cnf_formula.to_text();
        let n_vars = self.get_count();

        let n_clauses = cnf_formula.num_clauses();
        let (problem_line, weight_start) = match format {
            CNFFormat::MC21 => (
                format!("p cnf {} {}\n", n_vars, n_clauses),
                String::from("c p weight"),
            ),
            CNFFormat::MCC => (
                format!("p wcnf {} {}\n", n_vars, n_clauses),
                String::from("w"),
            ),
        };

        let (gate_weights, be_weights) = self.get_weights(weight_start, timebound);
        let be_str = be_weights.join("\n");
        let gate_str = gate_weights.join("\n");

        let mut formula_str = text_formula
            .replace(" ∧ ", " 0 \n")
            .replace(" V ", " ")
            .replace("(", "")
            .replace(")", "");
        formula_str.push_str(" 0 \n");

        let mut result = "".to_owned();
        result.push_str(&problem_line);
        result.push_str(&formula_str);
        result.push_str(&be_str);
        result.push_str(&"\n");
        result.push_str(&gate_str);

        result
    }

    /// Gives the weights in DIMACS format for the Gates and of the BE respectively.
    fn get_weights(&self, weight_start: String, timebound: f64) -> (Vec<String>, Vec<String>) {
        let n_vars = self.get_count();
        // We see which NodeId do not have the weights and set them to 1.
        // Only the BasicEvents Nodes have a fixed weight.
        let used_ids = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| match n.kind {
                NodeType::BasicEvent(_, _, _) => Some(i),
                _ => None,
            })
            .collect_vec();

        // Obtain weights of the gates
        let gate_weights: Vec<String> = (0..n_vars)
            .into_iter()
            .filter_map(|i| {
                if !used_ids.contains(&i) {
                    Some(format!(
                        "{} {} 1 0\n{} -{} 1 0",
                        weight_start,
                        i + 1,
                        weight_start,
                        i + 1
                    ))
                } else {
                    None
                }
            })
            .collect_vec();

        // Obtain weights of the basic Events.
        let be_weights = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| match &n.kind {
                NodeType::BasicEvent(_, method, value) => {
                    let prob: f64 = if method.to_lowercase().eq("prob") {
                        *value
                    } else if method.to_lowercase().eq("lambda") {
                        1.0 - (-value * timebound).exp()
                    } else {
                        panic!("Unsupported distribution of Basic Events. Try 'lambda' (Exponential) or 'prob' (Discrete).")
                    };
                    Some(format!(
                        "{} {} {} 0\n{} -{} {} 0",
                        weight_start,
                        i + 1,
                        prob,
                        weight_start,
                        i + 1,
                        1.0 - prob,
                    ))  
                }
                _ => None,
            })
            .collect_vec();
        (gate_weights, be_weights)
    }
}
