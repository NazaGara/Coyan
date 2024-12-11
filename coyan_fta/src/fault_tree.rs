use index_vec::IndexVec;
use indicatif::ParallelProgressIterator;
// use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use nodes::{Node, NodeId, NodeType};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::AtomicUsize;

use crate::fault_tree_normalizer::FaultTreeNormalizer;
use crate::formula::{CNFFormat, Formula};
use crate::modularizer::get_modules;
use crate::nodes;
use crate::preproc::*; //{BPlusE, Preprocessor, PMC};
use crate::solver::Solver;

impl<T> From<FaultTreeNormalizer<T>> for FaultTree<T> {
    fn from(ft_norm: FaultTreeNormalizer<T>) -> Self {
        FaultTree::<T> {
            nodes: ft_norm.nodes,
            root_id: ft_norm.root_id,
            node_counter: ft_norm.node_counter,
            negate_or: false,
        }
    }
}

impl<T> Clone for FaultTree<T>
where
    T: Clone,
{
    fn clone(&self) -> Self {
        FaultTree {
            nodes: self.nodes.clone(),
            root_id: self.root_id,
            node_counter: AtomicUsize::new(
                self.node_counter.load(std::sync::atomic::Ordering::Relaxed),
            ),
            negate_or: self.negate_or,
        }
    }
}

/// A Fault Tree representation in Coyan.
/// Can be created empty or read from a file using FT normalizer.
/// They have some extra details:
/// - Handle the logic of the Tseitin Encoding
/// - Do not have information about names of nodes
/// - Do not have VOT gates
/// - Do not have negations in arguments, but in separated gates.
pub struct FaultTree<T> {
    pub nodes: IndexVec<NodeId, Node<T>>,
    pub root_id: NodeId,
    node_counter: AtomicUsize,
    negate_or: bool,
}

impl FaultTree<String> {
    pub fn empty() -> Self {
        FaultTree {
            nodes: IndexVec::new(),
            root_id: NodeId::new(0),
            node_counter: AtomicUsize::new(0),
            negate_or: false,
        }
    }

    /// Generate a FT from a dft file.
    pub fn new_from_file(filename: &str, simplify: bool, negate_or: bool) -> Self {
        let mut ft_norm = FaultTreeNormalizer::new();
        ft_norm.read_from_file(filename, simplify);
        let mut ft = FaultTree::from(ft_norm);
        ft.negate_or = negate_or;
        ft
    }

    /// Internal method, changes the id of the root node.
    fn _set_root(&mut self, new_root_id: NodeId) {
        self.root_id = new_root_id;
    }

    /// Add the node to the fields of the struct.
    pub fn add_node(&mut self, node: Node<String>, nid: NodeId) {
        self.nodes.insert(nid, node)
    }

    fn new_id(&self) -> NodeId {
        NodeId::new(
            self.node_counter
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed),
        )
    }

    /// Update the roots of a Node in the nodes fields.
    pub fn update_roots(&mut self, new_node: Node<String>, nid: NodeId) {
        self.nodes.remove(nid);
        self.nodes.insert(nid, new_node);
    }

    /// Creates a new Fault Tree that is a submodule from the orignal one.
    /// It only contains the nodes that are used in the subtree
    pub fn subtree_with_root(&self, new_root_id: NodeId) -> FaultTree<String> {
        // Create New Fault Tree
        let mut sub_ft = FaultTree::empty();

        // Create a mapper to relate the old NodeIds to the new ones.
        let mut new_id_mapper = HashMap::new();

        // Create new ids for each children and add them to the subtree.
        let mut to_process = vec![new_root_id];
        while !to_process.is_empty() {
            let nid = to_process.pop().unwrap();
            let node = self.nodes[nid].clone();
            let children = node.get_children_ids();
            let new_nid = sub_ft.new_id();

            new_id_mapper.insert(nid, new_nid);
            sub_ft.add_node(node, new_nid);
            to_process.append(
                &mut children
                    .into_iter()
                    .filter(|c_id| !new_id_mapper.contains_key(c_id))
                    .collect(),
            );
        }

        // For each node, use the mapper to replace the node arguments.
        let _ = sub_ft
            .nodes
            .iter_mut()
            .map(|node| node.map_to_args(&new_id_mapper))
            .collect_vec();

        sub_ft
    }

    /// Return number of nodes in the tree.
    pub fn get_count(&self) -> usize {
        self.node_counter.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Get useful information of the Fault Tree.
    pub fn get_info(&self, _preprocess: Option<String>) -> (usize, usize, usize) {
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

    /// Apply the tseitin transformation to all the nodes in the tree.
    pub fn apply_tseitin(&self) -> Formula<NodeId> {
        let mut args = if self.nodes[self.root_id].kind.is_or() && self.negate_or {
            vec![Formula::Not(Box::new(Formula::Atom(self.root_id)))]
        } else {
            vec![Formula::Atom(self.root_id)]
        };

        for (nid, node) in self.nodes.iter_enumerated() {
            match node.tseitin_transformation(nid) {
                Formula::And(or_args) => args.extend(or_args),
                Formula::Or(literals) => args.push(Formula::Or(literals)),
                Formula::_True => {}
                _ => panic!("Something went wrong while translating the Tseitin transformation."),
            }
        }

        Formula::And(args)
    }

    /// Apply the tseitin transformation to all the nodes in the tree.
    pub fn apply_tseitin_used_only(&self) -> Formula<NodeId> {
        let mut args = if self.nodes[self.root_id].kind.is_or() && self.negate_or {
            vec![Formula::Not(Box::new(Formula::Atom(self.root_id)))]
        } else {
            vec![Formula::Atom(self.root_id)]
        };

        let mut to_process = vec![self.root_id];
        let mut seen = vec![self.root_id];
        while !to_process.is_empty() {
            let nid = to_process.pop().unwrap();
            let node = self.nodes[nid].clone();
            let mut children = node.get_children_ids();

            seen.append(&mut children);

            let mut unseen_children = children
                .into_iter()
                .filter(|cid| seen.contains(cid))
                .collect_vec();
            to_process.append(&mut unseen_children);

            match node.tseitin_transformation(nid) {
                Formula::And(or_args) => args.extend(or_args),
                Formula::Or(literals) => args.push(Formula::Or(literals)),
                Formula::_True => {}
                _ => panic!("Something went wrong while translating the Tseitin transformation."),
            }
        }

        Formula::And(args)
    }

    /// Save the fault tree CNF formula into a .wcnf o .cnf file depending on the format.
    pub fn dump_cnf_to_file(
        &self,
        filename: String,
        format: CNFFormat,
        timepoint: f64,
        w_file: Option<String>,
        preprocess: Option<String>,
    ) {
        let (formula_cnf, weights) = self.implicit_formula(format, timepoint, preprocess);

        let mut f = File::create(filename).expect("unable to create file");
        f.write_all(&formula_cnf.as_bytes())
            .expect("Error writing the formula to file");
        match w_file {
            None => {
                f.write_all(&weights.as_bytes())
                    .expect("Error writing weights to file");
            }
            Some(w_filename) => {
                let mut w_f =
                    File::create(format!("{}.w", w_filename)).expect("unable to create file");
                w_f.write_all(&weights.as_bytes())
                    .expect("Error writing the BE weights to file");
            }
        }
    }

    /// Internal method for code reading.
    /// Produces both the implicit boolean formula the weights in the specified format.
    fn implicit_formula(
        &self,
        format: CNFFormat,
        timepoint: f64,
        preprocess: Option<String>,
    ) -> (String, String) {
        // let cnf_formula = self.apply_tseitin();
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

        let (gate_weights, be_weights) = self.get_weights(weight_start, timepoint);
        let be_weights = be_weights.join("\n");
        let gate_weights = gate_weights.join("\n");

        let mut formula_str = text_formula
            .replace(" ∧ ", " 0 \n")
            .replace(" V ", " ")
            .replace("(", "")
            .replace(")", "");
        formula_str.push_str(" 0 \n");

        let formula_cnf = match preprocess {
            Some(preprocessor_path) => {
                let preprocessor: Box<dyn Preprocessor> =
                    get_preprocessor_from_path(&preprocessor_path);
                preprocessor.execute(&problem_line, &formula_str)
            }
            None => format!("{}\n{}\n", problem_line, formula_str),
        };

        let formula_cnf = if formula_cnf
            .split("\n")
            .any(|l| l.starts_with("c Solved by preprocessing"))
        {
            format!("{}\n{}\n", problem_line, formula_str)
        } else {
            formula_cnf
        };

        let weights = format!("{}\n{}", be_weights, gate_weights);

        (formula_cnf, weights)
    }

    /// Dump the implicit formula in CNF format to a String.
    pub fn dump_cnf(
        &self,
        format: CNFFormat,
        timepoint: f64,
        preprocess: Option<String>,
    ) -> String {
        let (mut formula_cnf, weights) = self.implicit_formula(format, timepoint, preprocess);
        formula_cnf.push_str(&weights);

        formula_cnf
    }

    /// Gives the weights in DIMACS format for the Gates and of the BE respectively.
    fn get_weights(&self, weight_start: String, timepoint: f64) -> (Vec<String>, Vec<String>) {
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
                        1.0 - (-value * timepoint).exp()
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

    /// Compute the Importance measures, the Birnbaum Measure, the Improvement Potential and the Criticality Measure.
    pub fn importance_measures(
        &self,
        solver: &Box<dyn Solver + Sync>,
        format: CNFFormat,
        timepoint: f64,
        negate_or: bool,
    ) -> HashMap<String, (String, String, String)> {
        let true_tep = solver.compute_probabilty(&self, format, timepoint, 300, None, negate_or);
        let be_lookup_table: HashMap<String, NodeId> = self
            .nodes
            .iter_enumerated()
            .filter_map(|(nid, n)| match &n.kind {
                NodeType::BasicEvent(name, _, _) => Some((name.to_owned(), nid.to_owned())),
                _ => None,
            })
            .collect::<HashMap<String, NodeId>>();

        let be_list = self
            .nodes
            .iter()
            .filter_map(|n| match &n.kind {
                NodeType::BasicEvent(name, _, _) => Some(name.to_owned()),
                _ => None,
            })
            .collect_vec();

        // This has to be a Vec<HashMap<String, (String, String)>>
        be_list
            .par_iter()
            .map(|be_name| {
                let mut ft = self.clone();
                (
                    be_name.to_owned(),
                    ft.measure_be(
                        be_name.clone(),
                        &solver,
                        &be_lookup_table,
                        format,
                        timepoint,
                        negate_or,
                    ),
                )
            })
            .collect::<HashMap<String, (f64, f64, f64)>>()
            .into_iter()
            .map(|(be_name, (ib, perf_tep, prob))| {
                (
                    be_name,
                    (
                        format!("{}", ib),
                        format!("{}", perf_tep - true_tep),
                        format!("{}", ib * (prob / true_tep)),
                    ),
                )
            })
            .collect::<HashMap<String, (String, String, String)>>()
    }

    /// Method called by [self] in the importance_measures method to compute each measure for a specific basic event.
    fn measure_be(
        &mut self,
        comp_name: String,
        solver: &Box<dyn Solver + Sync>,
        lookup_table: &HashMap<String, NodeId>,
        format: CNFFormat,
        timepoint: f64,
        negate_or: bool,
    ) -> (f64, f64, f64) {
        let nid = lookup_table
            .get(&comp_name)
            .expect("The name of the component is not a leaf in the Tree")
            .clone();

        let kind = &self.nodes.get(nid).unwrap().kind;
        let method = kind.get_method();
        let og_prob: f64 = if method.to_lowercase().eq("prob") {
            kind.get_prob()
        } else {
            1.0 - (-kind.get_prob() * timepoint).exp()
        };

        let pos_node = Node::new(NodeType::BasicEvent(
            comp_name.to_owned(),
            String::from("prob"),
            1.0,
        ));
        self.update_root(pos_node, nid);
        let pos_tep = solver.compute_probabilty(&self, format, timepoint, 300, None, negate_or);

        let neg_node = Node::new(NodeType::BasicEvent(
            comp_name.to_owned(),
            String::from("prob"),
            0.0,
        ));
        self.update_root(neg_node, nid);
        let neg_tep = solver.compute_probabilty(&self, format, timepoint, 300, None, negate_or);

        // Revert changes. There is no need to revert, because there are different FTs.
        // let og_node = Node::new(
        //     NodeType::BasicEvent(comp_name.to_owned(), method.to_owned(), og_prob),
        // );
        // self.update_roots(og_node, nid);
        (pos_tep - neg_tep, pos_tep, og_prob)
    }

    /// Update a Node by replacing it with another one.
    pub fn update_root(&mut self, new_node: Node<String>, nid: NodeId) {
        self.nodes.remove(nid);
        self.nodes.insert(nid, new_node);
    }

    /// Call to the Modularization algorithm.
    pub fn modularize_ft(&mut self) -> Vec<NodeId> {
        get_modules(self)
    }

    /// Method to replace the computed modules (in the module_ids parameter) with basic events with the same probability of failure at the given timepoint.
    /// Be careful with the provided number of threads, for large models (~2000 basic events) is easy to run out of memory.
    pub fn replace_modules(
        &mut self,
        solver: &Box<dyn Solver + Sync>,
        module_ids: Vec<NodeId>,
        format: CNFFormat,
        timepoint: f64,
        timeout_s: u64,
        num_threads: usize,
        negate_or: bool,
        display: bool,
    ) {
        // Chunk size should be related to the FT, not to the #threads.
        // But, to exploit parallelism, it should hold that chunk_size > #num_threads
        let chunk_size = std::cmp::max(module_ids.len().div_ceil(num_threads), num_threads);
        if display {
            // Compute the modules by chunks, could be more efficient if we take consideration of depth
            for chunk in module_ids.chunks(chunk_size).into_iter() {
                let to_replace: Vec<(NodeId, Node<String>)> = chunk
                    .par_iter()
                    .panic_fuse()
                    .progress()
                    .map(|&mod_id| {
                        let mod_ft = self.subtree_with_root(mod_id);
                        let tep = solver.compute_probabilty(
                            &mod_ft, format, timepoint, timeout_s, None, negate_or,
                        );
                        let repl_node = Node::new(NodeType::BasicEvent(
                            format!("repl_node_{}", mod_id),
                            String::from("prob"),
                            tep,
                        ));
                        (mod_id, repl_node)
                    })
                    .collect();
                for (mod_id, repl_node) in to_replace.into_iter() {
                    self.update_root(repl_node, mod_id);
                }
            }
        } else {
            for chunk in module_ids.chunks(chunk_size).into_iter() {
                let to_replace: Vec<(NodeId, Node<String>)> = chunk
                    .par_iter()
                    .map(|&mod_id| {
                        let mod_ft = self.subtree_with_root(mod_id);
                        let tep = solver.compute_probabilty(
                            &mod_ft, format, timepoint, timeout_s, None, negate_or,
                        );
                        let repl_node = Node::new(NodeType::BasicEvent(
                            format!("repl_node_{}", mod_id),
                            String::from("prob"),
                            tep,
                        ));
                        (mod_id, repl_node)
                    })
                    .collect();
                for (mod_id, repl_node) in to_replace.into_iter() {
                    self.update_root(repl_node, mod_id);
                }
            }
        }
    }
}
