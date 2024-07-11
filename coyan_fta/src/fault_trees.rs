use index_vec::IndexVec;
use itertools::Itertools;
use nodes::{Node, NodeId, NodeType};
use std::fs::File;
use std::io::Write;
use std::sync::atomic::AtomicUsize;
use std::{collections::HashMap, fs::read_to_string};

use crate::formula::{CNFFormat, Formula};
use crate::nodes;

/// Helper reader function.
fn _read_lines(filename: &str) -> Vec<String> {
    read_to_string(filename)
        .unwrap()
        .lines()
        .map(String::from)
        .collect()
}

/// A Fault Tree representation in Rust.
pub struct FaultTree<T> {
    pub lookup_table: HashMap<String, NodeId>,
    pub nodes: IndexVec<NodeId, Node<T>>,
    pub root_id: NodeId,
    node_counter: AtomicUsize,
}

impl Clone for FaultTree<String> {
    fn clone(&self) -> Self {
        let count = self.node_counter.load(std::sync::atomic::Ordering::Relaxed);
        let node_counter = AtomicUsize::new(count);

        FaultTree {
            lookup_table: self.lookup_table.clone(),
            nodes: self.nodes.clone(),
            root_id: self.root_id.clone(),
            node_counter,
        }
    }
}

/// Implementation of the Fault Tree.
impl FaultTree<String> {
    pub fn new() -> Self {
        FaultTree {
            lookup_table: HashMap::new(),
            nodes: IndexVec::new(),
            root_id: NodeId::new(0),
            node_counter: AtomicUsize::new(0),
        }
    }

    pub fn new_id(&self) -> NodeId {
        NodeId::new(
            self.node_counter
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed),
        )
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

    /// Method that reads the file, and create a node for each of the lines in the file.
    /// It generates Basic Events and Placeholders for the gates.
    /// Keeps track of the AND/OR/VOT gates with only one root, so it can then be replaced.
    fn read_file(&mut self, filename: &str, simplify: bool) -> String {
        let lines = _read_lines(filename);
        let mut root_name: String = "System".to_owned();
        let mut replace_mapper: HashMap<String, String> = HashMap::new();

        for l in lines.clone() {
            match &l.split_whitespace().map(str::to_string).collect_vec()[..] {
                [comment, ..] if comment.eq("//") || comment.starts_with("//") => {}
                [toplevel, name, ..] if toplevel.to_lowercase().as_str() == "toplevel" => {
                    root_name = name.replace("\"", "").replace(";", "").to_string();
                }
                [name, op, _args @ ..] if op.as_str().to_lowercase() == "not" => {
                    let name = name.replace("\"", "").replace(";", "");
                    if self.lookup_table.contains_key(&name) {
                        panic!("Name of Gate {} already in use.", name)
                    }
                    let args = _args
                        .iter()
                        .filter_map(|a| {
                            if a.eq(";") {
                                None
                            } else {
                                Some(a.replace("\"", "").replace(";", ""))
                            }
                        })
                        .collect_vec();
                    let nid = self.new_id();
                    let node = Node::new(
                        NodeType::PlaceHolder(name.to_owned(), op.to_lowercase().to_string(), args)
                    );
                    self.add_node(name.to_string(), node, nid);
                }
                [_name, op, _args @ ..]
                    if op.as_str().to_lowercase() == "csp"
                        || op.as_str().to_lowercase() == "wsp"
                        || op.as_str().to_lowercase() == "hsp"
                        || op.as_str().to_lowercase() == "pand"
                        || op.as_str().to_lowercase() == "seq"
                        || op.as_str().to_lowercase() == "fdep" =>
                {
                    panic!(
                        "Unsupported type of gate: {}. Is {} a Static FT?",
                        op.as_str(),
                        filename
                    )
                }
                [name, op, _args @ ..]
                    if op.as_str().to_lowercase() == "or"
                        || op.as_str().to_lowercase() == "and"
                        || op.as_str().to_lowercase() == "xor"
                        || op.as_str().contains("of") =>
                {
                    let name = name.replace("\"", "").replace(";", "");
                    if self.lookup_table.contains_key(&name) {
                        panic!("Name of Gate {} already in use.", name)
                    }
                    let args = _args
                        .iter()
                        .filter_map(|a| {
                            if a.eq(";") {
                                None
                            } else {
                                Some(a.replace("\"", "").replace(";", ""))
                            }
                        })
                        .collect_vec();

                    if simplify {
                        if args.len() == 1 {
                            replace_mapper.insert(name.clone(), args.first().unwrap().to_string());
                            if root_name == name {
                                root_name = args.first().unwrap().to_string();
                            }
                        } else {
                            let nid = self.new_id();
                            let node = Node::new(
                                NodeType::PlaceHolder(
                                    name.to_owned(),
                                    op.to_lowercase().to_string(),
                                    args,
                                )
                            );
                            self.add_node(name.to_string(), node, nid);
                        }
                    } else {
                        let nid = self.new_id();
                        let node = Node::new(
                            NodeType::PlaceHolder(
                                name.to_owned(),
                                op.to_lowercase().to_string(),
                                args,
                            )
                            
                        );
                        self.add_node(name.to_string(), node, nid);
                    }
                }
                [name, args, _others @ ..] => {
                    let name = name.replace("\"", "").replace(";", "");
                    if self.lookup_table.contains_key(&name) {
                        panic!("Name of Basic Event {} already in use.", name)
                    }
                    let (method, value) = args.split("=")
                    .collect_tuple()
                    .expect("Check if word 'lambda' or 'prob' is attached to the equal sign and its the first parameter of the basic event.");
                    let prob: f64 = value
                        .replace(";", "")
                        .parse()
                        .expect("Error while parsing value for basic event");
                    let nid = self.new_id();
                    let node = Node::new(
                        NodeType::BasicEvent(name.to_owned(), method.to_owned(), prob)
                    );
                    self.add_node(name.to_string(), node, nid);
                }
                _ => {}
            };
        }
        if simplify {
            self.preprocess_placeholders(replace_mapper)
        };
        root_name
    }

    /// Method to make a preprocess of the placeholders, updates the nodes that point to
    /// unnecesary gates.
    fn preprocess_placeholders(&mut self, mapper: HashMap<String, String>) {
        let mut corrected_mapper: HashMap<String, String> = mapper
            .clone()
            .into_iter()
            .filter(|(_k, v)| !mapper.contains_key(v))
            .collect();

        for (k, v) in mapper.iter() {
            if !corrected_mapper.contains_key(k) {
                let mut mid_value = v.clone();
                while mapper.contains_key(&mid_value) {
                    mid_value = mapper.get(&mid_value).unwrap().to_owned();
                }
                corrected_mapper.insert(k.to_owned(), mid_value);
            }
        }

        let placeholder_nids = self
            .nodes
            .indices()
            .filter(|i| {
                let n = self.nodes.get(i.raw()).unwrap();
                match n.kind {
                    NodeType::PlaceHolder(_, _, _) => true,
                    _ => false,
                }
            })
            .collect_vec();

        placeholder_nids
            .iter()
            .map(|nid| {
                let n = self.nodes.get(nid.to_owned()).unwrap();
                match &n.kind {
                    NodeType::PlaceHolder(name, op, args) => {
                        let args_replaced = args
                            .iter()
                            .map(|a| {
                                if corrected_mapper.contains_key(a) {
                                    corrected_mapper.get(a).unwrap().to_owned()
                                } else {
                                    a.to_owned()
                                }
                            })
                            .collect_vec();
                        let node = Node::new(
                            NodeType::PlaceHolder(
                                name.to_owned(),
                                op.to_lowercase().to_string(),
                                args_replaced,
                            )
                        );

                        self.update_roots(node, *nid)
                    }
                    _ => panic!(""),
                }
            })
            .collect_vec();
    }

    /// Checks all the placeholders on the nodes Vector, then replaces each one
    /// with the correct node.
    pub fn fill_placeholders(&mut self, set_formula: bool, keep_vot: bool) {
        let placeholder_nids = self
            .nodes
            .indices()
            .filter(|i| {
                let n = self.nodes.get(i.raw()).unwrap();
                match n.kind {
                    NodeType::PlaceHolder(_, _, _) => true,
                    _ => false,
                }
            })
            .collect_vec();

        for nid in placeholder_nids {
            let nn = self.nodes.get(nid).unwrap();
            match &nn.kind {
                NodeType::PlaceHolder(_name, op, args) => {
                    let args_ids = args
                        .iter()
                        .map(|a| {
                            // println!("{:?}", a);
                            *self.lookup_table.get(a).expect("Cant find argument")
                        })
                        .collect_vec();

                    let mut node: Node<String> = if op.eq("or") {
                        Node::new(NodeType::Or(args_ids))
                    } else if op.eq("xor") {
                        Node::new(NodeType::Xor(args_ids))
                    } else if op.eq("and") {
                        Node::new(NodeType::And(args_ids))
                    } else if op.eq("not") {
                        assert_eq!(
                            args_ids.len(),
                            1,
                            "Something is wrong with the negated gate"
                        );
                        Node::new(
                            NodeType::Not(args_ids.first().unwrap().to_owned())
                        )
                    } else if op.contains("of") {
                        let split = op.split("of").collect_vec();
                        let k = match split[0].parse::<usize>() {
                            Err(e) => panic!("Error {} on the VOT gate {}.", e, op),
                            Ok(v) => v,
                        };
                        let n = match split[1].parse::<usize>() {
                            Err(e) => panic!("Error {} on the VOT gate {}.", e, op),
                            Ok(v) => v,
                        };

                        if args.len() != n || k < 1 || k > n {
                            panic!("Error on the VOT gate {}. Check K and N in the definition of the Gate (<K>of<N>).", op)
                        }

                        if n == k {
                            Node::new(NodeType::And(args_ids))
                        } else if k == 1 {
                            Node::new(NodeType::Or(args_ids))
                        } else {
                            let clause_size = n - k;

                            if keep_vot {
                                Node::new(NodeType::Vot(k as i64, args_ids))
                            } else {
                                let mut roots = args.clone();
                                let mut aux_ids = vec![];
                                let mut new_args;
                                while roots.len() > clause_size {
                                    let elem = roots.pop().unwrap();
                                    for subset in roots.iter().combinations(clause_size) {
                                        new_args =
                                            vec![self.lookup_table.get(&elem).unwrap().clone()];
                                        let idk = subset
                                            .iter()
                                            .map(|s| {
                                                self.lookup_table.get(s.to_owned()).unwrap().clone()
                                            })
                                            .collect_vec();
                                        new_args.extend(idk);
                                        let aux_gid = self.new_id();
                                        aux_ids.push(aux_gid);
                                        let mut aux_gate =
                                            Node::new(NodeType::Or(new_args));
                                        if set_formula {
                                            aux_gate.set_formula(&self.nodes)
                                        };
                                        self.add_node(
                                            format!("aux_gate_{}", aux_gid),
                                            aux_gate,
                                            aux_gid,
                                        )
                                    }
                                }
                                Node::new(NodeType::And(aux_ids))
                            }
                        }
                    } else {
                        panic!("Something went wrong while processing the gates.")
                    };
                    if set_formula {
                        node.set_formula(&self.nodes)
                    };
                    self.update_roots(node, nid);
                }
                _ => panic!("This should not happen"),
            }
        }
    }

    /// Method that calls to the file reader method, then, depending ong the flag,
    /// will call to compute the tseitin. Finishes updating the struct.
    pub fn read_from_file(&mut self, filename: &str, simplify: bool) {
        let root_name: String = self.read_file(filename, simplify);
        self.fill_placeholders(false, false);
        let root_id = self.lookup_table.get(&root_name).unwrap();
        self.root_id = *root_id;
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

    /// Add the node to the fields of the struct.
    pub fn add_node(&mut self, name: String, node: Node<String>, nid: NodeId) {
        self.lookup_table.insert(name.clone(), nid);
        self.nodes.insert(nid, node);
    }

    /// Update the roots of a Node in the nodes fields.
    pub fn update_roots(&mut self, new_node: Node<String>, nid: NodeId) {
        self.nodes.remove(nid);
        self.nodes.insert(nid, new_node);
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
