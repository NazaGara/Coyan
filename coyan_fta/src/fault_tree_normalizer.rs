use index_vec::IndexVec;
use itertools::Itertools;
use nodes::{Node, NodeId, NodeType};
use std::sync::atomic::AtomicUsize;
use std::{collections::HashMap, fs::read_to_string};

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
pub struct FaultTreeNormalizer<T> {
    pub lookup_table: HashMap<String, NodeId>,
    pub nodes: IndexVec<NodeId, Node<T>>,
    pub root_id: NodeId,
    pub node_counter: AtomicUsize,
}

impl Clone for FaultTreeNormalizer<String> {
    fn clone(&self) -> Self {
        let count = self.node_counter.load(std::sync::atomic::Ordering::Relaxed);
        let node_counter = AtomicUsize::new(count);

        FaultTreeNormalizer {
            lookup_table: self.lookup_table.clone(),
            nodes: self.nodes.clone(),
            root_id: self.root_id.clone(),
            node_counter,
        }
    }
}

/// Fault Tree implementation to normalize the input.
impl FaultTreeNormalizer<String> {
    pub fn new() -> Self {
        FaultTreeNormalizer {
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
                    let node = Node::new(NodeType::PlaceHolder(
                        name.to_owned(),
                        op.to_lowercase().to_string(),
                        args,
                    ));
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
                            let node = Node::new(NodeType::PlaceHolder(
                                name.to_owned(),
                                op.to_lowercase().to_string(),
                                args,
                            ));
                            self.add_node(name.to_string(), node, nid);
                        }
                    } else {
                        let nid = self.new_id();
                        let node = Node::new(NodeType::PlaceHolder(
                            name.to_owned(),
                            op.to_lowercase().to_string(),
                            args,
                        ));
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
                    let node = Node::new(NodeType::BasicEvent(
                        name.to_owned(),
                        method.to_owned(),
                        prob,
                    ));
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
                        let node = Node::new(NodeType::PlaceHolder(
                            name.to_owned(),
                            op.to_lowercase().to_string(),
                            args_replaced,
                        ));

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
                            *self
                                .lookup_table
                                .get(a)
                                .expect(&format!("Cant find argument {}", a))
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
                        Node::new(NodeType::Not(args_ids.first().unwrap().to_owned()))
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
                                        let mut aux_gate = Node::new(NodeType::Or(new_args));
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
}
