use index_vec::IndexVec;
use itertools::Itertools;
use nodes::{Node, NodeId};
use std::sync::atomic::AtomicUsize;
use std::{collections::HashMap, fs::read_to_string};

use crate::nodes::{self, BasicEvent, RepairMode};

/// Helper reader function.
fn _read_lines(filename: &str) -> Vec<String> {
    read_to_string(filename)
        .unwrap()
        .lines()
        .map(String::from)
        .collect()
}

/// Normalizer struct for FTs.
/// It handles all the not so nice parsing and reading of the FT.
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
            root_id: self.root_id,
            node_counter,
        }
    }
}

impl Default for FaultTreeNormalizer<String> {
    fn default() -> Self {
        FaultTreeNormalizer {
            lookup_table: HashMap::new(),
            nodes: IndexVec::new(),
            root_id: NodeId::new(0),
            node_counter: AtomicUsize::new(0),
        }
    }
}

fn parse_basic_event(name: &str, args: &[String]) -> (String, BasicEvent) {
    let name = name.replace("\"", "").replace(";", "");
    let mut params = HashMap::new();

    for item in args {
        let (key, value) = item.split("=").collect_tuple().unwrap();
        let value = value
            .replace(";", "")
            .parse::<f64>()
            .unwrap_or_else(|_| panic!("Could not parse number {value}."));
        params.insert(key, value);
    }

    let be = if params.contains_key("prob") {
        BasicEvent::new_with_prob(*params.get("prob").unwrap())
    } else {
        let mut be = BasicEvent::new_with_rate(*params.get("lambda").expect(
            "Basic Event must have either a discrete or continuous distribution function.",
        ));
        if params.contains_key("repair") {
            be.with_repair_mode(RepairMode::Monitored(*params.get("repair").unwrap()));
        }
        be
    };

    (name, be)
}

impl FaultTreeNormalizer<String> {
    pub fn new_id(&self) -> NodeId {
        NodeId::new(
            self.node_counter
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed),
        )
    }

    /// Method that reads the file, and create a node for each of the lines in the file.
    /// Only create Basic Events and Placeholders.
    /// Keeps track of the gates with only one root (expect NOT), so later it can then be simplified.
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
                [name, op, args @ ..] if op.as_str().to_lowercase() == "not" => {
                    let name = name.replace("\"", "").replace(";", "");
                    if self.lookup_table.contains_key(&name) {
                        panic!("Name of Gate {} already in use.", name)
                    }
                    let args = args
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
                    let node =
                        Node::PlaceHolder(name.to_owned(), op.to_lowercase().to_string(), args);
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
                [name, op, args @ ..]
                    if op.as_str().to_lowercase() == "or"
                        || op.as_str().to_lowercase() == "and"
                        || op.as_str().to_lowercase() == "xor"
                        || op.as_str().contains("of") =>
                {
                    let name = name.replace("\"", "").replace(";", "");
                    if self.lookup_table.contains_key(&name) {
                        panic!("Name of Gate '{}' already in use.", name)
                    }
                    let args = args
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
                            let node = Node::PlaceHolder(
                                name.to_owned(),
                                op.to_lowercase().to_string(),
                                args,
                            );
                            self.add_node(name.to_string(), node, nid);
                        }
                    } else {
                        let nid = self.new_id();
                        let node =
                            Node::PlaceHolder(name.to_owned(), op.to_lowercase().to_string(), args);
                        self.add_node(name.to_string(), node, nid);
                    }
                }
                [name, args @ ..] => {
                    let (name, be) = parse_basic_event(name, args);
                    let nid = self.new_id();
                    let node = Node::BasicEvent(name.to_owned(), be);
                    self.add_node(name.to_string(), node, nid);
                }
                _ => {}
            };
        }
        if simplify {
            self.preprocess_placeholders(replace_mapper);
        };
        root_name
    }

    /// Method to make a preprocess of the placeholders, updating the nodes that point to
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
            .iter_enumerated()
            .filter_map(|(nid, node)| {
                if matches!(node, Node::PlaceHolder(_, _, _)) {
                    Some(nid)
                } else {
                    None
                }
            })
            .collect_vec();

        placeholder_nids
            .iter()
            .map(|nid| {
                let n = self.nodes.get(nid.to_owned()).unwrap();
                match &n {
                    Node::PlaceHolder(name, op, args) => {
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
                        let node = Node::PlaceHolder(
                            name.to_owned(),
                            op.to_lowercase().to_string(),
                            args_replaced,
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
    pub fn fill_placeholders(&mut self, keep_vot: bool) {
        let placeholder_nids = self
            .nodes
            .iter_enumerated()
            .filter_map(|(nid, node)| {
                if matches!(node, Node::PlaceHolder(_, _, _)) {
                    Some(nid)
                } else {
                    None
                }
            })
            .collect_vec();

        for nid in placeholder_nids {
            let nn = self.nodes.get(nid).unwrap();
            match &nn {
                Node::PlaceHolder(_name, op, args) => {
                    let args_ids = args
                        .iter()
                        .map(|a| {
                            *self
                                .lookup_table
                                .get(a)
                                .unwrap_or_else(|| panic!("Cant find argument {}", a))
                        })
                        .collect_vec();

                    let node: Node<String> = if op.eq("or") {
                        Node::Or(args_ids)
                    } else if op.eq("xor") {
                        Node::Xor(args_ids)
                    } else if op.eq("and") {
                        Node::And(args_ids)
                    } else if op.eq("not") {
                        assert_eq!(
                            args_ids.len(),
                            1,
                            "Something is wrong with the negated gate"
                        );
                        Node::Not(args_ids.first().unwrap().to_owned())
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
                            Node::And(args_ids)
                        } else if k == 1 {
                            Node::Or(args_ids)
                        } else {
                            let clause_size = n - k;

                            if keep_vot {
                                Node::Vot(k as i64, args_ids)
                            } else {
                                let mut roots = args.clone();
                                let mut aux_ids = vec![];
                                let mut new_args;
                                while roots.len() > clause_size {
                                    let elem = roots.pop().unwrap();
                                    for subset in roots.iter().combinations(clause_size) {
                                        new_args = vec![*self.lookup_table.get(&elem).unwrap()];
                                        new_args.extend(subset.iter().map(|s| {
                                            *self.lookup_table.get(s.to_owned()).unwrap()
                                        }));
                                        let aux_gid = self.new_id();
                                        aux_ids.push(aux_gid);
                                        let aux_gate = Node::Or(new_args);
                                        self.add_node(
                                            format!("aux_gate_{}", aux_gid),
                                            aux_gate,
                                            aux_gid,
                                        )
                                    }
                                }
                                Node::And(aux_ids)
                            }
                        }
                    } else {
                        panic!("Something went wrong while processing the gates.")
                    };
                    self.update_roots(node, nid);
                }
                _ => panic!("This should not happen"),
            }
        }
    }

    /// Public method to read the FT from the File, apply simplifications and replace the placeholders.
    pub fn read_from_file(&mut self, filename: &str, simplify: bool) {
        let root_name: String = self.read_file(filename, simplify);
        self.fill_placeholders(false);
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
