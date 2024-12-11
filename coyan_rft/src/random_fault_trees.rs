use coyan_fta::{fault_tree::*, fault_tree_normalizer::*, nodes::*};
use itertools::Itertools;
use rand::{
    rngs::StdRng,
    seq::{IteratorRandom, SliceRandom},
    Rng, SeedableRng,
};
use std::{collections::HashMap, fs::File, io::Write, ops::Index};
const EPSILON: f64 = f64::EPSILON; // 2.2204460492503131E-16f64

/// Configuration for the random FT. Each value represent the proportion of each type of node.
/// 0st value for Basic Events
/// 1st-2nd value for AND gate, and OR gate.
/// 3rd value for Vot gates.
#[derive(Debug)]
pub struct RFTConfig(f64, f64, f64, f64);

impl RFTConfig {
    pub fn from_vec(args: Vec<f64>) -> Self {
        assert!(args.len() == 4);
        let r_be = args[0];
        let r_and = args[1];
        let r_or = args[2];
        let r_vot = args[3];
        assert!(
            r_and + r_or + r_vot >= 1.0 - EPSILON && r_and + r_or + r_vot <= 1.0 + EPSILON,
            "Check the gates rates, make sure that SUM(gate_rates) = 1"
        );
        RFTConfig(r_be, r_and, r_or, r_vot)
    }
}

pub struct RFaultTree<T> {
    // ft: FaultTree<T>,
    ft: FaultTreeNormalizer<T>,
    _n_nodes: usize,
    _config: RFTConfig,
}

impl RFaultTree<String> {
    /// Method to extract the fault tree. Notice that consumes the object.
    pub fn extract_ft(self) -> FaultTree<String> {
        FaultTree::from(self.ft)
    }
    pub fn _number_nodes(&self) -> usize {
        self._n_nodes
    }

    /// Creates a new Random Fault Tree. Uses a custom method
    pub fn new_random(
        n_nodes: usize,
        config: RFTConfig,
        p_multipler: f64,
        perc_last: f64,
        seed: u64,
        max_number_children: usize,
    ) -> Self {
        let mut rng = StdRng::seed_from_u64(seed);
        assert!(config.0 < 1.0, "The rate of basic events can't be 1.");
        let n_be = (config.0 * n_nodes as f64) as usize;
        assert!(n_be > 1, "We need at least more than 2 Basic Events.");
        let n_gates = (n_nodes - n_be) - 1;
        let (p_and, p_or, _p_vot) = (config.1, config.2, config.3);

        // Create FT, generate BE and Gates.
        // let mut ft = FaultTree::new();
        let mut ft_norm = FaultTreeNormalizer::new();
        let basic_events = (0..n_be)
            .into_iter()
            .map(|i| format!("x{}", i))
            .collect_vec();

        let gates = (0..n_gates)
            .into_iter()
            .map(|i| format!("g{}", i))
            .collect_vec();

        // Copy gates, take first gates so root can have it as children.
        let copy = gates.clone();
        let elems = rng.gen_range(2..=4);
        let first_gates = gates.clone()[0..6].to_vec();
        let roots = first_gates.into_iter().choose_multiple(&mut rng, elems);

        // Create root node.
        let root_name = "root".to_owned();
        let val: f64 = rng.gen();
        let mut root_node: Node<String> = if val >= p_and {
            Node::new(NodeType::PlaceHolder(
                root_name.clone(),
                "and".to_owned(),
                roots,
            ))
        } else {
            Node::new(NodeType::PlaceHolder(
                root_name.clone(),
                "or".to_owned(),
                roots,
            ))
        };

        root_node.set_formula(&ft_norm.nodes);

        let nid = ft_norm.new_id();
        ft_norm.root_id = nid;
        ft_norm.add_node(root_name.to_string(), root_node, nid);

        let mut used_be = vec![];
        // For each gate, we take as roots other gates with the bigger id.
        // If id index is too large, fills with basic events.
        for (i, g_name) in copy.into_iter().enumerate() {
            let nid = ft_norm.new_id();
            let k = rng.gen_range(3..=max_number_children);
            let val: f64 = rng.gen();
            let ahead: usize = max_number_children.max(8);

            let mut numbers: Vec<usize> = (1..=ahead).collect(); // Indicates how much 'ahead' I can take a gate.
            numbers.shuffle(&mut rng);
            let idxs = numbers[0..k].into_iter().map(|j| j + i).collect_vec(); //take K index from numbers

            let roots = idxs
                .iter()
                .map(|idx| {
                    if *idx >= n_gates {
                        let be = basic_events.choose(&mut rng).unwrap().to_owned();
                        used_be.push(be.clone());
                        be
                    } else {
                        gates.index(*idx).to_owned()
                    }
                })
                .collect_vec();

            let mut gate = if val <= p_and {
                Node::new(NodeType::PlaceHolder(
                    g_name.to_owned(),
                    "and".to_owned(),
                    roots,
                ))
            } else if val <= p_and + p_or {
                Node::new(NodeType::PlaceHolder(
                    g_name.to_owned(),
                    "or".to_owned(),
                    roots,
                ))
            } else {
                let choose_k = rng.gen_range(2..roots.len());
                Node::new(NodeType::PlaceHolder(
                    g_name.to_owned(),
                    format!("{}of{}", choose_k, roots.len()),
                    roots,
                ))
            };
            gate.set_formula(&ft_norm.nodes);
            ft_norm.add_node(g_name.to_string(), gate, nid);
        }

        // Set probability for the basic events. Currently using discrete probabilities.
        let _ = basic_events
            .iter()
            .map(|be| {
                let nid = ft_norm.new_id();
                let p: f64 = rng.gen();
                let mut node = Node::new(NodeType::BasicEvent(
                    be.to_string(),
                    "prob".to_owned(),
                    p * p_multipler,
                ));
                node.set_formula(&ft_norm.nodes);
                ft_norm.add_node(be.to_string(), node, nid);
            })
            .collect_vec();

        // Take all the unused basic events, and put them from the (1-PERC_LAST)%
        // of gates.
        let unused_be = basic_events
            .iter()
            .filter_map(|be| {
                if used_be.contains(be) {
                    None
                } else {
                    Some(be.to_owned())
                }
            })
            .collect_vec();

        // Take last (1-PERC_LAST)% of gates
        let lasts_gates = gates[(n_gates as f64 * perc_last) as usize..].to_vec();

        // Put the unused basic events as children of these lasts gates.
        let _ = unused_be
            .iter()
            .map(|be| {
                let mut new_roots = vec![be.to_string()];
                let g = lasts_gates.iter().choose(&mut rng).unwrap();
                let nid = ft_norm.lookup_table.get(g).unwrap().to_owned();
                let gate = ft_norm.nodes.get(nid).unwrap();

                let op = match &gate.kind {
                    NodeType::PlaceHolder(_, op, r) => {
                        new_roots.extend(r.to_vec());
                        op.to_owned()
                    }
                    _ => panic!("This should not happen"),
                };
                let op = if op.contains("of") {
                    let (choose_k, n) = op.split("of").collect_tuple().unwrap();
                    let n: usize = n
                        .parse()
                        .expect("Something went wrong when parsing VOT gate.");
                    format!("{}of{}", choose_k, n + 1)
                } else {
                    op
                };
                let mut new_node = Node::new(NodeType::PlaceHolder(g.to_owned(), op, new_roots));
                new_node.set_formula(&ft_norm.nodes);
                ft_norm.update_roots(new_node, nid);
            })
            .collect_vec();

        // Fill placeholders rearrenges the gates and set the correct types.
        ft_norm.fill_placeholders(true, true);

        RFaultTree {
            ft: ft_norm,
            _n_nodes: n_nodes,
            _config: config,
        }
    }

    /// Save the fault tree CNF formula into a .dft.
    pub fn save_to_dft(&self, filename: String) {
        let reverse_lookup_table: HashMap<NodeId, String> = self
            .ft
            .lookup_table
            .iter()
            .map(|(k, v)| (v.clone(), k.clone()))
            .collect();
        let top_line = format!(
            "toplevel {};",
            reverse_lookup_table.get(&self.ft.root_id).unwrap()
        );

        let gates = self
            .ft
            .nodes
            .iter_enumerated()
            .filter_map(|(i, n)| match &n.kind {
                NodeType::And(_) => Some(format!(
                    "{} {};",
                    reverse_lookup_table.get(&i).unwrap(),
                    n.get_formula()._reduce_formula()._formula_to_dft()
                )),
                NodeType::Or(_) => Some(format!(
                    "{} {};",
                    reverse_lookup_table.get(&i).unwrap(),
                    n.get_formula()._reduce_formula()._formula_to_dft()
                )),
                NodeType::Vot(_, _) => Some(format!(
                    "{} {};",
                    reverse_lookup_table.get(&i).unwrap(),
                    n.get_formula()._reduce_formula()._formula_to_dft()
                )),
                _ => None,
            })
            .join("\n");

        let be = self
            .ft
            .nodes
            .iter()
            .filter_map(|n| match &n.kind {
                NodeType::BasicEvent(name, method, prob) => {
                    Some(format!("{} {}={};", name, method, prob))
                }
                _ => None,
            })
            .join("\n");

        let mut f = File::create(filename).expect("unable to create file");
        f.write_all(&top_line.as_bytes())
            .expect("Error writing problem line to file");
        f.write_all("\n".as_bytes())
            .expect("Error writing the formula to file");
        f.write_all(&gates.as_bytes())
            .expect("Error writing the Gate weights to file");
        f.write_all("\n".as_bytes())
            .expect("Error writing . to file");
        f.write_all(&be.as_bytes())
            .expect("Error writing the BE weights to file");
    }
}
