use crate::formula::Formula;
use index_vec::IndexVec;
use itertools::Itertools;
use std::fmt::{self, Debug, Display};

index_vec::define_index_type! {
    pub struct NodeId = usize;
    MAX_INDEX = usize::MAX;
    DISABLE_MAX_INDEX_CHECK = cfg!(not(debug_assertions));

}

impl Display for NodeId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self._raw + 1)
    }
}

/// Enum representing the Types of nodes we can have.
#[derive(Debug, Clone)]
pub enum NodeType<T> {
    /// Arguments: Name, method, prob
    BasicEvent(T, String, f64),
    /// Arguments: args of gate
    Not(NodeId), // Not(Box<Self>),
    /// Arguments: args of gate
    And(Vec<NodeId>),
    /// Arguments: args of gate
    Or(Vec<NodeId>),
    /// Arguments: args of gate
    Xor(Vec<NodeId>),
    /// Arguments: k, args of gate
    Vot(i64, Vec<NodeId>),
    /// Arguments: name, op, args of gate
    PlaceHolder(T, String, Vec<T>),
}

/// Implementation of the nodes.
impl<T> NodeType<T>
where
    T: Clone + Debug,
{
    pub fn get_method(&self) -> String {
        match self {
            NodeType::BasicEvent(_, m, _) => m.to_owned(),
            _ => panic!("Get Method only available for Basic Events"),
        }
    }
    pub fn get_prob(&self) -> f64 {
        match self {
            NodeType::BasicEvent(_, _, p) => p.clone(),
            _ => panic!("Get Method only available for Basic Events"),
        }
    }

    pub fn is_or(&self) -> bool {
        match self {
            NodeType::Or(_) => true,
            _ => false,
        }
    }

    /// Allows to obtain a formula from the nodetype by checking the vec of node ids.
    pub fn _get_formula(&self, nodes: &IndexVec<NodeId, Node<T>>) -> Formula<T> {
        if nodes.is_empty() {
            match self {
                NodeType::BasicEvent(name, _, _) => return Formula::Atom(name.to_owned()),
                NodeType::PlaceHolder(name, _, _) => return Formula::Atom(name.to_owned()),
                _ => panic!("Start creating filling the nodes with the Basic Events first."),
            }
        }
        match self {
            NodeType::PlaceHolder(name, _, _) => Formula::Atom(name.to_owned()),
            NodeType::BasicEvent(name, _, _) => Formula::Atom(name.to_owned()),
            NodeType::Not(arg) => Formula::Not(Box::new(nodes.get(*arg).unwrap().get_formula())),
            NodeType::And(args) => {
                if args.is_empty() {
                    panic!("Empty args")
                };
                Formula::And(
                    args.iter()
                        .map(|id| nodes.get(*id).unwrap().get_formula())
                        .collect(),
                )
            }
            NodeType::Or(args) => {
                if args.is_empty() {
                    panic!("Empty args")
                };
                Formula::Or(
                    args.iter()
                        .map(|id| nodes.get(*id).unwrap().get_formula())
                        .collect(),
                )
            }
            NodeType::Vot(k, args) => {
                if args.is_empty() {
                    panic!("Empty args")
                };
                Formula::Vot(
                    *k,
                    args.iter()
                        .map(|id| nodes.get(*id).unwrap().get_formula())
                        .collect(),
                )
            }
            NodeType::Xor(args) => {
                if args.is_empty() {
                    panic!("Empty args")
                };
                Formula::Or(
                    args.iter()
                        .map(|id| nodes.get(*id).unwrap().get_formula())
                        .collect(),
                )
            }
        }
    }
}

/// Struct for Nodes of the tree, each Node contains the type of Node and the formula.
/// The formula is used when creating random fault trees.
#[derive(Debug, Clone)]
pub struct Node<T> {
    pub kind: NodeType<T>,
    formula: Formula<T>,
}

/// Implementation of the Node Struct.
impl<T> Node<T>
where
    T: Debug + Clone,
{
    pub fn new(kind: NodeType<T>) -> Self {
        Node {
            kind,
            formula: Formula::_True,
        }
    }

    pub fn set_formula(&mut self, nodes: &IndexVec<NodeId, Node<T>>) {
        self.formula = self.kind._get_formula(nodes);
    }

    pub fn get_formula(&self) -> Formula<T> {
        self.formula.clone()
    }

    /// Apply the Tseitin rule for the AND NodeType.
    fn tseitin_and(&self, self_id: NodeId, args: &Vec<NodeId>) -> Formula<NodeId> {
        let mut clauses: Vec<Formula<_>> = args
            .iter()
            .map(|nid| {
                Formula::Or(vec![
                    Formula::Not(Box::new(Formula::Atom(self_id.clone()))),
                    Formula::Atom(nid.clone()),
                ])
            })
            .collect();
        let mut other: Vec<Formula<NodeId>> = args
            .iter()
            .map(|nid| Formula::Not(Box::new(Formula::Atom(nid.clone()))))
            .collect();
        other.push(Formula::Atom(self_id.clone()));
        clauses.push(Formula::Or(other));
        Formula::And(clauses)
    }

    /// Apply the Tseitin rule for the OR NodeType.
    fn tseitin_or(&self, self_id: NodeId, args: &Vec<NodeId>) -> Formula<NodeId> {
        let mut clauses: Vec<Formula<_>> = args
            .iter()
            .map(|nid| {
                Formula::Or(vec![
                    Formula::Atom(self_id.clone()),
                    Formula::Not(Box::new(Formula::Atom(nid.clone()))),
                ])
            })
            .collect();
        let mut other: Vec<Formula<NodeId>> =
            args.iter().map(|nid| Formula::Atom(nid.clone())).collect();
        other.push(Formula::Not(Box::new(Formula::Atom(self_id.clone()))));
        clauses.push(Formula::Or(other));
        Formula::And(clauses)
    }

    /// Apply the Tseitin rule for the XOR NodeType.
    fn tseitin_xor(&self, self_id: NodeId, args: &Vec<NodeId>) -> Formula<NodeId> {
        let mut clause_g_neg = args
            .iter()
            .map(|nid| Formula::Atom(nid.clone()))
            .collect_vec();
        clause_g_neg.push(Formula::Not(Box::new(Formula::Atom(self_id.clone()))));

        let mut clause_all_neg = args
            .iter()
            .map(|nid| Formula::Not(Box::new(Formula::Atom(nid.clone()))))
            .collect_vec();
        clause_all_neg.push(Formula::Not(Box::new(Formula::Atom(self_id.clone()))));

        let mut other_clauses = (0..args.len())
            .into_iter()
            .map(|_| {
                let mut c = args.iter().map(|nid| Formula::Atom(*nid)).collect_vec();
                c.push(Formula::Atom(self_id.clone()));
                c
            })
            .collect_vec();

        for (i, c) in other_clauses.iter_mut().enumerate() {
            let f = c.remove(i).negate();
            c.insert(i, f);
        }

        let mut formula = other_clauses
            .into_iter()
            .map(|c| Formula::Or(c))
            .collect_vec();
        formula.push(Formula::Or(clause_g_neg));
        formula.push(Formula::Or(clause_all_neg));

        Formula::And(formula)
    }

    /// Apply the Tseitin rule for the NOT NodeType.
    fn tseitin_not(&self, self_id: NodeId, arg: NodeId) -> Formula<NodeId> {
        Formula::And(vec![
            Formula::Or(vec![
                Formula::Atom(self_id.clone()),
                Formula::Atom(arg.clone()),
            ]),
            Formula::Or(vec![
                Formula::Not(Box::new(Formula::Atom(self_id.clone()))),
                Formula::Not(Box::new(Formula::Atom(arg.clone()))),
            ]),
        ])
    }

    /// Apply the Tseitin transformation for the Node, depending on the type of node
    /// will use different rules.
    /// Look that the output type is a Formula of NodeIds, ready for be used in the CNF.
    pub fn tseitin_transformation(&self, self_id: NodeId) -> Formula<NodeId> {
        match &self.kind {
            NodeType::PlaceHolder(_, _, _) => {
                panic!("Cant apply Tseitin transform to placeholder node.")
            }
            NodeType::BasicEvent(_, _, _) => Formula::_True,
            NodeType::Not(arg) => self.tseitin_not(self_id, arg.to_owned()),
            NodeType::And(args) => self.tseitin_and(self_id, args),
            NodeType::Or(args) => self.tseitin_or(self_id, args),
            NodeType::Xor(args) => self.tseitin_xor(self_id, args),
            // VOT case is handled when reading the file.
            NodeType::Vot(k, args) => panic!("Unprocessed VOT gate {:?} {:?}", k, args),
        }
    }
}
