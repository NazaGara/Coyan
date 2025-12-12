use crate::formula::Formula;
use itertools::Itertools;
use std::{
    collections::HashMap,
    fmt::{self, Debug, Display},
};

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
#[derive(Debug, Clone)]
pub enum Distribution {
    /// Parameter Lambda of the exponential governing the distribution.
    // TODO: Add usize parameter for phases if Erlang.
    Continuous(f64),
    /// Discrete probability.
    Discrete(f64),
}
#[derive(Debug, Clone)]
pub enum RepairMode {
    /// Parameter is the rate governing the repair time.
    Monitored(f64),
    /// Parameters are: the time interval, and the avg repair time
    PeriodicallyTested(f64, f64),
}

#[derive(Debug, Clone)]
pub struct BasicEvent {
    dist: Distribution,
    repair_mode: Option<RepairMode>,
}

impl Default for BasicEvent {
    fn default() -> Self {
        BasicEvent {
            dist: Distribution::Discrete(0.0),
            repair_mode: None,
        }
    }
}

impl Display for BasicEvent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (&self.dist, &self.repair_mode) {
            (Distribution::Discrete(prob), _) => write!(f, "prob={prob}"),
            (Distribution::Continuous(lambda), None) => write!(f, "lambda={lambda}"),
            (Distribution::Continuous(lambda), Some(RepairMode::Monitored(t_d))) => {
                write!(f, "lambda={lambda} repair={t_d}")
            }
            (Distribution::Continuous(lambda), Some(RepairMode::PeriodicallyTested(t, t_r))) => {
                write!(f, "lambda={lambda} interval={t} repair={t_r}")
            }
        }
    }
}

impl BasicEvent {
    pub fn const_true() -> Self {
        Self {
            dist: Distribution::Discrete(1.0),
            repair_mode: None,
        }
    }

    pub fn const_false() -> Self {
        Self {
            dist: Distribution::Discrete(0.0),
            repair_mode: None,
        }
    }

    pub fn new_with_prob(prob: f64) -> Self {
        Self {
            dist: Distribution::Discrete(prob),
            repair_mode: None,
        }
    }

    pub fn with_repair_mode(&mut self, rm: RepairMode) {
        self.repair_mode = Some(rm);
    }

    pub fn new_with_rate(prob: f64) -> Self {
        Self {
            dist: Distribution::Continuous(prob),
            repair_mode: None,
        }
    }

    pub fn unreliability(&self, timepoint: f64) -> f64 {
        match &self.dist {
            Distribution::Discrete(prob) => *prob,
            Distribution::Continuous(lambda) => 1.0 - (-lambda * timepoint).exp(),
        }
    }

    /// As taken from Table XI-2 of W. E. Vesely, F. F. Goldberg, N. H. Roberts, and D. F. Haasl,
    /// Fault Tree Handbook. U.S. Nuclear Regulatory Commission, 1981.
    pub fn unavailability(&self, timepoint: f64) -> f64 {
        match (&self.dist, &self.repair_mode) {
            (_, None) => self.unreliability(timepoint),
            (Distribution::Discrete(prob), _) => *prob,
            (Distribution::Continuous(lambda), Some(RepairMode::Monitored(l_d))) => {
                (lambda * (1.0 / l_d)) / (1.0 + (lambda * (1.0 / l_d)))
            }
            (Distribution::Continuous(lambda), Some(RepairMode::PeriodicallyTested(t, t_r))) => {
                ((lambda * t) / 2.0) + (lambda * t_r)
            }
        }
    }
}

/// Enum representing the Types of nodes we can have.
#[derive(Debug, Clone)]
pub enum Node<T> {
    BasicEvent(T, BasicEvent),
    Not(NodeId),
    And(Vec<NodeId>),
    Or(Vec<NodeId>),
    Xor(Vec<NodeId>),
    Vot(i64, Vec<NodeId>),
    PlaceHolder(T, String, Vec<T>),
}

/// Implementation of the nodes.
impl<T> Node<T>
where
    T: Clone + Debug,
{
    pub fn unavailability(&self, timepoint: f64) -> Option<f64> {
        match self {
            Node::BasicEvent(_, be) => Some(be.unavailability(timepoint)),
            _ => None,
        }
    }

    pub fn unreliability(&self, timepoint: f64) -> Option<f64> {
        match self {
            Node::BasicEvent(_, be) => Some(be.unreliability(timepoint)),
            _ => None,
        }
    }

    pub fn is_or(&self) -> bool {
        matches!(self, Node::Or(_))
    }

    pub fn gate_type(&self) -> String {
        match self {
            Node::BasicEvent(_, _) => String::from("basic event"),
            Node::Not(_) => String::from("not"),
            Node::And(_) => String::from("and"),
            Node::Or(_) => String::from("or"),
            Node::Vot(_, _) => String::from("vot"),
            Node::Xor(_) => String::from("xor"),
            Node::PlaceHolder(_, _, _) => String::from("placeholder"),
        }
    }

    pub fn map_to_args(&mut self, mapper: &HashMap<NodeId, NodeId>) {
        *self = match &self {
            Node::PlaceHolder(_, _, _) => return,
            Node::BasicEvent(_, _) => return,
            Node::Not(arg) => Node::Not(mapper.get(arg).unwrap().to_owned()),
            Node::And(args) => Node::And(
                args.iter()
                    .map(|a| mapper.get(a).unwrap().to_owned())
                    .collect_vec(),
            ),
            Node::Or(args) => Node::Or(
                args.iter()
                    .map(|a| mapper.get(a).unwrap().to_owned())
                    .collect_vec(),
            ),
            Node::Xor(args) => Node::Xor(
                args.iter()
                    .map(|a| mapper.get(a).unwrap().to_owned())
                    .collect_vec(),
            ),
            Node::Vot(k, args) => Node::Vot(
                k.to_owned(),
                args.iter()
                    .map(|a| mapper.get(a).unwrap().to_owned())
                    .collect_vec(),
            ),
        };
    }

    pub fn is_gate(&self) -> bool {
        match &self {
            Node::PlaceHolder(_, _, _) => false,
            Node::BasicEvent(_, _) => false,
            Node::Not(_) => true,
            Node::And(_) => true,
            Node::Or(_) => true,
            Node::Xor(_) => true,
            Node::Vot(_, _) => true,
        }
    }

    pub fn children(&self) -> Vec<NodeId> {
        match &self {
            Node::PlaceHolder(_, _, _) => {
                panic!("Cant apply Tseitin transform to placeholder node.")
            }
            Node::BasicEvent(_, _) => vec![],
            Node::Not(arg) => vec![*arg],
            Node::And(args) => args.to_vec(),
            Node::Or(args) => args.to_vec(),
            Node::Xor(args) => args.to_vec(),
            Node::Vot(_, args) => args.to_vec(),
        }
    }

    /// Apply the Tseitin rule for the AND NodeType.
    fn tseitin_and(&self, self_id: NodeId, args: &[NodeId]) -> Formula<NodeId> {
        let mut clauses: Vec<Formula<_>> = args
            .iter()
            .map(|nid| {
                Formula::Or(vec![
                    Formula::Not(Box::new(Formula::Atom(self_id))),
                    Formula::Atom(*nid),
                ])
            })
            .collect();
        let mut other: Vec<Formula<NodeId>> = args
            .iter()
            .map(|nid| Formula::Not(Box::new(Formula::Atom(*nid))))
            .collect();
        other.push(Formula::Atom(self_id));
        clauses.push(Formula::Or(other));
        Formula::And(clauses)
    }

    /// Apply the Tseitin rule for the OR NodeType.
    fn tseitin_or(&self, self_id: NodeId, args: &[NodeId]) -> Formula<NodeId> {
        let mut clauses: Vec<Formula<_>> = args
            .iter()
            .map(|nid| {
                Formula::Or(vec![
                    Formula::Atom(self_id),
                    Formula::Not(Box::new(Formula::Atom(*nid))),
                ])
            })
            .collect();
        let mut other: Vec<Formula<NodeId>> = args.iter().map(|nid| Formula::Atom(*nid)).collect();
        other.push(Formula::Not(Box::new(Formula::Atom(self_id))));
        clauses.push(Formula::Or(other));
        Formula::And(clauses)
    }

    /// Apply the Tseitin rule for the XOR NodeType.
    fn tseitin_xor(&self, self_id: NodeId, args: &[NodeId]) -> Formula<NodeId> {
        let mut clause_g_neg = args.iter().map(|nid| Formula::Atom(*nid)).collect_vec();
        clause_g_neg.push(Formula::Not(Box::new(Formula::Atom(self_id))));

        let mut clause_all_neg = args
            .iter()
            .map(|nid| Formula::Not(Box::new(Formula::Atom(*nid))))
            .collect_vec();
        clause_all_neg.push(Formula::Not(Box::new(Formula::Atom(self_id))));

        let mut other_clauses = (0..args.len())
            .map(|_| {
                let mut c = args.iter().map(|nid| Formula::Atom(*nid)).collect_vec();
                c.push(Formula::Atom(self_id));
                c
            })
            .collect_vec();

        for (i, c) in other_clauses.iter_mut().enumerate() {
            let f = c.remove(i).negate();
            c.insert(i, f);
        }

        let mut formula = other_clauses.into_iter().map(Formula::Or).collect_vec();
        formula.push(Formula::Or(clause_g_neg));
        formula.push(Formula::Or(clause_all_neg));

        Formula::And(formula)
    }

    /// Apply the Tseitin rule for the NOT NodeType.
    fn tseitin_not(&self, self_id: NodeId, arg: NodeId) -> Formula<NodeId> {
        Formula::And(vec![
            Formula::Or(vec![Formula::Atom(self_id), Formula::Atom(arg)]),
            Formula::Or(vec![
                Formula::Not(Box::new(Formula::Atom(self_id))),
                Formula::Not(Box::new(Formula::Atom(arg))),
            ]),
        ])
    }

    /// Apply the Tseitin transformation for the Node, depending on the type of node
    /// will use different rules.
    /// The output type is a Formula of NodeIds, ready to be used in the CNF.
    pub fn tseitin_transformation(&self, self_id: NodeId) -> Formula<NodeId> {
        match &self {
            Node::PlaceHolder(_, _, _) => {
                panic!("Cant apply Tseitin transform to placeholder node.")
            }
            Node::BasicEvent(_, _) => Formula::True,
            Node::Not(arg) => self.tseitin_not(self_id, arg.to_owned()),
            Node::And(args) => self.tseitin_and(self_id, args),
            Node::Or(args) => self.tseitin_or(self_id, args),
            Node::Xor(args) => self.tseitin_xor(self_id, args),
            // VOT case is handled when reading the file.
            Node::Vot(k, args) => panic!("Unprocessed VOT gate {:?} {:?}", k, args),
        }
    }
}
