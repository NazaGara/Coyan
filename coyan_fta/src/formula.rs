use itertools::Itertools;
use std::{
    fmt::{Debug, Display},
    str::FromStr,
};

#[derive(Clone, Copy)]
pub enum CNFFormat {
    MCC,
    MC21,
}
impl FromStr for CNFFormat {
    type Err = ();

    fn from_str(input: &str) -> Result<CNFFormat, Self::Err> {
        match input {
            "MCC" => Ok(CNFFormat::MCC),
            "mcc" => Ok(CNFFormat::MCC),
            "MC21" => Ok(CNFFormat::MC21),
            "mc21" => Ok(CNFFormat::MC21),
            _ => Err(()),
        }
    }
}

/// A propositional logic formula.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Formula<A> {
    Atom(A),
    Not(Box<Self>),
    And(Vec<Self>),
    Or(Vec<Self>),
    Xor(Vec<Self>),
    Vot(i64, Vec<Self>),
    True,
    False,
}

impl<A> Formula<A>
where
    A: Clone + Display + Debug + PartialEq,
{
    pub fn negate(&self) -> Formula<A> {
        match self {
            Formula::Atom(a) => Formula::Not(Box::new(Formula::Atom(a.clone()))),
            Formula::Not(b) => *b.to_owned(),
            Formula::False => Formula::True,
            Formula::True => Formula::False,
            Formula::And(args) => Formula::Or(args.iter().map(|a| a.negate()).collect_vec()),
            Formula::Or(args) => Formula::And(args.iter().map(|a| a.negate()).collect_vec()),
            Formula::Xor(args) => {
                // All arguments "as is"
                let pos_args = Formula::And(args.clone());
                // All arguments negated
                let neg_args = Formula::And(args.iter().map(|a| a.negate()).collect_vec());
                Formula::Or(vec![pos_args, neg_args])
            }
            Formula::Vot(_k, _args) => panic!("Negation of VOT gates is not allowed."),
        }
    }

    /// Translate the formula to a GALILEO style of formula
    pub fn formula_to_dft(&self) -> String {
        let mut out = "".to_owned();
        match self {
            Formula::Atom(name) => out += &name.to_string(),
            Formula::Vot(k, args) => {
                let text = format!("{}of{} ", k, args.len());
                let rec = args.iter().map(|a| a.formula_to_dft()).join(" ");
                out += &text;
                out += &rec;
            }
            Formula::And(args) => {
                let rec = args.iter().map(|a| a.formula_to_dft()).join(" ");
                out.push_str("and ");
                out += &rec;
            }
            Formula::Or(args) => {
                let rec = args.iter().map(|a| a.formula_to_dft()).join(" ");
                out.push_str("or ");
                out += &rec;
            }
            Formula::Xor(args) => {
                let rec = args.iter().map(|a| a.formula_to_dft()).join(" ");
                out.push_str("xor ");
                out += &rec;
            }
            Formula::Not(x) => {
                let inner = *x.to_owned();
                out.push('-');
                out += &inner.formula_to_dft();
            }
            _ => {}
        }
        out
    }

    /// Returns the number of clauses only for CNF formulas.
    pub fn num_clauses(&self) -> Option<usize> {
        match self {
            Formula::And(clauses) => Some(clauses.len()),
            _ => None,
        }
    }

    /// Translate the Formula to text.
    pub fn to_text(&self) -> String {
        let mut sep: String = String::from("");
        let mut out: String = String::from("(");
        let mut arguments: Vec<Formula<A>> = vec![];
        match self {
            Formula::Vot(k, args) => {
                sep.push_str("VOT(");
                sep.push_str(&k.to_string());
                sep.push_str(") ");
                arguments = args.to_vec();
            }
            Formula::And(args) => {
                arguments = args.to_vec();
                sep = String::from(" ∧ ")
            }
            Formula::Or(args) => {
                arguments = args.to_vec();
                sep = String::from(" V ")
            }
            Formula::Xor(args) => {
                arguments = args.to_vec();
                sep = String::from(" ⊕ ")
            }
            Formula::Not(x) => {
                let mut text = String::from("-");
                text += &x.to_text();
                return text;
            }
            Formula::Atom(x) => return x.to_string(),
            _ => {}
        }
        out += &arguments
            .iter()
            .map(|a| a.to_text())
            .collect::<Vec<String>>()
            .join(&sep);
        out.push(')');
        out
    }
}
