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
// #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Debug, Clone, PartialEq)]
pub enum Formula<A> {
    Atom(A),
    Not(Box<Self>),
    And(Vec<Self>),
    Or(Vec<Self>),
    Vot(i64, Vec<Self>),
    _True,
    _False,
}

impl<A> Formula<A>
where
    A: Clone + Display + Debug + PartialEq,
{
    pub fn negate(&self) -> Formula<A> {
        match self {
            Formula::Atom(a) => Formula::Not(Box::new(Formula::Atom(a.clone()))),
            Formula::Not(b) => *b.to_owned(),
            Formula::_False => Formula::_True,
            Formula::_True => Formula::_False,
            Formula::And(args) => Formula::Or(args.iter().map(|a| a.negate()).collect_vec()),
            Formula::Or(args) => Formula::And(args.iter().map(|a| a.negate()).collect_vec()),
            Formula::Vot(_k, _args) => todo!("Negations of VOT is not allowed at the moment."),
        }
    }
    /// Translate the formula to a GALILEO style of formula
    pub fn _formula_to_dft(&self) -> String {
        let mut out = "".to_owned();
        match self {
            Formula::Atom(name) => out += &name.to_string(),
            Formula::Vot(k, args) => {
                let text = format!("{}of{} ", k, args.len());
                let rec = args.iter().map(|a| a._formula_to_dft()).join(" ");
                out += &text;
                out += &rec;
            }
            Formula::And(args) => {
                let rec = args.iter().map(|a| a._formula_to_dft()).join(" ");
                out.push_str("and ");
                out += &rec;
            }
            Formula::Or(args) => {
                let rec = args.iter().map(|a| a._formula_to_dft()).join(" ");
                out.push_str("or ");
                out += &rec;
            }
            Formula::Not(x) => {
                let inner = *x.to_owned();
                out.push_str("-");
                out += &inner._formula_to_dft();
            }
            _ => {}
        }
        out
    }

    /// Method that reduces the size of the formula a bit.
    pub fn _reduce_formula(&self) -> Formula<A> {
        match self {
            Formula::And(args) => Formula::And(
                args.iter()
                    .flat_map(|a| match a {
                        Formula::And(inner_args) => inner_args.to_owned(),
                        f => vec![f.to_owned()],
                    })
                    .collect_vec(),
            ),
            Formula::Or(args) => Formula::Or(
                args.iter()
                    .flat_map(|a| match a {
                        Formula::Or(inner_args) => inner_args.to_owned(),
                        f => vec![f.to_owned()],
                    })
                    .collect_vec(),
            ),
            _ => self.to_owned(),
        }
    }

    /// Returns the number of clauses of a CNF formula.
    pub fn num_clauses(&self) -> usize {
        match self {
            Formula::And(clauses) => clauses.len(),
            _ => {
                panic!("Formula is not a AND Formula, therefore is not a CNF formula.");
            }
        }
    }

    pub fn _get_args(&self) -> Option<Vec<Self>> {
        match self {
            Formula::And(args) => Some(args.to_vec()),
            Formula::Or(args) => Some(args.to_vec()),
            Formula::Vot(_, args) => Some(args.to_vec()),
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
        out.push_str(")");
        out
    }
}
