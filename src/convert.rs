
use ordered_float::OrderedFloat;

#[derive(Debug)]
#[derive(Clone)]
pub struct FloatConstraint {
    pub rhs: f64,
    pub terms: Vec<(usize, f64)>,
}

#[derive(Debug)]
pub struct IntConstraint {
    pub rhs: i32,
    pub terms: Vec<(usize, i32)>,
}

pub fn round_constraint(c: &FloatConstraint) -> IntConstraint {
    IntConstraint {
        rhs: c.rhs.round() as i32,
        terms: c
            .terms
            .iter()
            .map(|(v, c)| (*v, c.round() as _))
            .filter(|(v, c)| *c != 0)
            .collect(),
    }
}

pub fn invert_constraint(c: &FloatConstraint) -> FloatConstraint {
    FloatConstraint {
        rhs: -c.rhs,
        terms: c.terms.iter().map(|(v, c)| (*v, -*c)).collect(),
    }
}

pub fn round_err(x: f64) -> f64 {
    (x - x.round()).abs() / (x.abs().round() + 1.)
}

pub fn constraint_round_err(c: &FloatConstraint) -> f64 {
    std::iter::once(round_err(c.rhs))
        .chain(c.terms.iter().map(|(_, c)| round_err(*c)))
        .map(OrderedFloat)
        .max()
        .map(OrderedFloat::into_inner)
        .unwrap_or(0.)
}

pub fn rescale(c: FloatConstraint) -> FloatConstraint {
    let smallest = c
        .terms
        .iter()
        .map(|(_, c)| OrderedFloat(c.abs()))
        .min()
        .unwrap()
        .0;
    if smallest > 1.0e-4 && smallest < 1.0 - 1e-4 {
        let c = FloatConstraint {
            rhs: c.rhs / smallest,
            terms: c.terms.iter().map(|(v, c)| (*v, *c / smallest)).collect(),
        };
        c
    } else {
        c
    }
}

pub fn get_constraints(mps: &mpsparser::MPSInstance) -> Vec<(i32, Vec<(usize, i32)>)> {
    let mut constraints = Vec::new();
    for constraint in mps.constraints.iter() {
        // println!("{:?}", constraint);
        match constraint.rowtype {
            mpsparser::RowType::None => continue,
            mpsparser::RowType::Gte => {
                constraints.push((
                    constraint
                        .rhs
                        .map(|r| -r.as_f64().round() as i32)
                        .unwrap_or(0),
                    constraint
                        .cells
                        .iter()
                        .map(|c| (c.var, -c.coeff.as_f64().round() as i32))
                        .collect(),
                ));
            }
            mpsparser::RowType::Lte => {
                constraints.push((
                    constraint
                        .rhs
                        .map(|r| r.as_f64().round() as i32)
                        .unwrap_or(0),
                    constraint
                        .cells
                        .iter()
                        .map(|c| (c.var, c.coeff.as_f64().round() as i32))
                        .collect(),
                ));
            }
            mpsparser::RowType::Equal => {
                constraints.push((
                    constraint
                        .rhs
                        .map(|r| r.as_f64().round() as i32)
                        .unwrap_or(0),
                    constraint
                        .cells
                        .iter()
                        .map(|c| (c.var, c.coeff.as_f64().round() as i32))
                        .collect(),
                ));
                constraints.push((
                    constraint
                        .rhs
                        .map(|r| -r.as_f64().round() as i32)
                        .unwrap_or(0),
                    constraint
                        .cells
                        .iter()
                        .map(|c| (c.var, -c.coeff.as_f64().round() as i32))
                        .collect(),
                ));
            }
        }
    }
    constraints
}
