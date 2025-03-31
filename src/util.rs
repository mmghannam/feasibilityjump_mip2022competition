use log::info;
use mpsparser::Number;

use crate::UNBOUNDED;

pub fn check_domains(mps: &mpsparser::MPSInstance) {
    let mut n_continuous = 0;
    let mut n_integer = 0;
    let mut n_binary = 0;
    let mut n_binary_convert = 0;
    let mut integer_domain_size = 0;
    let mut integer_unbounded = 0;
    for var in mps.variables.iter() {
        match var.var_type {
            mpsparser::VarType::Continuous => {
                n_continuous += 1;
            }
            mpsparser::VarType::Integer => {
                let int_lb = var.lb.map(Number::as_f64).unwrap_or(0.).ceil() as i32;
                let int_ub = var.ub.map(Number::as_f64).unwrap_or(UNBOUNDED).floor() as i32;
                let domain_size = int_ub - int_lb;
                if domain_size == 1 {
                    if int_lb == 0 {
                        assert!(int_ub == 1);
                        n_binary += 1;
                    } else {
                        n_binary_convert += 1;
                    }
                } else if domain_size >= UNBOUNDED as _ {
                    n_integer += 1;
                    integer_unbounded += 1;
                } else {
                    n_integer += 1;
                    integer_domain_size += domain_size;
                }
            }
        }
    }
    info!(
        "  cont:{} int:{} bin:{} bin_conv:{} domain_size:{} unbounded:{}",
        n_continuous, n_integer, n_binary, n_binary_convert, integer_domain_size, integer_unbounded
    );
    // assert!(n_binary_convert == 0);
}

pub fn check_constraints(mps: &mpsparser::MPSInstance) {
    let constraints = crate::convert::get_constraints(mps);
    let mut n_clause = 0;
    let mut n_atmost1 = 0;
    let mut n_card = 0;
    let mut n_pb = 0;
    let mut n_general = 0;
    for constraint in &constraints {
        if constraint_is_clause(mps, constraint) {
            n_clause += 1;
        } else if constraint_is_atmost1(mps, constraint) {
            n_atmost1 += 1;
        } else if constraint_is_card(mps, constraint) {
            n_card += 1;
        } else if constraint_is_pb(mps, constraint) {
            n_pb += 1;
        } else {
            n_general += 1;
            // println!("general {:?}", constraint);
        }
    }
    info!(
        "  clauses:{} atm1:{} card:{} pb:{} general:{}",
        n_clause, n_atmost1, n_card, n_pb, n_general
    );
}


fn constraint_is_pb(mps: &mpsparser::MPSInstance, constraint: &(i32, Vec<(usize, i32)>)) -> bool {
    for (v, _) in constraint.1.iter() {
        let var = &mps.variables[*v];
        let int_lb = var.lb.map(Number::as_f64).unwrap_or(0.).ceil() as i32;
        let int_ub = var.ub.map(Number::as_f64).unwrap_or(UNBOUNDED).floor() as i32;

        if int_lb != 0 || int_ub != 1 {
            return false;
        }
    }
    true
}

fn constraint_is_card(
    mps: &mpsparser::MPSInstance,
    (mut rhs, terms): &(i32, Vec<(usize, i32)>),
) -> bool {
    for (var, coeff) in terms.iter() {
        let var = &mps.variables[*var];
        let int_lb = var.lb.map(Number::as_f64).unwrap_or(0.).ceil() as i32;
        let int_ub = var.ub.map(Number::as_f64).unwrap_or(UNBOUNDED).floor() as i32;

        if int_lb != 0 || int_ub != 1 || coeff.abs() != 1 {
            return false;
        }

        if *coeff == -1 {
            rhs -= 1;
        }
    }

    rhs >= 1
}

fn constraint_is_atmost1(
    mps: &mpsparser::MPSInstance,
    (mut rhs, terms): &(i32, Vec<(usize, i32)>),
) -> bool {
    // at-most-one constraints are x1 + x2 + x3 <= 1
    // if some are negated, we have  x1 + x2 + (1 - x3) + (1 - x4) <= 1
    // x1 + x2 - x3 - x4 <= 1 - k
    // -x1 - x2 + x3 + x4 >= k - 1

    for (var, coeff) in terms.iter() {
        let var = &mps.variables[*var];
        let int_lb = var.lb.map(Number::as_f64).unwrap_or(0.).ceil() as i32;
        let int_ub = var.ub.map(Number::as_f64).unwrap_or(UNBOUNDED).floor() as i32;

        if int_lb != 0 || int_ub != 1 || coeff.abs() != 1 {
            return false;
        }

        if *coeff == -1 {
            rhs -= 1;
        }
    }

    rhs == 1
}

fn constraint_is_clause(
    mps: &mpsparser::MPSInstance,
    (mut rhs, terms): &(i32, Vec<(usize, i32)>),
) -> bool {
    // x + y >= 1
    // -x -y <= -1

    for (var, coeff) in terms.iter() {
        let var = &mps.variables[*var];
        let int_lb = var.lb.map(Number::as_f64).unwrap_or(0.).ceil() as i32;
        let int_ub = var.ub.map(Number::as_f64).unwrap_or(UNBOUNDED).floor() as i32;

        if int_lb != 0 || int_ub != 1 || coeff.abs() != 1 {
            return false;
        }

        if *coeff == 1 {
            rhs -= 1;
        }
    }

    rhs == -1
}
