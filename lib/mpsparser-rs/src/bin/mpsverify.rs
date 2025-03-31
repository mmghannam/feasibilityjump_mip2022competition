use std::collections::HashMap;

use mpsparser::{check_values, parse_solution, DEFAULT_EQ_TOLERANCE, DEFAULT_INT_TOLERANCE};

pub fn main() -> Result<(), &'static str> {
    let args = std::env::args().collect::<Vec<_>>();
    if args.len() != 3 {
        return Err("Usage: mpsverify <MPSFILE> <SOLFILE>");
    }
    let mps_filename = std::path::Path::new(&args[1]);
    let sol_filename = std::path::Path::new(&args[2]);

    let mps = mpsparser::parse(
        std::fs::read(mps_filename)
            .map_err(|_| "Could not read MPS file.")?
            .as_slice(),
    )
    .map_err(|_| "Could not parse MPS file.")?;
    let sol_contents =
        std::fs::read_to_string(sol_filename).map_err(|_| "Could not read SOL file.")?;
    let sol = parse_solution(&sol_contents)
        .map_err(|_| "Could not read SOL file.")?
        .into_iter()
        .collect::<HashMap<_, _>>();

    let mut var_values = Vec::new();
    for (idx,v) in mps.variables.iter().enumerate() {
        if let Some(x) = sol.get(&v.name.as_str()) {
            var_values.push(*x as f64);
        } else {
            println!("No value for variable idx={} name={}", idx, v.name);
            panic!();
        }
    }

    check_values(
        &mps,
        &var_values,
        DEFAULT_INT_TOLERANCE,
        DEFAULT_EQ_TOLERANCE,
    )?;
    let objective = match mps.objective() {
        Some((constant, obj)) => {
            let constant_term = constant.map(|n| n.as_f64()).unwrap_or(0.);
            -constant_term
                + obj
                    .iter()
                    .map(|c| c.coeff.as_f64() * var_values[c.var])
                    .sum::<f64>()
        }
        None => 0.,
    } as f32;
    println!("Objective value: {}", objective);
    Ok(())
}
