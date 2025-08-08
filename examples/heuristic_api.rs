use lsmip::heuristic::*;
use mpsparser::{RowType, VarType};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Simple 2-variable problem: minimize x + y subject to x + y >= 1, x,y integer 0-1

    let variables = vec![
        Variable {
            name: Some("x".to_string()),
            var_type: VarType::Integer,
            lower_bound: 0.0,
            upper_bound: 1.0,
            objective_coeff: 1.0,
        },
        Variable {
            name: Some("y".to_string()),
            var_type: VarType::Integer,
            lower_bound: 0.0,
            upper_bound: 1.0,
            objective_coeff: 1.0,
        },
    ];

    let constraints = vec![Constraint {
        name: Some("sum_constraint".to_string()),
        sense: RowType::Gte,
        rhs: 1.0,
        coefficients: vec![(0, 1.0), (1, 1.0)],
    }];

    let heuristic = HeuristicBuilder::new()
        .time_limit(5.0)
        .threads(4) // Test multi-threading
        .seed(12345) // Test seed consistency
        .weight_decay(0.8) // Test weight decay parameter
        .build();

    let solution = heuristic.solve(&variables, &constraints);

    println!("Solution status: {:?}", solution.status);

    if solution.is_feasible() {
        println!("Objective value: {}", solution.objective_value);
        println!("Variable values: {:?}", solution.variable_values);
    } else {
        println!("No feasible solution found");
    }

    println!("Solve time: {:.3}s", solution.solve_time);
    println!("Iterations: {}", solution.iterations);

    Ok(())
}
