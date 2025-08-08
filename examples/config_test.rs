use lsmip::heuristic::*;
use mpsparser::{RowType, VarType};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // More complex problem to test configuration parameters
    println!("Testing Feasibility Jump Heuristic Configuration");

    let variables = vec![
        Variable {
            name: Some("x1".to_string()),
            var_type: VarType::Integer,
            lower_bound: 0.0,
            upper_bound: 10.0,
            objective_coeff: 1.0,
        },
        Variable {
            name: Some("x2".to_string()),
            var_type: VarType::Integer,
            lower_bound: 0.0,
            upper_bound: 8.0,
            objective_coeff: 2.0,
        },
        Variable {
            name: Some("x3".to_string()),
            var_type: VarType::Continuous,
            lower_bound: 0.0,
            upper_bound: 5.0,
            objective_coeff: 1.5,
        },
    ];

    let constraints = vec![
        Constraint {
            name: Some("constraint1".to_string()),
            sense: RowType::Lte,
            rhs: 15.0,
            coefficients: vec![(0, 2.0), (1, 1.0), (2, 0.5)],
        },
        Constraint {
            name: Some("constraint2".to_string()),
            sense: RowType::Gte,
            rhs: 3.0,
            coefficients: vec![(0, 1.0), (2, 1.0)],
        },
    ];

    // Test 1: Single thread with specific configuration
    println!("\n=== Test 1: Single Thread ===");
    let heuristic1 = HeuristicBuilder::new()
        .time_limit(2.0)
        .threads(1)
        .seed(42)
        .weight_decay(0.9)
        .relax_continuous(true)
        .build();

    let start = std::time::Instant::now();
    let solution1 = heuristic1.solve(&variables, &constraints);
    let elapsed1 = start.elapsed();

    println!("Config: 1 thread, 2s limit, seed=42, decay=0.9");
    println!("Status: {:?}", solution1.status);
    println!("Time: {:.3}s", elapsed1.as_secs_f64());
    if solution1.is_feasible() {
        println!("Objective: {}", solution1.objective_value);
        println!("Variables: {:?}", solution1.variable_values);
    }
    println!("Iterations: {}", solution1.iterations);

    // Test 2: Multi-thread with different configuration
    println!("\n=== Test 2: Multi-Thread ===");
    let heuristic2 = HeuristicBuilder::new()
        .time_limit(3.0)
        .threads(4)
        .seed(123)
        .weight_decay(0.7)
        .relax_continuous(false)
        .build();

    let start = std::time::Instant::now();
    let solution2 = heuristic2.solve(&variables, &constraints);
    let elapsed2 = start.elapsed();

    println!("Config: 4 threads, 3s limit, seed=123, decay=0.7");
    println!("Status: {:?}", solution2.status);
    println!("Time: {:.3}s", elapsed2.as_secs_f64());
    if solution2.is_feasible() {
        println!("Objective: {}", solution2.objective_value);
        println!("Variables: {:?}", solution2.variable_values);
    }
    println!("Iterations: {}", solution2.iterations);

    // Test 3: Deterministic behavior with same seed
    println!("\n=== Test 3: Deterministic Behavior ===");
    let heuristic3a = HeuristicBuilder::new()
        .time_limit(1.0)
        .threads(1)
        .seed(999)
        .build();

    let heuristic3b = HeuristicBuilder::new()
        .time_limit(1.0)
        .threads(1)
        .seed(999)
        .build();

    let solution3a = heuristic3a.solve(&variables, &constraints);
    let solution3b = heuristic3b.solve(&variables, &constraints);

    println!("Same seed (999) should give similar results:");
    println!(
        "Run A - Objective: {:.6}, Status: {:?}",
        solution3a.objective_value, solution3a.status
    );
    println!(
        "Run B - Objective: {:.6}, Status: {:?}",
        solution3b.objective_value, solution3b.status
    );

    if solution3a.is_feasible() && solution3b.is_feasible() {
        let obj_diff = (solution3a.objective_value - solution3b.objective_value).abs();
        println!("Objective difference: {:.6}", obj_diff);
        if obj_diff < 1e-6 {
            println!("✓ Deterministic behavior confirmed!");
        } else {
            println!("⚠ Results differ - might indicate thread-specific behavior");
        }
    }

    Ok(())
}
