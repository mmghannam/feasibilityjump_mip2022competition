// Public API for using Feasibility Jump as a MIP heuristic

use mpsparser::{RowType, VarType};

/// Configuration for the Feasibility Jump heuristic
#[derive(Debug, Clone)]
pub struct HeuristicConfig {
    /// Number of threads to use (default: 1 for embedded use)
    pub threads: usize,
    /// Time limit in seconds (default: 30.0)
    pub time_limit: f64,
    /// Random seed (default: None for random)
    pub seed: Option<u64>,
    /// Weight decay factor (default: 1.0)
    pub weight_decay: f64,
    /// Whether to relax continuous variables (default: true)
    pub relax_continuous: bool,
    /// Iteration limit
    pub iteration_limit: Option<usize>,
}

impl Default for HeuristicConfig {
    fn default() -> Self {
        Self {
            threads: 1, // Conservative for embedded use
            time_limit: 30.0,
            seed: None,
            weight_decay: 1.0,
            relax_continuous: true,
            iteration_limit: None, // No limit by default
        }
    }
}

/// A variable in the optimization problem
#[derive(Debug, Clone)]
pub struct Variable {
    pub var_type: VarType,
    pub lower_bound: f64,
    pub upper_bound: f64,
    pub objective_coeff: f64,
    pub name: Option<String>,
}

/// A constraint in the optimization problem
#[derive(Debug, Clone)]
pub struct Constraint {
    pub sense: RowType,
    pub rhs: f64,
    pub coefficients: Vec<(usize, f64)>, // (variable_index, coefficient)
    pub name: Option<String>,
}

/// Solution status from the heuristic
#[derive(Debug, Clone, PartialEq)]
pub enum SolutionStatus {
    /// Found a feasible solution
    Feasible,
    /// No feasible solution found within time limit
    NoSolutionFound,
    /// Heuristic was interrupted or stopped early
    Interrupted,
    /// Time limit reached (may have found a solution)
    TimeLimitReached,
}

/// Solution returned by the heuristic
#[derive(Debug, Clone)]
pub struct Solution {
    pub status: SolutionStatus,
    pub objective_value: f64,
    pub variable_values: Vec<f64>,
    pub solve_time: f64,
    pub iterations: usize,
}

impl Solution {
    /// Check if a solution was found
    pub fn is_feasible(&self) -> bool {
        matches!(
            self.status,
            SolutionStatus::Feasible | SolutionStatus::TimeLimitReached
        ) && self.objective_value.is_finite()
    }

    /// Get the value of a specific variable
    pub fn get_variable_value(&self, var_index: usize) -> Option<f64> {
        self.variable_values.get(var_index).copied()
    }
}

/// Main heuristic solver for embedding in MIP solvers
pub struct FeasibilityJumpHeuristic {
    config: HeuristicConfig,
}

impl FeasibilityJumpHeuristic {
    /// Create a new heuristic with default configuration
    pub fn new() -> Self {
        Self {
            config: HeuristicConfig::default(),
        }
    }

    /// Create a new heuristic with custom configuration
    pub fn with_config(config: HeuristicConfig) -> Self {
        Self { config }
    }

    /// Solve the given problem and return the best solution found
    pub fn solve(&self, variables: &[Variable], constraints: &[Constraint]) -> Solution {
        self.solve_with_initial_solution(variables, constraints, None)
    }

    /// Solve with an initial solution (warm start)
    pub fn solve_with_initial_solution(
        &self,
        variables: &[Variable],
        constraints: &[Constraint],
        initial_solution: Option<&[f64]>,
    ) -> Solution {
        use std::time::Instant;
        let start_time = Instant::now();

        // If only one thread, solve directly without crossbeam overhead
        if self.config.threads == 1 {
            self.solve_single_thread(variables, constraints, initial_solution, start_time)
        } else {
            self.solve_multi_thread(variables, constraints, initial_solution, start_time)
        }
    }

    /// Solve with a callback for intermediate solutions
    pub fn solve_with_callback<F>(
        &self,
        variables: &[Variable],
        constraints: &[Constraint],
        mut callback: F,
    ) -> Solution
    where
        F: FnMut(&Solution) -> bool, // Return false to stop early
    {
        // Implementation would call the callback with intermediate solutions
        // For now, just solve normally
        let solution = self.solve(variables, constraints);
        callback(&solution);
        solution
    }

    /// Try to improve an existing solution
    pub fn improve_solution(
        &self,
        variables: &[Variable],
        constraints: &[Constraint],
        current_solution: &[f64],
        improvement_time_limit: f64,
    ) -> Solution {
        let mut config = self.config.clone();
        config.time_limit = improvement_time_limit;
        let heuristic = Self::with_config(config);
        heuristic.solve_with_initial_solution(variables, constraints, Some(current_solution))
    }
}

impl Default for FeasibilityJumpHeuristic {
    fn default() -> Self {
        Self::new()
    }
}

// Implementation methods (to be implemented)
impl FeasibilityJumpHeuristic {
    fn solve_single_thread(
        &self,
        variables: &[Variable],
        constraints: &[Constraint],
        initial_solution: Option<&[f64]>,
        start_time: std::time::Instant,
    ) -> Solution {
        Self::solve_single_thread_with_id(
            &self.config,
            0, // Thread ID 0 for single-threaded execution
            variables,
            constraints,
            initial_solution,
            start_time,
        )
    }

    fn solve_single_thread_with_id(
        config: &HeuristicConfig,
        thread_id: usize,
        variables: &[Variable],
        constraints: &[Constraint],
        initial_solution: Option<&[f64]>,
        start_time: std::time::Instant,
    ) -> Solution {
        use crate::solver::Solver;

        // Create the internal solver with proper thread ID
        let seed = config.seed.unwrap_or_else(|| {
            use std::time::{SystemTime, UNIX_EPOCH};
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_secs()
                + thread_id as u64
        });

        let mut solver = Solver::with_seed(thread_id, seed as u8, config.weight_decay, config.iteration_limit.unwrap_or(usize::MAX));

        // Add variables to solver
        for var in variables {
            solver.add_var(
                var.var_type,
                var.lower_bound as f32,
                var.upper_bound as f32,
                var.objective_coeff as f32,
            );
        }

        // Set initial values if provided
        if let Some(init_sol) = initial_solution {
            for (var_idx, &value) in init_sol.iter().enumerate() {
                if var_idx < variables.len() {
                    solver.set_initial_value(var_idx, value as f32);
                }
            }
        }

        // Add constraints to solver
        for constraint in constraints {
            if matches!(constraint.sense, RowType::None) {
                continue;
            }

            let coeffs: Vec<(usize, f32)> = constraint
                .coefficients
                .iter()
                .map(|(var_idx, coeff)| (*var_idx, *coeff as f32))
                .collect();

            solver.add_constraint(
                constraint.sense,
                constraint.rhs as f32,
                &coeffs,
                config.relax_continuous,
            );
        }

        // Solve with a time limit
        let mut best_solution: Option<Vec<f64>> = None;
        let mut best_objective = f64::INFINITY;
        let mut iterations = 0;

        let time_limit = std::time::Duration::from_secs_f64(config.time_limit);

        solver.solve(|solution| {
            iterations += 1;

            // Check time limit
            if start_time.elapsed() >= time_limit {
                return; // Stop solving
            }

            // Calculate objective
            let objective: f64 = variables
                .iter()
                .enumerate()
                .map(|(i, var)| var.objective_coeff * solution[i] as f64)
                .sum();

            if objective < best_objective {
                best_objective = objective;
                best_solution = Some(solution.iter().map(|&x| x as f64).collect());
            }
        });

        let solve_time = start_time.elapsed().as_secs_f64();
        let time_limit_reached =
            start_time.elapsed() >= std::time::Duration::from_secs_f64(config.time_limit);

        match best_solution {
            Some(solution) => {
                let status = if time_limit_reached {
                    SolutionStatus::TimeLimitReached
                } else {
                    SolutionStatus::Feasible
                };

                Solution {
                    status,
                    objective_value: best_objective,
                    variable_values: solution,
                    solve_time,
                    iterations,
                }
            }
            None => Solution {
                status: SolutionStatus::NoSolutionFound,
                objective_value: f64::INFINITY,
                variable_values: vec![0.0; variables.len()],
                solve_time,
                iterations,
            },
        }
    }

    fn solve_multi_thread(
        &self,
        variables: &[Variable],
        constraints: &[Constraint],
        initial_solution: Option<&[f64]>,
        start_time: std::time::Instant,
    ) -> Solution {
        use std::sync::{mpsc, Arc};
        use std::thread;

        let (tx, rx) = mpsc::channel();
        let variables = Arc::new(variables.to_vec());
        let constraints = Arc::new(constraints.to_vec());
        let initial_solution = initial_solution.map(|s| s.to_vec());

        let mut handles = Vec::new();

        for thread_id in 0..self.config.threads {
            let tx = tx.clone();
            let variables = Arc::clone(&variables);
            let constraints = Arc::clone(&constraints);
            let initial_solution = initial_solution.clone();
            let mut config = self.config.clone();

            // Vary the seed per thread
            if let Some(base_seed) = config.seed {
                config.seed = Some(base_seed + thread_id as u64);
            } else {
                // Generate different seeds for each thread
                use std::time::{SystemTime, UNIX_EPOCH};
                let base_seed = SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_secs();
                config.seed = Some(base_seed + thread_id as u64);
            }

            let handle = thread::spawn(move || {
                // Create solver directly with proper thread ID
                let solution = Self::solve_single_thread_with_id(
                    &config,
                    thread_id,
                    &variables,
                    &constraints,
                    initial_solution.as_deref(),
                    start_time,
                );
                let _ = tx.send(solution);
            });

            handles.push(handle);
        }

        drop(tx); // Close the sender

        // Collect results from all threads
        let mut best_solution = Solution {
            status: SolutionStatus::NoSolutionFound,
            objective_value: f64::INFINITY,
            variable_values: vec![0.0; variables.len()],
            solve_time: 0.0,
            iterations: 0,
        };

        let mut total_iterations = 0;

        for solution in rx {
            total_iterations += solution.iterations;
            if solution.is_feasible() && solution.objective_value < best_solution.objective_value {
                best_solution = solution;
            }
        }

        // Wait for all threads to complete
        for handle in handles {
            let _ = handle.join();
        }

        best_solution.solve_time = start_time.elapsed().as_secs_f64();
        best_solution.iterations = total_iterations; // Sum of all thread iterations
        best_solution
    }
}

/// Builder for easy heuristic configuration
pub struct HeuristicBuilder {
    config: HeuristicConfig,
}

impl HeuristicBuilder {
    pub fn new() -> Self {
        Self {
            config: HeuristicConfig::default(),
        }
    }

    pub fn threads(mut self, n: usize) -> Self {
        self.config.threads = n;
        self
    }

    pub fn time_limit(mut self, seconds: f64) -> Self {
        self.config.time_limit = seconds;
        self
    }

    pub fn seed(mut self, seed: u64) -> Self {
        self.config.seed = Some(seed);
        self
    }

    pub fn weight_decay(mut self, decay: f64) -> Self {
        self.config.weight_decay = decay;
        self
    }

    pub fn relax_continuous(mut self, relax: bool) -> Self {
        self.config.relax_continuous = relax;
        self
    }

    pub fn build(self) -> FeasibilityJumpHeuristic {
        FeasibilityJumpHeuristic::with_config(self.config)
    }
}

impl Default for HeuristicBuilder {
    fn default() -> Self {
        Self::new()
    }
}
