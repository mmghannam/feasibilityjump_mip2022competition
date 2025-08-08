mod interval;
#[cfg(feature = "gurobi")]
pub mod lp;
pub mod solver;

// Public heuristic API
pub mod heuristic;

// Re-export the main heuristic types for convenience
pub use heuristic::{
    Constraint, FeasibilityJumpHeuristic, HeuristicBuilder, HeuristicConfig, Solution,
    SolutionStatus, Variable,
};

// Re-export variable and constraint types from mpsparser for convenience
pub use mpsparser::{RowType, VarType};

pub const UNBOUNDED: f64 = 1.0e5;

pub const USE_EXPONENTIAL_WEIGHT_DECAY: bool = true;
pub const USE_DAMPING: bool = false;
pub const USE_SIZE_WEIGHTING: bool = false;
