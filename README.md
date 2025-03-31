# Feasibility Jump experimental version 

This is the code for Feasibility Jump that was used as an entry into the MIP
2022 Computational Competition.

## Usage

 * Install a version of Gurobi supported by the `grb` library for Rust (9.5?)
 * Install Rust
 * Build `cargo build --release`
 * Copy the executable `cp target/release/lsmip .`
 * Enter a Python virtual environment with the `gurobipy` package: `python -m venv env; source env/bin/activate; pip install gurobipy`.
 * Run the solver `./solve.sh qap10.mps`.

