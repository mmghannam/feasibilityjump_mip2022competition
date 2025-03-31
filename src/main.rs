use std::{collections::HashMap, path::PathBuf, sync::RwLock};

use log::info;
use lsmip::{UNBOUNDED, USE_EXPONENTIAL_WEIGHT_DECAY};
use mpsparser::{MPSInstance, Number};
use structopt::StructOpt;

fn presolve(
    path: &std::path::Path,
    python_path: &Option<String>,
    convert_script_path: &str,
    temp_path: &str,
) -> (std::path::PathBuf, MPSInstance) {
    println!("instance:{}", path.file_name().unwrap().to_string_lossy());

    let mut filename = path.file_name().unwrap().to_string_lossy().to_string();
    filename.push_str("_presolved.mps");
    let mut other_path = PathBuf::from(temp_path);
    other_path.push(filename);

    info!("Presolving");
    let python_path = if let Some(python_path) = python_path {
        PathBuf::from(python_path)
    } else {
        which::which("python").unwrap()
    };

    let output = std::process::Command::new(python_path)
        .args([
            convert_script_path,
            "--presolve",
            &path.to_string_lossy(),
            &other_path.to_string_lossy(),
        ])
        .output()
        .expect("failed to execute process");
    println!("python output: {:?}", output);
    info!("Presolved. File exists: {}", other_path.exists());

    let mps_content = std::fs::read(&other_path).unwrap();
    (
        other_path,
        mpsparser::parse(mps_content.as_slice()).unwrap(),
    )
}

#[derive(Debug, StructOpt)]
#[structopt()]
struct Opts {
    #[structopt(name = "INPUT")]
    input_file: String,

    #[structopt(name = "OUTPUT")]
    output_file: String,

    #[structopt(long)]
    python_path: Option<String>,

    #[structopt(long)]
    convert_script_path: String,

    #[structopt(long)]
    temp_path: String,
}

fn main() {
    pretty_env_logger::env_logger::Builder::from_env(
        pretty_env_logger::env_logger::Env::default().default_filter_or("trace"),
    )
    .init();

    let _p = hprof::enter("solver");
    // let opt = Opts::from_args();
    // println!("{:?}", opt);

    mpsparser::cli::with_cli_mpsinstance_path(|path,mps,output| {

        // let path = std::path::Path::new(&opt.input_file);

        println!("instance:{}", path.file_name().unwrap().to_string_lossy());
        let original_mps_obj = {
            let mps_content = std::fs::read(path).unwrap();
            mpsparser::parse(mps_content.as_slice()).unwrap()
        };
        let original_mps = &original_mps_obj;
    
        let original_path = path;
        // let output_file = &opt.output_file; //original_path.to_string_lossy().replace(".mps", ".sol");
                                            // if !output_file.to_lowercase().ends_with(".sol") {
                                            //     panic!("Unsupported output file extension.");
                                            // }
                                            // info!("checking domains and constraints... ");
                                            // lsmip::util::check_domains(original_mps);
                                            // lsmip::util::check_constraints(original_mps);
    
        let presolved: RwLock<Option<(PathBuf, MPSInstance)>> = RwLock::new(None);
        let presolved_borrowed = &presolved;
        let mut presolved_output_path = original_path.to_owned();
    
        let mut presolved_filename = original_path
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string();
        presolved_filename.push_str("_presolved.mps");
        presolved_output_path.set_file_name(presolved_filename);
    
        let mut presolved_solution_path = presolved_output_path.to_owned();
        presolved_solution_path.set_extension("sol");
    
        let (sol_tx, sol_rx) = std::sync::mpsc::channel::<(bool, Vec<f64>)>();
    
        crossbeam::scope(|s| {
            s.spawn(move |_| {
                let mut mps = presolved_borrowed.write().unwrap();
                let (presolved_path, presolved_mps) = presolve(
                    path,
                    &None,
                    "scripts/convert.py",
                    "tmp",
                );
                // info!("Presolved domains and constraints:");
                // lsmip::util::check_domains(&presolved_mps);
                // lsmip::util::check_constraints(&presolved_mps);
                *mps = Some((presolved_path, presolved_mps));
            });
    
            for thread_idx in [1, 2, 3, 4, 5, 6, 7] {
                // Thread settings.
    
                let seed = thread_idx as u8;
                let relax_continuous = thread_idx % 2 == 0;
                let use_presolved = (thread_idx / 4) % 2 == 1;
                let decay_factor = if (thread_idx / 2) % 2 == 0 || !USE_EXPONENTIAL_WEIGHT_DECAY {
                    1.0
                } else {
                    0.9999
                };
    
                let sol_tx = sol_tx.clone();
                s.spawn(move |_| {
                    let lock = use_presolved.then(|| presolved_borrowed.read().unwrap());
    
                    let mps = if use_presolved {
                        &lock.as_ref().unwrap().as_ref().unwrap().1
                    } else {
                        original_mps
                    };
    
                    let objective_coeffs = mps
                        .objective()
                        .iter()
                        .flat_map(|(_, v)| v.iter())
                        .map(|c| (c.var, c.coeff.as_f64() as f32))
                        .collect::<HashMap<usize, f32>>();
    
                    let mut solver =
                        lsmip::solver::Solver::with_seed(thread_idx as usize, seed, decay_factor);
                    for (var_idx, var) in mps.variables.iter().enumerate() {
                        let lb = var.lb.map(Number::as_f64).unwrap_or(0.) as f32;
                        let ub = var.ub.map(Number::as_f64).unwrap_or(UNBOUNDED) as f32; // TODO unbounded vars?
                        solver.add_var(
                            var.var_type,
                            lb,
                            ub,
                            objective_coeffs.get(&var_idx).copied().unwrap_or(0.),
                        );
                    }
    
                    for constraint in mps.constraints.iter() {
                        if matches!(constraint.rowtype, mpsparser::RowType::None) {
                            continue;
                        }
                        solver.add_constraint(
                            constraint.rowtype,
                            constraint.rhs.map(Number::as_f64).unwrap_or(0.) as f32,
                            &constraint
                                .cells
                                .iter()
                                .map(|c| (c.var, c.coeff.as_f64() as f32))
                                .collect::<Vec<_>>(),
                            relax_continuous,
                        );
                    }
    
                    solver.solve(|solution| {
                        let solution = solution.iter().map(|v| *v as f64).collect::<Vec<_>>();
                        sol_tx.send((use_presolved, solution)).unwrap();
                    });
                });
            }
    
            drop(sol_tx);
    
            let mut presolved_lock = None;
    
            let mut best_objective = f64::INFINITY;
            for (use_presolved, solution) in sol_rx.into_iter() {
                if use_presolved && presolved_lock.is_none() {
                    presolved_lock = Some(presolved.read().unwrap());
                }
    
                let mps = if use_presolved {
                    &presolved_lock.as_ref().unwrap().as_ref().unwrap().1
                } else {
                    original_mps
                };
    
                let path = if use_presolved {
                    &presolved_lock.as_ref().unwrap().as_ref().unwrap().0
                } else {
                    path
                };
    
                lsmip::lp::repair_continuous_and_save(
                    path,
                    original_path,
                    &mut |obj,model| {
                        let vars = model.get_vars().unwrap().iter().map(|v| {
                            let name = model.get_obj_attr(grb::attr::VarName, v).unwrap();
                            let value = model.get_obj_attr(grb::attr::X, v).unwrap();
                            (name, value as f32)
                        }).collect::<Vec<_>>();
                        output(obj, &mut vars.iter().map(|(n,v)| (n.as_str(),*v)));
                    },
                    &presolved_output_path,
                    &presolved_solution_path,
                    &mps.variables,
                    &solution,
                    &mut best_objective,
                );
    
                // if mpsparser::check_values(
                //     mps,
                //     &solution,
                //     mpsparser::DEFAULT_INT_TOLERANCE,
                //     mpsparser::DEFAULT_EQ_TOLERANCE,
                // )
                // .is_err()
                // {
                // }
    
                // if mpsparser::check_values(
                //     mps,
                //     &solution,
                //     mpsparser::DEFAULT_INT_TOLERANCE,
                //     mpsparser::DEFAULT_EQ_TOLERANCE,
                // )
                // .is_err()
                // {
                //     warn!("Invalid solution received.");
                //     continue;
                // }
    
                // let objective = match mps.objective() {
                //     Some((constant, obj)) => {
                //         let constant_term = constant.map(|n| n.as_f64()).unwrap_or(0.);
                //         -constant_term
                //             + obj
                //                 .iter()
                //                 .map(|c| c.coeff.as_f64() * solution[c.var])
                //                 .sum::<f64>()
                //     }
                //     None => 0.,
                // } as f32;
    
                // if objective < best_objective {
                //     best_objective = objective;
                //     output_solution(
                //         objective,
                //         &mut solution
                //             .iter()
                //             .enumerate()
                //             .map(|(idx, x)| (mps.variables[idx].name.as_str(), *x as f32)),
                //     );
                // }
            }
        })
        .unwrap();
    
    });

    hprof::profiler().print_timing();
}
