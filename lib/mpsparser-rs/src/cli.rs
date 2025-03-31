use std::io::Write;
use std::path::Path;

use log::info;

use crate::{parse, MPSInstance};

pub fn with_cli_mpsinstance(
    f: impl FnOnce(&MPSInstance, &mut dyn FnMut(f32, &mut dyn Iterator<Item = (&str, f32)>)),
) {
    with_cli_mpsinstance_path(|_p,mps,inner_f| f(mps, inner_f))
}

pub fn with_cli_mpsinstance_path(
    f: impl FnOnce(&Path, &MPSInstance, &mut dyn FnMut(f32, &mut dyn Iterator<Item = (&str, f32)>)),
) {
    let _p = hprof::enter("read and parse");
    let args = std::env::args().collect::<Vec<_>>();
    assert!(args.len() == 3);
    let filename = Path::new(&args[1]);
    let outdir = Path::new(&args[2]);
    assert!(filename.is_file());
    assert!(outdir.is_dir());

    println!("instance:{}", filename.file_name().unwrap().to_string_lossy());
    let mps = if filename.to_string_lossy().ends_with(".mps.gz") {
        let file = std::fs::File::open(filename).unwrap();
        let decoder = flate2::read::GzDecoder::new(file);
        parse(std::io::BufReader::new(decoder)).unwrap()
    } else if filename.to_string_lossy().ends_with(".mps") {
        let mps_content = std::fs::read(filename).unwrap();
        parse(mps_content.as_slice()).unwrap()
    } else {
        panic!("Unsupported file extension.");
    };

    let t0 = std::time::Instant::now();

    drop(_p);

    let mut output_solution = |obj: f32, values: &mut dyn Iterator<Item = (&str, f32)>| {
        let _p = hprof::enter("write solution file");
        let t = (std::time::Instant::now() - t0).as_secs_f32();
        let solution_filename = format!(
            "{:.2}_{}_{}.sol",
            t,
            obj,
            filename.file_name().unwrap().to_string_lossy()
        );
        let solution_path = outdir.join(solution_filename);
        info!("Outputting to {:?}", solution_path);
        let mut file = std::fs::File::create(&solution_path).expect("Unable to create file");
        for (name,val) in values {
            writeln!(file, "{} {}", name, val).unwrap();
        }
        info!("Saved solution to file {:?}", solution_path);
    };

    f(filename, &mps, &mut output_solution);
}
