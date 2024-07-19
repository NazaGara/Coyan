use coyan_fta::{self, solver::get_solver_from_path};
const EPSILON: f64 = f64::EPSILON; // 2.2204460492503131E-16f64
#[test]
fn xor() {
    let solver_cmd = String::from("../solvers/gpmc");
    let filename = "../tests/xor.dft";
    let ft = coyan_fta::fault_tree::FaultTree::new_from_file(&filename, true);
    let solver = get_solver_from_path(&solver_cmd);
    let tep = solver.compute_probabilty(&ft, coyan_fta::formula::CNFFormat::MC21, 1.0, 100);
    let true_tep = 0.6925;
    assert!(f64::abs(true_tep - tep) < EPSILON);
    // assert_eq!(true_tep, tep);
}

#[test]
fn not() {
    let solver_cmd = String::from("../solvers/gpmc");
    let filename = "../tests/not.dft";
    let ft = coyan_fta::fault_tree::FaultTree::new_from_file(filename, true);
    let solver = get_solver_from_path(&solver_cmd);
    let tep = solver.compute_probabilty(&ft, coyan_fta::formula::CNFFormat::MC21, 1.0, 100);
    let true_tep = 1.0 - 0.25;
    assert!(f64::abs(true_tep - tep) < EPSILON);
    // assert_eq!(true_tep, tep);
}
#[test]
fn and() {
    let solver_cmd = String::from("../solvers/gpmc");
    let filename = "../tests/and.dft";
    let ft = coyan_fta::fault_tree::FaultTree::new_from_file(filename, true);
    let solver = get_solver_from_path(&solver_cmd);
    let tep = solver.compute_probabilty(&ft, coyan_fta::formula::CNFFormat::MC21, 1.0, 100);
    let true_tep = 0.25 * (0.35 * 0.45);
    assert!(f64::abs(true_tep - tep) < EPSILON);
    // assert_eq!(true_tep, tep);
}

#[test]
fn or() {
    let solver_cmd = String::from("../solvers/gpmc");
    let filename = "../tests/or.dft";
    let ft = coyan_fta::fault_tree::FaultTree::new_from_file(filename, true);
    let solver = get_solver_from_path(&solver_cmd);
    let tep = solver.compute_probabilty(&ft, coyan_fta::formula::CNFFormat::MC21, 1.0, 100);
    let true_tep = 1.0 - ((1.0 - 0.25) * (1.0 - 0.35) * (1.0 - 0.45));
    assert!(f64::abs(true_tep - tep) < EPSILON);
    // assert_eq!(true_tep, tep);
}

#[test]
fn vot() {
    let solver_cmd = String::from("../solvers/gpmc");
    let filename = "../tests/3of5.dft";
    let ft = coyan_fta::fault_tree::FaultTree::new_from_file(filename, true);
    let solver = get_solver_from_path(&solver_cmd);
    let tep = solver.compute_probabilty(&ft, coyan_fta::formula::CNFFormat::MC21, 1.0, 100);
    let true_tep = 0.403040625; //Obtained from Storm-DFT
    assert!(f64::abs(true_tep - tep) < EPSILON);
    // assert_eq!(true_tep, tep);
}

#[test]
fn ffort_sample() {
    let solver_cmd = String::from("../solvers/gpmc");
    let solver = get_solver_from_path(&solver_cmd);
    let filename0 = "../tests/ogpf.dft";
    let filename1 = "../tests/pt.dft";
    let filename2 = "../tests/rbc.dft";
    let ft0 = coyan_fta::fault_tree::FaultTree::new_from_file(filename0, true);
    let ft1 = coyan_fta::fault_tree::FaultTree::new_from_file(filename1, true);
    let ft2 = coyan_fta::fault_tree::FaultTree::new_from_file(filename2, true);
    let tep = (
        solver.compute_probabilty(&ft0, coyan_fta::formula::CNFFormat::MC21, 1.0, 100),
        solver.compute_probabilty(&ft1, coyan_fta::formula::CNFFormat::MC21, 1.0, 100),
        solver.compute_probabilty(&ft2, coyan_fta::formula::CNFFormat::MC21, 1.0, 100),
    );
    let true_tep = (0.9769512022, 3.501334916e-05, 2.231940945e-10);
    // True TEP obtained from Storm-DFT, epsilon lower because Coyan handles more precision digits.
    assert!(f64::abs(true_tep.0 - tep.0) < 1e-5);
    assert!(f64::abs(true_tep.1 - tep.1) < 1e-5);
    assert!(f64::abs(true_tep.2 - tep.2) < 1e-5);
}
