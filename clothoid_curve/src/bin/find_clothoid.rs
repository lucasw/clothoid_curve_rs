use clothoid_curve::fit::find_clothoid;

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let theta0: f64 = args[1].parse().unwrap();
    let curvature0: f64 = args[2].parse().unwrap();
    let target_theta: f64 = args[3].parse().unwrap();
    let target_curvature: f64 = args[3].parse().unwrap();
    let _rv = find_clothoid(0.0, 0.0, theta0, curvature0, target_theta, target_curvature);
}
