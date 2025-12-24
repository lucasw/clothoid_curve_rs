pub use clothoid_bezier::f64::ClothoidBezierApproximation;

fn main() {
    let mut cd = test_no_std::Clothoid::default();

    cd.x0 = 1.0 as test_no_std::Float;

    let mut cd = ClothoidBezierApproximation {
        curvature: 1.0,
        curvature_rate: 0.1,
        length: 10.0,
    };
    let clothoid = cd.to_clothoid();
}
