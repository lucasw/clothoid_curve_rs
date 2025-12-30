pub use clothoid_bezier::f64::ClothoidBezierApproximation;

use uom::si::{f32::Length, length::meter};

fn main() {
    let mut cd = test_no_std::Clothoid::default();

    cd.xy0.x = Length::new::<meter>(1.0 as test_no_std::Float);

    let cd = ClothoidBezierApproximation {
        curvature: 1.0,
        curvature_rate: 0.1,
        length: 10.0,
    };
    let _clothoid = cd.to_clothoid();
}
