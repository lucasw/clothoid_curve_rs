// TODO(lucasw) maybe redundant with lib.rs
#![no_std]

use clothoid_bezier::f64::ClothoidBezierApproximation;
use clothoid_curve::f64::{Clothoid, Position};

use uom::num_traits::Zero;
use uom::si::{f64::Length, length::meter};

fn main() {
    let mut cd = Clothoid::default();

    cd.xy0.x = Length::new::<meter>(1.0);

    let cd = ClothoidBezierApproximation {
        curvature: 1.0,
        curvature_rate: 0.1,
        length: 10.0,
    };
    let clothoid = cd.to_clothoid();

    let pt = Position { x: Length::new::<meter>(2.0), y: Length::new::<meter>(1.0) };
    let (_pos, _distance, _iterations) = clothoid.get_nearest(&pt, Length::zero());
}
