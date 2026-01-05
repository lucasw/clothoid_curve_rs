/// Copyright 2024 Lucas Walter
///
///
use clothoid_curve::f64::{
    AngleCosSin, Clothoid, CurvaturePerLength, Float, Position, angle_unwrap,
};
use clothoid_util::fit::find_clothoid;

use uom::num_traits::Zero;
use uom::si::{
    angle::radian,
    curvature::radian_per_meter,
    // area::square_meter,
    f64::{Angle, Curvature, Length},
    length::meter,
};

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let theta0: Float = args[1].parse().unwrap();
    let curvature0: Float = args[2].parse().unwrap();
    let target_theta: Float = args[3].parse().unwrap();
    let target_curvature: Float = args[4].parse().unwrap();

    let theta0 = AngleCosSin::from_angle(angle_unwrap(Angle::new::<radian>(theta0)));
    let target_theta = AngleCosSin::from_angle(angle_unwrap(Angle::new::<radian>(target_theta)));
    let target_curvature = Curvature::new::<radian_per_meter>(target_curvature);

    let length_guess: Length = {
        if target_curvature != Curvature::zero() {
            // TODO(lucasw) The real quantity of clothoid length is radian meter
            ((target_theta.angle - theta0.angle) / target_curvature) / Angle::new::<radian>(1.0)
        } else {
            Length::new::<meter>(1.0)
        }
    };

    let curvature0 = Curvature::new::<radian_per_meter>(curvature0);
    let curvature_rate_guess: CurvaturePerLength =
        ((target_curvature - curvature0) / length_guess).into();
    println!(
        "length guess: {:?}, curvature_rate: {:?}",
        length_guess, curvature_rate_guess
    );

    let clothoid0 = Clothoid::create(
        Position::default(),
        theta0,
        curvature0,
        curvature_rate_guess,
        length_guess,
    );
    let _rv = find_clothoid(clothoid0, target_theta, target_curvature);
}
