/// Copyright 2024 Lucas Walter
///
///
use clothoid_curve::clothoid::{angle_unwrap, Clothoid};
use clothoid_curve::fit::find_clothoid;

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let theta0: f64 = angle_unwrap(args[1].parse().unwrap());
    let curvature0: f64 = args[2].parse().unwrap();
    let target_theta: f64 = angle_unwrap(args[3].parse().unwrap());
    let target_curvature: f64 = args[4].parse().unwrap();

    let curvature_rate_guess = target_curvature - curvature0;
    let length_guess = 1.0;

    // TODO(lucasw) 'create' isn't that helpful vs. building the object directly, the latter
    // is better because it shows which parameter is which
    let clothoid0 = Clothoid::create(
        0.0,
        0.0,
        theta0,
        curvature0,
        curvature_rate_guess,
        length_guess,
    );
    let _rv = find_clothoid(clothoid0, target_theta, target_curvature);
}
