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

    let length_guess = {
        if target_curvature != 0.0 {
            (target_theta - theta0) / target_curvature
        } else {
            1.0
        }
    };
    let curvature_rate_guess = (target_curvature - curvature0) / length_guess;
    println!("length guess: {length_guess:0.3}, curvature_rate: {curvature_rate_guess:0.3}");

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
