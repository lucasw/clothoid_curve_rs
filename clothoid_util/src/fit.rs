// Copyright 2018-2024 argmin developers
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

#[allow(unused_imports)]
use argmin::{
    core::{CostFunction, Error, Executor, Gradient, observers::ObserverMode},
    solver::{
        gradientdescent::SteepestDescent,
        linesearch::{HagerZhangLineSearch, MoreThuenteLineSearch},
    },
};
#[allow(unused_imports)]
use argmin_observer_slog::SlogLogger;
use clothoid_curve::f64::{
    Clothoid, CurvaturePerLength, Position, angle_unwrap, curvature_per_meter,
    curvature_per_meter_float,
};
use finitediff::FiniteDiff;

use uom::num_traits::Zero;
use uom::si::{
    angle::{degree, radian},
    area::square_meter,
    curvature::radian_per_meter,
    f64::{Angle, Curvature, Length},
    length::meter,
};

#[derive(Clone)]
struct Curve {
    xy0: Position,
    theta0: Angle,
    curvature0: Curvature,
    target_theta: Angle,
    target_curvature: Curvature,
}

impl Curve {
    fn from_clothoid(
        clothoid0: &Clothoid,
        target_theta: Angle,
        target_curvature: Curvature,
    ) -> Self {
        Self {
            xy0: clothoid0.xy0.clone(),
            theta0: clothoid0.theta0,
            curvature0: clothoid0.curvature(),
            target_theta,
            target_curvature,
        }
    }

    fn to_clothoid(&self, curvature_rate: CurvaturePerLength, length: Length) -> Clothoid {
        Clothoid::create(
            self.xy0.x,
            self.xy0.y,
            self.theta0,
            self.curvature0,
            curvature_rate,
            length,
        )
    }

    pub fn cost0(&self, curvature_rate: CurvaturePerLength, length: Length) -> Result<f64, Error> {
        // TODO(lucasw) need to limit curvature_rate to a maximum
        let curve0 = self.to_clothoid(curvature_rate, length);
        let curve_s = curve0.get_clothoid(length);
        // println!("curvature rates {} -> {}",  curve0.curvature_rate(), curve_s.curvature_rate());
        // println!("curvature {} -> {}",  curve0.curvature(), curve_s.curvature());
        let d_yaw = self.target_theta - curve_s.get_start_theta();
        // println!("d_yaw: {d_yaw}");
        let d_yaw = angle_unwrap(d_yaw).get::<radian>();
        // println!("unwrapped d_yaw: {d_yaw}");
        let d_curvature = (self.target_curvature - curve_s.curvature()).get::<radian_per_meter>();
        // let curvature_rate = curvature_rate.get::<reciprocal_square_meter>();
        let curvature_rate = curvature_per_meter_float(curvature_rate);
        // println!("d_curvature: {d_curvature} = {} - {}, init curvature: {}", self.target_curvature, curve_s.curvature(), curve0.curvature());
        let mut residual = 0.5 * (d_yaw * d_yaw)
            + 4.0 * (d_curvature * d_curvature)
            + 0.1 * (curvature_rate * curvature_rate);
        // TODO(lucasw) length needs to be > 0
        if length < Length::zero() {
            residual += (length * length).get::<square_meter>();
        }
        Ok(residual)
    }
}

impl CostFunction for Curve {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        let curvature_rate = curvature_per_meter(p[0]);
        // CurvaturePerLength::new::<reciprocal_square_meter>(p[0]),
        self.cost0(curvature_rate, Length::new::<meter>(p[1]))
    }
}

impl Gradient for Curve {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    fn gradient(&self, p: &Self::Param) -> Result<Self::Gradient, Error> {
        // Ok(p.central_diff(&|x| self.cost(x).unwrap()))
        Ok(p.forward_diff(&|x| self.cost(x).unwrap()))
    }
}

pub fn find_clothoid(
    clothoid0: Clothoid, // initial conditions, curvature_rate and length are initial guess
    target_theta: Angle,
    target_curvature: Curvature,
) -> Result<Clothoid, Error> {
    // Define cost function (must implement `CostFunction` and `Gradient`)
    let cost = Curve::from_clothoid(&clothoid0, target_theta, target_curvature);

    let curvature_rate_guess = clothoid0.curvature_rate();
    let length_guess = clothoid0.length;

    let verbose = false;

    if verbose {
        println!(
            "initial cost {:?}",
            cost.cost0(curvature_rate_guess, length_guess)
        );
    }

    // if true {
    //    return Ok(clothoid0.clone());
    // }
    // Define initial parameter vector
    // easy case
    // TODO(lucasw) If the initial curvature rate is opposite target direction the solution
    // will loop around the longer way- would have to test both to see which is lower cost
    // should make it so doing a loop is almost always the wrong thing
    let init_param: Vec<f64> = vec![
        curvature_per_meter_float(curvature_rate_guess),
        length_guess.get::<meter>(),
    ];
    if verbose {
        println!("initial cost {:?}", cost.cost(&init_param));
    }
    // tough case
    // let init_param: Vec<f64> = vec![-1.2, 1.0];

    // Pick a line search.
    // let linesearch = HagerZhangLineSearch::new();
    let linesearch = MoreThuenteLineSearch::new();

    // Set up solver
    let solver = SteepestDescent::new(linesearch);

    // Run solver
    let res = Executor::new(cost.clone(), solver)
        .configure(|state| state.param(init_param).target_cost(0.01).max_iters(80))
        // .add_observer(SlogLogger::term(), ObserverMode::Always)
        .run()?;

    // print result
    if verbose {
        println!("{res}");
    }

    let param = res.state.param.unwrap();
    let curvature_rate = curvature_per_meter(param[0]);
    // TODO(lucasw) length needs to be > 0
    let length = Length::new::<meter>(param[1]);
    let curve_solution = Clothoid::create(
        cost.xy0.x,
        cost.xy0.y,
        cost.theta0,
        cost.curvature0,
        curvature_rate,
        length,
    );
    // let curve_end = curve_solution.get_clothoid(length);

    if verbose {
        // println!("solution start: {curve_solution:?}");
        // println!("solution end: {curve_end:?}");
        println!(
            "target theta {:0.3} ({:0.3}°), target_curvature {:0.3}",
            target_theta.get::<radian>(),
            target_theta.get::<degree>(),
            target_curvature.get::<radian_per_meter>(),
        );
    }
    Ok(curve_solution)
}
