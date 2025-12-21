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
use clothoid_curve::f64::{Clothoid, angle_unwrap};
#[allow(unused_imports)]
use argmin_observer_slog::SlogLogger;
use finitediff::FiniteDiff;

#[derive(Clone)]
struct Curve {
    x0: f64,
    y0: f64,
    theta0: f64,
    curvature0: f64,
    target_theta: f64,
    target_curvature: f64,
}

impl Curve {
    fn from_clothoid(clothoid0: &Clothoid, target_theta: f64, target_curvature: f64) -> Self {
        Self {
            x0: clothoid0.x0,
            y0: clothoid0.y0,
            theta0: clothoid0.theta0,
            curvature0: clothoid0.curvature(),
            target_theta,
            target_curvature,
        }
    }

    fn to_clothoid(&self, curvature_rate: f64, length: f64) -> Clothoid {
        Clothoid::create(
            self.x0,
            self.y0,
            self.theta0,
            self.curvature0,
            curvature_rate,
            length,
        )
    }

    pub fn cost0(&self, curvature_rate: f64, length: f64) -> Result<f64, Error> {
        // TODO(lucasw) need to limit curvature_rate to a maximum
        let curve0 = self.to_clothoid(curvature_rate, length);
        let curve_s = curve0.get_clothoid(length);
        // println!("curvature rates {} -> {}",  curve0.curvature_rate(), curve_s.curvature_rate());
        // println!("curvature {} -> {}",  curve0.curvature(), curve_s.curvature());
        let d_yaw = self.target_theta - curve_s.get_start_theta();
        // println!("d_yaw: {d_yaw}");
        let d_yaw = angle_unwrap(d_yaw);
        // println!("unwrapped d_yaw: {d_yaw}");
        let d_curvature = self.target_curvature - curve_s.curvature();
        // println!("d_curvature: {d_curvature} = {} - {}, init curvature: {}", self.target_curvature, curve_s.curvature(), curve0.curvature());
        let mut residual = 0.5 * (d_yaw * d_yaw)
            + 4.0 * (d_curvature * d_curvature)
            + 0.1 * (curvature_rate * curvature_rate);
        // TODO(lucasw) length needs to be > 0
        if length < 0.0 {
            residual += length * length;
        }
        Ok(residual)
    }
}

impl CostFunction for Curve {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        self.cost0(p[0], p[1])
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
    target_theta: f64,
    target_curvature: f64,
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
    let init_param: Vec<f64> = vec![curvature_rate_guess, length_guess];
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
    let curvature_rate = param[0];
    // TODO(lucasw) length needs to be > 0
    let length = param[1];
    let curve_solution = Clothoid::create(
        cost.x0,
        cost.y0,
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
            target_theta,
            target_theta.to_degrees(),
            target_curvature
        );
    }
    Ok(curve_solution)
}
