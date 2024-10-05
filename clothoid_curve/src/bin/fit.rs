// Copyright 2018-2024 argmin developers
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

#[allow(unused_imports)]
use argmin::{
    core::{observers::ObserverMode, CostFunction, Error, Executor, Gradient},
    solver::{
        gradientdescent::SteepestDescent,
        linesearch::{HagerZhangLineSearch, MoreThuenteLineSearch},
    },
};
use argmin_observer_slog::SlogLogger;
use clothoid_curve::Clothoid;
use finitediff::FiniteDiff;
use std::f64::consts::{FRAC_PI_2, PI};

#[derive(Clone)]
struct Curve {
    x0: f64,
    y0: f64,
    theta0: f64,
    curvature0: f64,
    target_theta: f64,
}

impl CostFunction for Curve {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        // need to limited this to a maximum
        let curvature_rate = p[0];
        // TODO(lucasw) length needs to be > 0
        let length = p[1];
        let curve0 = Clothoid::create(
            self.x0,
            self.y0,
            self.theta0,
            self.curvature0,
            curvature_rate,
            length,
        );
        let curve_s = curve0.get_clothoid(length);
        let d_yaw = self.target_theta - curve_s.get_start_theta();
        let d_yaw = (d_yaw + FRAC_PI_2) % PI - FRAC_PI_2;
        let mut error = d_yaw * d_yaw + 0.1 * curvature_rate * curvature_rate;
        if length < 0.0 {
            error += length * length;
        }
        Ok(error)
    }
}

impl Gradient for Curve {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    fn gradient(&self, p: &Self::Param) -> Result<Self::Gradient, Error> {
        Ok(p.central_diff(&|x| self.cost(x).unwrap()))
    }
}

fn run() -> Result<(), Error> {
    // Define cost function (must implement `CostFunction` and `Gradient`)
    // Test making a left 90 degree turn
    let cost = Curve {
        x0: 0.0,
        y0: 0.0,
        theta0: 0.0,
        curvature0: 0.0,
        target_theta: -FRAC_PI_2,
    };

    // Define initial parameter vector
    // easy case
    // TODO(lucasw) If the initial curvature rate is opposite target direction the solution
    // will loop around the longer way- would have to test both to see which is lower cost
    // should make it so doing a loop is almost always the wrong thing
    let init_param: Vec<f64> = vec![0.01, 1.0];
    // tough case
    // let init_param: Vec<f64> = vec![-1.2, 1.0];

    // Pick a line search.
    // let linesearch = HagerZhangLineSearch::new();
    let linesearch = MoreThuenteLineSearch::new();

    // Set up solver
    let solver = SteepestDescent::new(linesearch);

    // Run solver
    let res = Executor::new(cost.clone(), solver)
        .configure(|state| state.param(init_param).max_iters(100))
        .add_observer(SlogLogger::term(), ObserverMode::Always)
        .run()?;

    // print result
    println!("{res}");

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
    let curve_end = curve_solution.get_clothoid(length);

    println!("{curve_solution:?}");
    println!("{curve_end:?}");
    Ok(())
}

fn main() {
    if let Err(ref e) = run() {
        println!("{e}");
        std::process::exit(1);
    }
}
