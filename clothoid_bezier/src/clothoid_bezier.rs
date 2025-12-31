// Lucas Walter
// December 2025

#[cfg(feature = "argmin_fit")]
#[allow(unused_imports)]
use argmin::{
    core::{CostFunction, Error, Executor, Gradient, observers::ObserverMode},
    solver::{
        gradientdescent::SteepestDescent,
        linesearch::{HagerZhangLineSearch, MoreThuenteLineSearch},
    },
};
#[cfg(feature = "argmin_fit")]
#[allow(unused_imports)]
use argmin_observer_slog::SlogLogger;
#[cfg(feature = "argmin_fit")]
use finitediff::FiniteDiff;

use uom::num_traits::Zero;
use uom::si::{
    curvature::radian_per_meter,
    length::meter,
};

/// approximate a clothoid with a cubic bezier
/// not using uom types here so this can be used directly in egui
/// TODO(lucasw) make a uom-egui adapter
pub struct ClothoidBezierApproximation {
    pub curvature: NativeFloat,
    pub curvature_rate: NativeFloat,
    pub length: NativeFloat,
}

#[cfg(feature = "argmin_fit")]
extern crate alloc;
#[cfg(feature = "argmin_fit")]
use alloc::vec::Vec;
#[cfg(feature = "argmin_fit")]
use alloc::vec;

#[cfg(feature = "argmin_fit")]
impl CostFunction for &ClothoidBezierApproximation {
    type Param = Vec<NativeFloat>;
    type Output = NativeFloat;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        let handle_length0 = Length::new::<meter>(p[0]);
        let handle_length1 = Length::new::<meter>(p[1]);
        self.cost0(handle_length0, handle_length1)
    }
}

#[cfg(feature = "argmin_fit")]
impl Gradient for &ClothoidBezierApproximation {
    type Param = Vec<NativeFloat>;
    type Gradient = Vec<NativeFloat>;

    fn gradient(&self, p: &Self::Param) -> Result<Self::Gradient, Error> {
        // Ok(p.central_diff(&|x| self.cost(x).unwrap()))
        Ok(p.forward_diff(&|x| self.cost(x).unwrap()))
    }
}

impl ClothoidBezierApproximation {
    pub fn to_clothoid(&self) -> Clothoid {
        let x0 = 0.0;
        let y0 = 0.0;
        Clothoid::create(
            Length::new::<meter>(x0),
            Length::new::<meter>(y0),
            Angle::zero(),
            Curvature::new::<radian_per_meter>(self.curvature),
            curvature_per_meter(self.curvature_rate),
            Length::new::<meter>(self.length),
        )
    }

    pub fn get_bezier(
        clothoid_start: &Clothoid,
        handle0: Length,
        handle1: Length,
    ) -> CubicBezier2 {
        let clothoid_end = clothoid_start.get_end_clothoid();
        // the bezier end points
        let start = &clothoid_start.xy0;
        let end = &clothoid_end.xy0;

        // find where the handles are
        let clothoid_start_straight = clothoid_start.zero_curvature();
        let ctrl1 = clothoid_start_straight.get_xy(handle0);

        let clothoid_end_straight = clothoid_end.zero_curvature();
        let ctrl2 = clothoid_end_straight.get_xy(-handle1);

        CubicBezier2::new(start, &ctrl1, &ctrl2, end)
    }

    #[cfg(feature = "argmin_fit")]
    fn cost0(&self, handle0: Length, handle1: Length) -> Result<NativeFloat, Error> {
        let clothoid = self.to_clothoid();
        let bezier = Self::get_bezier(&clothoid, handle0, handle1);
        let curvature0 = clothoid.curvature();
        let delta0 = (curvature0 - bezier.curvature(ParametricTFrac::start())).get::<radian_per_meter>();

        let clothoid_end = clothoid.get_end_clothoid();
        let curvature1 = clothoid_end.curvature();
        let delta1 = (curvature1 - bezier.curvature(ParametricTFrac::end())).get::<radian_per_meter>();

        let length_delta = self.length - bezier.arclen(32).get::<meter>();

        let mut residual = delta0 * delta0 + delta1 * delta1 + 0.1 * length_delta * length_delta;

        if handle0 < Length::zero() {
            residual += -handle0.get::<meter>();
        }
        if handle1 < Length::zero() {
            residual += -handle1.get::<meter>();
        }

        /*
        let max = self.length * 4.0;

        let m0 = handle_length0 - max;
        if m0 > 0.0 {
            residual += m0;
        }
        let m1 = handle_length1 - max;
        if m1 > 0.0 {
            residual += m1;
        }
        info!("{handle_length0} {handle_length1} {residual}");
        */
        Ok(residual)
    }

    // TODO(lucasw) maybe argmin can be no_std
    #[cfg(feature = "argmin_fit")]
    pub fn find_handles(
        &self,
        handle0_guess: Length,
        handle1_guess: Length,
    ) -> Result<(NativeFloat, NativeFloat), Error> {
        let verbose = true;

        let init_param: Vec<NativeFloat> = vec![handle0_guess.get::<meter>(), handle1_guess.get::<meter>()];
        if verbose {
            /*
            info!(
                "initial cost {:?} -> {:?}",
                init_param,
                self.cost(&init_param)
            );
            */
        }

        // Pick a line search.
        // let linesearch = HagerZhangLineSearch::new();
        let linesearch = MoreThuenteLineSearch::new();
        let solver = SteepestDescent::new(linesearch);

        let res = Executor::new(self, solver)
            .configure(|state| state.param(init_param).target_cost(0.001).max_iters(150))
            // .add_observer(SlogLogger::term(), ObserverMode::Always)
            .run()?;

        // print result
        if verbose {
            /*
            info!("{res}");
            */
        }

        let param = res.state.param.unwrap();
        Ok((param[0], param[1]))
    }
}

#[cfg(test)]
mod clothoid_bezier_tests {
    use super::*;

    #[cfg(feature = "argmin_fit")]
    #[test]
    fn test_circle_fit() {
        for radius in [0.5, 1.0, 1.2] {
            let cba = ClothoidBezierApproximation {
                curvature: 1.0 / radius,
                curvature_rate: 0.0,
                length: radius * PI / 2.0,
            };

            // TODO(lucasw) this fails if too high or too low
            let scs = [0.45, 0.5, 0.6, 0.7, 0.8, 1.0];
            for sc0 in &scs {
                for sc1 in &scs {
                    let handle0_guess = Length::new::<meter>(radius * sc0);
                    let handle1_guess = Length::new::<meter>(radius * sc1);
                    let (handle0, handle1) = cba.find_handles(handle0_guess, handle1_guess).unwrap();
                    let expected = radius * 0.55;
                    let threshold = 0.02 * radius;
                    assert!((handle0 - expected).abs() < threshold, "sc {sc0}, r {radius}, guess {handle0_guess:?} -> {handle0:?} {expected:?}");
                    assert!((handle1 - expected).abs() < threshold, "sc {sc1}, r {radius}, guess {handle1_guess:?} -> {handle1:?} {expected:?}");
                }
            }
        }
    }
}
