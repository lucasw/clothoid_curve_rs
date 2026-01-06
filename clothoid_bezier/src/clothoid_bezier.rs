// Lucas Walter
// December 2025

use clothoid_curve::f64::curvature_per_meter_float;

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
#[cfg(feature = "std")]
use tracing::info;
#[cfg(not(feature = "std"))]
macro_rules! info {
    ($($tt:tt)*) => {};
}

use uom::num_traits::Zero;
use uom::si::{
    // angle::degree,
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
            Position {
                x: Length::new::<meter>(x0),
                y: Length::new::<meter>(y0),
            },
            AngleCosSin::default(),
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

/// approximate a bezier with a clothoid
/// only using curvature rate and length with error defined only by the
/// distance from the end of the clothoid to the end of the bezier curve
/// TODO(lucasw) also the end curvature needs to be constrained?
/// Or that can be a two-clothoid process, split a bezier, don't constrain
/// curvature (or position?) in the middle just at the ends.
/// Would like to be able to join Clothoid segments together with
/// constrained positions and curvatures but that will require some
/// relaxation somewhere.
/// The initial position, angle is defined by by the bezier curve t = 0.0 start.
/// The target clothoid end position by the t=1.0 end
pub struct BezierToClothoid {
    /// store the bezier curve but only use the derived clothoids below during
    /// finding
    pub bezier: CubicBezier2,
    /// this is derived from the bezier, which
    pub start_clothoid: Clothoid,
    /// cached target end clothoid
    pub target_end_clothoid: Clothoid,
}

#[cfg(feature = "argmin_fit")]
impl CostFunction for &BezierToClothoid {
    type Param = Vec<NativeFloat>;
    type Output = NativeFloat;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        let curvature_per_length = curvature_per_meter(p[0]);
        let length = Length::new::<meter>(p[1]);
        self.cost0(curvature_per_length, length)
    }
}

#[cfg(feature = "argmin_fit")]
impl Gradient for &BezierToClothoid {
    type Param = Vec<NativeFloat>;
    type Gradient = Vec<NativeFloat>;

    fn gradient(&self, p: &Self::Param) -> Result<Self::Gradient, Error> {
        // Ok(p.central_diff(&|x| self.cost(x).unwrap()))
        Ok(p.forward_diff(&|x| self.cost(x).unwrap()))
    }
}

impl BezierToClothoid {
    pub fn from_bezier(bezier: CubicBezier2) -> Self {
        let curvature_per_length = CurvaturePerLength::zero();
        let length = bezier.arclen_castlejau();

        let start_clothoid = {
            let start_pos = bezier.start();
            // TODO(lucasw) bezier.derivative() is called redundantly in these calls
            let t = ParametricTFrac(0.0);
            let angle = bezier.angle(t);
            // TODO(lucasw) there may be situations where the initial curvature needs
            // to be overridden
            let curvature = bezier.curvature(t);

            Clothoid::create(
                start_pos,
                angle,
                curvature,
                curvature_per_length,
                length,
            )
        };

        let target_end_clothoid = {
            let end_pos = bezier.end();
            let t = ParametricTFrac(1.0);
            let angle = bezier.angle(t);
            // TODO(lucasw) likely will override this
            let curvature = bezier.curvature(t);

            Clothoid::create(
                end_pos,
                angle,
                curvature,
                curvature_per_length,
                length,
            )
        };
        info!("target end {target_end_clothoid:?}");

        Self {
            bezier,
            start_clothoid,
            target_end_clothoid,
        }
    }

    #[cfg(feature = "argmin_fit")]
    fn cost0(&self, curvature_per_length: CurvaturePerLength, length: Length) -> Result<NativeFloat, Error> {
        let mut clothoid = self.start_clothoid.clone();
        clothoid.set_curvature_rate(curvature_per_length);
        clothoid.length = length;
        let end_clothoid = clothoid.get_end_clothoid();
        let delta_pos = self.target_end_clothoid.xy0 - end_clothoid.xy0;
        let dx = delta_pos.x.get::<meter>();
        let dy = delta_pos.y.get::<meter>();
        // info!("{:?} - {:?} -> {dx} {dy}", self.target_end_clothoid.xy0, clothoid.xy0);
        let cost = dx * dx + dy * dy;

        // can't improve on angle given the constraints
        /*
        let delta_angle = (self.target_end_clothoid.theta0.angle - end_clothoid.theta0.angle)
            .get::<degree>().abs();
        let cost = cost + 0.01 * delta_angle;
        */

        Ok(cost)
    }

    #[cfg(feature = "argmin_fit")]
    pub fn find_clothoid(
        &self,
        curvature_per_length_guess: CurvaturePerLength,
        length_guess: Length,
    ) -> Result<Clothoid, Error> {
        let verbose = true;

        let init_param: Vec<NativeFloat> = vec![
            curvature_per_meter_float(curvature_per_length_guess),
            length_guess.get::<meter>(),
        ];
        if verbose {
            info!(
                "initial cost {:?} -> {:?}",
                init_param,
                self.cost(&init_param)
            );
        }

        // Pick a line search.
        // let linesearch = HagerZhangLineSearch::new();
        let linesearch = MoreThuenteLineSearch::new();
        let solver = SteepestDescent::new(linesearch);

        let res = Executor::new(self, solver)
            .configure(|state| state.param(init_param).target_cost(0.000001).max_iters(150))
            // .add_observer(SlogLogger::term(), ObserverMode::Always)
            .run()?;

        // print result
        if verbose {
            info!("{res}");
        }

        let param = res.state.param.unwrap();
        let curvature_per_length = curvature_per_meter(param[0]);
        let length = Length::new::<meter>(param[1]);
        Ok(Clothoid::create(
            self.start_clothoid.xy0,
            self.start_clothoid.theta0,
            self.start_clothoid.curvature(),
            curvature_per_length,
            length,
        ))
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
