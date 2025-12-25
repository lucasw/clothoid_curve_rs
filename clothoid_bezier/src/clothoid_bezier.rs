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

/// approximate a clothoid with a cubic bezier
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
        let handle_length0 = p[0];
        let handle_length1 = p[1];
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
            x0,
            y0,
            0.0,
            self.curvature,
            self.curvature_rate,
            self.length,
        )
    }

    pub fn get_bezier(
        clothoid: &Clothoid,
        handle_length0: NativeFloat,
        handle_length1: NativeFloat,
    ) -> CubicBezier2 {
        let clothoid_end = clothoid.get_end_clothoid();
        // the bezier end points
        let bz_pt0: [NativeFloat; 2] = [clothoid.x0, clothoid.y0];
        let bz_pt3: [NativeFloat; 2] = [clothoid_end.x0, clothoid_end.y0];

        // find where the handles are
        let bz_angle0 = clothoid.theta0;
        let bz_pt1 = [
            bz_pt0[0] + handle_length0 * cos(bz_angle0),
            bz_pt0[1] + handle_length0 * sin(bz_angle0),
        ];

        let bz_angle1 = clothoid_end.theta0;
        let bz_pt2 = [
            bz_pt3[0] - handle_length1 * cos(bz_angle1),
            bz_pt3[1] - handle_length1 * sin(bz_angle1),
        ];

        CubicBezier2::new(
            PointN::new(bz_pt0),
            PointN::new(bz_pt1),
            PointN::new(bz_pt2),
            PointN::new(bz_pt3),
        )
    }

    #[cfg(feature = "argmin_fit")]
    fn cost0(&self, handle_length0: NativeFloat, handle_length1: NativeFloat) -> Result<NativeFloat, Error> {
        let clothoid = self.to_clothoid();
        let bezier = Self::get_bezier(&clothoid, handle_length0, handle_length1);
        let curvature0 = clothoid.curvature();
        let delta0 = curvature0 - bezier.curvature(0.0);

        let clothoid_end = clothoid.get_end_clothoid();
        let curvature1 = clothoid_end.curvature();
        let delta1 = curvature1 - bezier.curvature(1.0);

        let length_delta = self.length - bezier.arclen(32);

        let mut residual = delta0 * delta0 + delta1 * delta1 + 0.1 * length_delta * length_delta;

        if handle_length0 < 0.0 {
            residual += -handle_length0;
        }
        if handle_length1 < 0.0 {
            residual += -handle_length1;
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
        handle0_guess: NativeFloat,
        handle1_guess: NativeFloat,
    ) -> Result<(NativeFloat, NativeFloat), Error> {
        let verbose = true;

        let init_param: Vec<NativeFloat> = vec![handle0_guess, handle1_guess];
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
