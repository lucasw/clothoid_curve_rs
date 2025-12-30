// this doesn't have to be no_std, but with no-default-features it should be no_std compatible
#![no_std]

use core::include;
pub mod f32 {
    use super::*;

    pub type Float = f32;
    use core::f32::consts::{FRAC_2_SQRT_PI, FRAC_PI_2, PI};
    use libm::cosf as cos;
    use libm::floorf as floor;
    use libm::sinf as sin;
    use libm::sqrtf as sqrt;

    use uom::si::f32::{Angle, Area, Length, V};

    include!("clothoid.rs");
}

pub mod f64 {
    use super::*;

    pub type Float = f64;
    use core::f64::consts::{FRAC_2_SQRT_PI, FRAC_PI_2, PI};
    use libm::cos;
    use libm::floor;
    use libm::sin;
    use libm::sqrt;

    use uom::si::f64::{Angle, Area, Length, V};

    /*
    use std::fmt;
    impl fmt::Debug for Clothoid {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let mut msg = "Clothoid: [".to_string()
                + &format!("\n\tx: {:0.3}, y: {:0.3}", self.x0, self.y0)
                + &format!(
                    "\n\ttheta0 (initial yaw/heading): {:0.3}radians ({:0.3}°)",
                    self.theta0,
                    self.theta0.to_degrees()
                )
                + &format!("\n\tkappa0 (curvature 1/r): {:0.3}", self.kappa0);
            if self.kappa0 != 0.0 {
                msg += &format!(" radius: {:.3}", 1.0 / self.kappa0);
            }
            let msg =
                msg + &format!(
                    "\n\tdk (curvature rate, curvature per unit length): {:0.6}",
                    self.dk
                ) + &format!("\n\tlength: {:0.3}", self.length)
                    + "\n]";
            write!(f, "{}", msg)
        }
    }
    */

    include!("clothoid.rs");
}

#[cfg(test)]
mod tests {
    use super::f64::Curvature;
    use super::*;
    use uom::si::f64::ReciprocalLength;
    use uom::si::reciprocal_length::reciprocal_meter;

    #[test]
    fn curvature_vs_reciprocal_length() {
        let c0 = Curvature::new::<reciprocal_meter>(0.5);
        let c1 = ReciprocalLength::new::<reciprocal_meter>(0.5);

        // TODO(lucasw) this passes, would need a newtype to prevent it
        assert_eq!(c0, c1);

        // are these basically type aliases?
        let _c2: Curvature = c1.clone();
        let _c3: ReciprocalLength = c0.clone();
    }
}
