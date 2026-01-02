/*!
Lucas Walter

December 2025

show how well a set of bezier curves can approximate a clothoid curve
*/
#![no_std]

use tinyvec::ArrayVec;
// export common types at crate root

pub mod f32 {
    use super::*;

    pub type NativeFloat = f32;
    const EPSILON: core::primitive::f32 = core::primitive::f32::EPSILON;

    // use clothoid_curve::f32::curvature_per_meter;
    use clothoid_curve::f32::Position;

    use core::f32::consts::PI;

    // TODO(lucasw) use std versions if available
    use libm::acosf as acos;
    use libm::atan2f as atan2;
    use libm::cosf as cos;
    use libm::powf as pow;
    // use libm::sinf as sin;
    use libm::sqrtf as sqrt;

    use uom::si::f32::{Angle, Curvature, Length};

    // specialized types
    pub mod cubic_bezier {
        include!("cubic_bezier.rs");
    }
    pub mod line {
        include!("line.rs");
    }
    pub mod quadratic_bezier {
        include!("quadratic_bezier.rs");
    }

    // Traits
    pub mod point {
        include!("point.rs");
    }

    mod roots {
        include!("roots.rs");
    }

    pub use cubic_bezier::{CubicBezier2, ParametricTFrac};
    pub use line::LineSegment;
    pub use point::Point;
    pub use point::PointN;
    pub use quadratic_bezier::QuadraticBezier;

    // TODO(lucasw) the argmin stuff doesn't work with f32, make it only use f64
    // use clothoid_curve::f32::Clothoid;
    // include!("clothoid_bezier.rs");
}

pub mod f64 {
    use super::*;

    pub type NativeFloat = f64;
    const EPSILON: core::primitive::f64 = core::primitive::f64::EPSILON;

    use clothoid_curve::f64::{Position, curvature_per_meter};

    use core::f64::consts::PI;
    use libm::acos;
    use libm::atan2;
    use libm::cos;
    use libm::pow;
    // use libm::sin;
    use libm::sqrt;

    use uom::si::f64::{Angle, Curvature, Length};

    // specialized types
    pub mod cubic_bezier {
        include!("cubic_bezier.rs");
    }
    pub mod line {
        include!("line.rs");
    }
    pub mod quadratic_bezier {
        include!("quadratic_bezier.rs");
    }

    // Traits
    pub mod point {
        include!("point.rs");
    }

    mod roots {
        include!("roots.rs");
    }

    pub use cubic_bezier::{CubicBezier2, ParametricTFrac};
    pub use line::LineSegment;
    pub use point::Point;
    pub use point::PointN;
    pub use quadratic_bezier::QuadraticBezier;

    use clothoid_curve::f64::Clothoid;
    // pub mod clothoid_bezier {
    include!("clothoid_bezier.rs");
    // }
}
