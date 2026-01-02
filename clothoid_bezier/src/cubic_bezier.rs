/// Adapted from stroke cubic_bezier
// use serde::{Deserialize, Serialize};
use super::point::Point;
use super::LineSegment;
use super::QuadraticBezier;
use super::*;
use uom::si::{
    curvature::radian_per_meter,
    length::meter,
};

fn position_to_point2(pos: &Position) -> PointN<2> {
    let xy: [NativeFloat; 2] = [pos.x.get::<meter>(), pos.y.get::<meter>()];
    PointN::new(xy)
}

fn point2_to_position(p: PointN<2>) -> Position {
    Position {
        x: Length::new::<meter>(p.axis(0)),
        y: Length::new::<meter>(p.axis(1)),
    }
}

// TODO(lucasw) can I avoid the redundant sizes?
type CubicBezier2Base = CubicBezier::<PointN<2>, 2>;
#[derive(Clone)]
pub struct CubicBezier2(CubicBezier2Base);

#[derive(Clone, Copy)]
pub struct ParametricTFrac(NativeFloat);

impl ParametricTFrac {
    pub fn start() -> Self {
        Self(0.0)
    }

    pub fn end() -> Self {
        Self(1.0)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
// #[repr(transparent)]
pub struct EuclideanTFrac(pub NativeFloat);

impl EuclideanTFrac {
    pub fn start() -> Self {
        Self(0.0)
    }

    pub fn end() -> Self {
        Self(1.0)
    }
}

impl CubicBezier2 {
    pub fn new(start: &Position, ctrl1: &Position, ctrl2: &Position, end: &Position) -> Self {
        // TODO(lucasw) optionally cache length?
        Self(CubicBezier2Base::new(
            position_to_point2(start),
            position_to_point2(ctrl1),
            position_to_point2(ctrl2),
            position_to_point2(end),
        ))
    }

    pub fn start(&self) -> Position {
        Position::from_array_meter([self.0.start.axis(0), self.0.start.axis(1)])
    }

    pub fn ctrl1(&self) -> Position {
        Position::from_array_meter([self.0.ctrl1.axis(0), self.0.ctrl1.axis(1)])
    }

    pub fn ctrl2(&self) -> Position {
        Position::from_array_meter([self.0.ctrl2.axis(0), self.0.ctrl2.axis(1)])
    }

    pub fn end(&self) -> Position {
        Position::from_array_meter([self.0.end.axis(0), self.0.end.axis(1)])
    }

    pub fn arclen(&self, nsteps: usize) -> Length {
        Length::new::<meter>(self.0.arclen(nsteps))
    }

    // TODO(lucasw) provide another fn with tolerance
    pub fn arclen_castlejau(&self) -> Length {
        Length::new::<meter>(self.0.arclen_castlejau(None))
    }

    /// require the length be computed elsehwere, though may want to optionally have a cached
    /// length to use after computing it the first time
    pub fn euclidean_to_parametric_t(&self, euclidean_t: EuclideanTFrac, bezier_length: Length) -> (Length, ParametricTFrac) {
        let desired_length = euclidean_t.0 * bezier_length;
        let (achieved_len, parametric_tfrac) = self.0.desired_len_to_parametric_t(desired_length.get::<meter>(), None);
        // TODO(lucasw) also return achieved_len / bezier_length as achieved_euclidean_t?
        (Length::new::<meter>(achieved_len), ParametricTFrac(parametric_tfrac))
    }

    pub fn eval(&self, t: ParametricTFrac) -> Position {
        point2_to_position(self.0.eval(t.0))
    }

    pub fn curvature(&self, t: ParametricTFrac) -> Curvature {
        Curvature::new::<radian_per_meter>(self.0.curvature(t.0))
    }
}

/// A 2D  cubic Bezier curve defined by four points: the starting point, two successive control points and the ending point.
///
/// The curve is defined by equation
/// : ```∀ t ∈ [0..1],  P(t) = (1 - t)³ * start + 3 * (1 - t)² * t * ctrl1 + 3 * t² * (1 - t) * ctrl2 + t³ * end```
// TODO(lucasw) the vim syntax highlighting doesn't like the triple slashes followed by triple
// back-tick above, but putting the ':' on the same line avoids it
#[derive(Copy, Clone, Debug, PartialEq)] // , Deserialize, Serialize)]
pub struct CubicBezier<P, const PDIM: usize> {
    pub start: P,
    pub ctrl1: P,
    pub ctrl2: P,
    pub end: P,
}

//#[allow(dead_code)]
impl<P, const PDIM: usize> CubicBezier<P, PDIM>
where
    P: Point,
{
    pub fn new(start: P, ctrl1: P, ctrl2: P, end: P) -> Self {
        CubicBezier {
            start,
            ctrl1,
            ctrl2,
            end,
        }
    }

    /// Evaluate a CubicBezier curve at t by direct evaluation of the polynomial (not numerically stable)
    pub fn eval(&self, t: NativeFloat) -> P {
        self.start * ((-t + 1.0) * (-t + 1.0) * (-t + 1.0))
            + self.ctrl1 * (t * (-t + 1.0) * (-t + 1.0) * 3.0)
            + self.ctrl2 * (t * t * (-t + 1.0) * 3.0)
            + self.end * (t * t * t)
    }

    /// Evaluate a CubicBezier curve at t using the numerically stable De Casteljau algorithm
    pub fn eval_casteljau(&self, t: NativeFloat) -> P {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl1 - self.start) * t;
        let ctrl_1bc = self.ctrl1 + (self.ctrl2 - self.ctrl1) * t;
        let ctrl_1cd = self.ctrl2 + (self.end - self.ctrl2) * t;
        // second iteration
        let ctrl_2ab = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;
        let ctrl_2bc = ctrl_1bc + (ctrl_1cd - ctrl_1bc) * t;
        // third iteration, return final point on the curve ctrl_3ab
        ctrl_2ab + (ctrl_2bc - ctrl_2ab) * t
    }

    pub fn control_points(&self) -> [P; 4] {
        [self.start, self.ctrl1, self.ctrl2, self.end]
    }

    /// Returns the x coordinate of the curve evaluated at t
    /// Convenience shortcut for bezier.eval(t).x()
    pub fn axis(&self, t: NativeFloat, axis: usize) -> NativeFloat {
        let t2 = t * t;
        let t3 = t2 * t;
        let one_t = -t + 1.0;
        let one_t2 = one_t * one_t;
        let one_t3 = one_t2 * one_t;

        one_t3 * self.start.axis(axis)
            + one_t2 * t * self.ctrl1.axis(axis) * 3.0
            + one_t * t2 * self.ctrl2.axis(axis) * 3.0
            + t3 * self.end.axis(axis)
    }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// Remember arclen also works by linear approximation, not the integral, so we have to accept error!
    /// This approximation is unfeasable if desired accuracy is greater than 2 decimal places
    pub fn arclen(&self, nsteps: usize) -> NativeFloat {
        let stepsize = 1.0 / (nsteps as NativeFloat);
        let mut arclen = NativeFloat::from(0.0);
        for t in 1..nsteps {
            let t = (t as NativeFloat) * 1.0 / (nsteps as NativeFloat);
            let p1 = self.eval_casteljau(t);
            let p2 = self.eval_casteljau(t + stepsize);

            arclen += sqrt((p1 - p2).squared_length());
        }
        arclen
    }

    /// from graphite bezier
    /// either compute the full length of the curve, or find the parametric t that results
    /// in the desired length
    /// return length, parametric_t
    fn recurse(a0: P, a1: P, a2: P, a3: P,
            desired_len: Option<NativeFloat>, tolerance: NativeFloat, level: u8, min_level: u8) -> (NativeFloat, NativeFloat) {
        // TODO(lucasw) when the line is straight, these values are the same so this exits
        // immediately, but in that case the handles need to be 1/3 the line length to work
        // properly
        // But it looks like forcing this to recurse 4 times avoids this issue, but only
        // do that when computing parametric t, don't use it for computing arclen
        let lower = a0.distance(&a3);
        let upper = a0.distance(&a1) + a1.distance(&a2) + a2.distance(&a3);
        let in_tolerance = (upper - lower) <= (2. * tolerance);
        let over_level_threshold = level >= 8;
        if (in_tolerance || over_level_threshold) && level >= min_level {
            let approx_len = (lower + upper) / 2.;
            let parametric_t = match desired_len {
                Some(desired_len) => {
                    desired_len / approx_len
                }
                None => {
                    1.0
                }
            };
            return (approx_len, parametric_t)
        }

        let b1 = (a0 + a1) * 0.5;
        let t0 = (a1 + a2) * 0.5;
        let c1 = (a2 + a3) * 0.5;
        let b2 = (b1 + t0) * 0.5;
        let c2 = (t0 + c1) * 0.5;
        let b3 = (b2 + c2) * 0.5;
        let (first_len, t) = Self::recurse(a0, b1, b2, b3, desired_len, 0.5 * tolerance, level + 1, min_level);
        let new_desired_len = {
            if let Some(desired_len) = desired_len {
                if first_len > desired_len {
                    return (first_len, t * 0.5)
                }

                Some(desired_len - first_len)
            } else {
                None
            }
        };
        let (second_len, t) = Self::recurse(b3, c2, c1, a3, new_desired_len, 0.5 * tolerance, level + 1, min_level);
        (first_len + second_len, t * 0.5 + 0.5)
    }

    /// Use Casteljau subdivision, noting that the length is more than the straight line distance from start to end but less than the straight line distance through the handles
    pub fn arclen_castlejau(&self, tolerance: Option<NativeFloat>) -> NativeFloat {
        let (approx_len, _) = Self::recurse(self.start, self.ctrl1, self.ctrl2, self.end, None, tolerance.unwrap_or_default(), 0, 0);
        approx_len
    }

    // return the achieved length (which may be slightly off desired) and the parametric t value
    // that resulted in it
    pub fn desired_len_to_parametric_t(&self, desired_len: NativeFloat, tolerance: Option<NativeFloat>) -> (NativeFloat, NativeFloat) {
        let start_level = 0;
        let min_level = 4;
        let (len, parametric_t) = Self::recurse(self.start, self.ctrl1, self.ctrl2, self.end, Some(desired_len), tolerance.unwrap_or_default(), start_level, min_level);
        (len, parametric_t)
    }

    pub fn split(&self, t: NativeFloat) -> (Self, Self) {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl1 - self.start) * t;
        let ctrl_1bc = self.ctrl1 + (self.ctrl2 - self.ctrl1) * t;
        let ctrl_1cd = self.ctrl2 + (self.end - self.ctrl2) * t;
        // second iteration
        let ctrl_2ab = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;
        let ctrl_2bc = ctrl_1bc + (ctrl_1cd - ctrl_1bc) * t;
        // third iteration, final point on the curve
        let ctrl_3ab = ctrl_2ab + (ctrl_2bc - ctrl_2ab) * t;

        (
            CubicBezier {
                start: self.start,
                ctrl1: ctrl_1ab,
                ctrl2: ctrl_2ab,
                end: ctrl_3ab,
            },
            CubicBezier {
                start: ctrl_3ab,
                ctrl1: ctrl_2bc,
                ctrl2: ctrl_1cd,
                end: self.end,
            },
        )
    }

    /// Return the derivative curve.
    /// The derivative is also a bezier curve but of degree n-1 (cubic->quadratic)
    /// Since it returns the derivative function, eval() needs to be called separately
    pub fn derivative(&self) -> QuadraticBezier<P, PDIM> {
        QuadraticBezier {
            start: (self.ctrl1 - self.start) * 3.0,
            ctrl: (self.ctrl2 - self.ctrl1) * 3.0,
            end: (self.end - self.ctrl2) * 3.0,
        }
    }

    /// Direct Derivative - Sample the axis coordinate at 'axis' of the curve's derivative at t
    /// without creating a new curve. This is a convenience function for .derivative().eval(t).axis(n)  
    /// Parameters:
    ///   t: the sampling parameter on the curve interval [0..1]
    ///   axis: the index of the coordinate axis [0..N]
    /// Returns:
    ///   Scalar value of the points own type type F  
    /// May be deprecated in the future.  
    /// This function can cause out of bounds panic when axis is larger than dimension of P
    pub fn dd(&self, t: NativeFloat, axis: usize) -> NativeFloat {
        let t2 = t * t;
        let c0 = t * -3.0 + t * 6.0 - 3.0;
        let c1 = t2 * 9.0 - t * 12.0 + 3.0;
        let c2 = t2 * -9.0 + t * 6.0;
        let c3 = t2 * 3.0;

        self.start.axis(axis) * c0
            + self.ctrl1.axis(axis) * c1
            + self.ctrl2.axis(axis) * c2
            + self.end.axis(axis) * c3
    }

    /// Return the tangent at position t
    pub fn tangent(&self, t: NativeFloat) -> P {
        let derivative = self.derivative().eval(t);
        let d_len = sqrt(derivative.squared_length());
        let scale = 1.0 / d_len;
        derivative * scale
    }

    // curvature is 1.0 / radius, or 0.0 when radius is close to zero
    pub fn curvature(&self, t: NativeFloat) -> NativeFloat {
        let d = self.derivative();
        let d_t = d.eval(t);
        let dx = d_t.axis(0);
        let dy = d_t.axis(1);

        let dd = d.derivative();
        let dd_t = dd.eval(t);
        let ddx = dd_t.axis(0);
        let ddy = dd_t.axis(1);

        let numerator = dx * ddy - dy * ddx;
        let denominator = pow(dx * dx + dy * dy, 1.5);
        // this is what graphite bezier solver does
        if denominator.abs() < 1e-3 {
            0.0
        } else {
            numerator / denominator
        }
    }

    // pub fn radius(&self, t: NativeFloat) -> F
    // where
    // F: NativeFloatloat,
    // P:  Sub<P, Output = P>
    //     + Add<P, Output = P>
    //     + Mul<F, Output = P>,
    // NativeFloat: Sub<F, Output = F>
    //     + Add<F, Output = F>
    //     + Mul<F, Output = F>
    //     + Float
    //     + Into
    // {
    //     return 1.0 / self.curvature(t)
    // }

    fn capture_closest(
        closest_point_on_curve: &mut P,
        tmin: &mut NativeFloat,
        dmin: &mut NativeFloat,
        candidate_on_curve: &P,
        t: NativeFloat,
        point: &P,
    ) {
        let distance = (*candidate_on_curve - *point).squared_length();
        if distance <= *dmin {
            *tmin = t;
            *closest_point_on_curve = *candidate_on_curve;
            *dmin = distance;
        }
    }

    /// Calculates the minimum distance between given 'point' and the curve, returns the closest point on
    /// the curve, the t-value, and the distance.
    /// Uses two passes with the same amount of steps in t:
    /// 1. coarse search over the whole curve
    /// 2. fine search around the minimum yielded by the coarse search
    pub fn closest_to_point(&self, point: P) -> (P, NativeFloat, NativeFloat) {
        // this gets to a good enough closeness to the true value,
        // though maybe some curves produce worse results, and of course
        // the caller may want better results or to not do so many
        // iterations
        let nsteps = 32;
        self.closest_to_point_with_nsteps(point, nsteps)
    }

    // TODO(lucasw) just use argmin to search the line?
    // TODO(lucasw) want a version of this that takes a starting position and only searches
    // forward a little- maybe make the fine pass a separate function and just call that
    /// brute force search along curve to find closest point- the less curvature there is
    /// the fewer steps needed- TODO(lucasw) ideally could figure that out, or analyze the curve
    /// ahead of time.
    /// Could keep a cache in the curve so eval(t) can return immediately if this is called
    /// repeatedly for the same curve (or fix nsteps ahead of time so it's certain new t values
    /// won't be generated here.
    /// Also probably could use a much coarser first pass and then use the two lowest sequential
    /// distances to determine a straight line segment, find the closest point on that segment then
    /// evaluate the real eval(t) of that t value and then refine a few more times around the found
    /// value.
    /// return the closest point on the line, the t-value (make that s-value), and the distance
    pub fn closest_to_point_with_nsteps(&self, point: P, nsteps: usize) -> (P, NativeFloat, NativeFloat) {
        let mut tmin: NativeFloat = 0.5;
        let mut dmin_squared: NativeFloat = 1e6; // 2.0 * (point - self.start).squared_length();
        let mut closest_point_on_curve = self.start;

        // 1. coarse pass
        // the coarse pass and fine pass combined
        {
            let step = 1.0 / nsteps as NativeFloat;
            for i in 0..nsteps {
                // calculate next step value
                let t = i as NativeFloat * step;
                // calculate distance to candidate
                let candidate = self.eval(t);
                Self::capture_closest(
                    &mut closest_point_on_curve,
                    &mut tmin,
                    &mut dmin_squared,
                    &candidate,
                    t,
                    &point,
                );
            }
        }
        self.closest_to_point_with_nsteps_fine_pass(point, nsteps, closest_point_on_curve, tmin, dmin_squared)
    }

    /// brute force search over subset of curve with set starting conditions, search forward and
    /// backward of the initial t-value
    pub fn closest_to_point_with_nsteps_fine_pass(&self, point: P, nsteps: usize, closest_point0: P, tmin0: NativeFloat, dmin_squared0: NativeFloat) -> (P, NativeFloat, NativeFloat) {
        let mut tmin = tmin0;
        let mut dmin_squared = dmin_squared0;
        let mut closest_point_on_curve = closest_point0;

        // TODO(lucasw) pass in the range and compute the fine_step from it
        let fine_step = 1.0 / (nsteps * nsteps) as NativeFloat;
        let half_range = nsteps as NativeFloat * fine_step / 2.0;
        let tmin_coarse = tmin;
        for i in 0..nsteps {
            // calculate next step value ( a 64th of a 64th from first step)
            let t = tmin_coarse + i as NativeFloat * fine_step - half_range;
            // calculate distance to candidate centered around tmin from before
            let candidate: P = self.eval(t);
            Self::capture_closest(
                &mut closest_point_on_curve,
                &mut tmin,
                &mut dmin_squared,
                &candidate,
                t,
                &point,
            );
        }
        (closest_point_on_curve, tmin, sqrt(dmin_squared))
    }

    pub fn distance_to_point(&self, point: P) -> NativeFloat {
        let (_, _, distance) = self.closest_to_point(point);
        distance
    }

    pub fn distance_to_point_with_nsteps(&self, point: P, nsteps: usize) -> NativeFloat {
        let (_, _, distance) = self.closest_to_point_with_nsteps(point, nsteps);
        distance
    }

    pub fn baseline(&self) -> LineSegment<P, PDIM> {
        LineSegment {
            start: self.start,
            end: self.end,
        }
    }

    pub fn is_linear(&self, tolerance: NativeFloat) -> bool {
        // if start and end are (nearly) the same
        if (self.start - self.end).squared_length() < EPSILON {
            return false;
        }
        // else check if ctrl points lie on baseline
        self.are_points_colinear(tolerance)
    }

    fn are_points_colinear(&self, tolerance: NativeFloat) -> bool {
        let line = self.baseline();
        line.distance_to_point(self.ctrl1) <= tolerance
            && line.distance_to_point(self.ctrl2) <= tolerance
    }

    // Returs if the whole set of control points can be considered one singular point
    // given some tolerance.
    // TODO use machine epsilon vs squared_length OK?
    pub fn is_a_point(&self, tolerance: NativeFloat) -> bool {
        let tolerance_squared = tolerance * tolerance;
        // Use <= so that tolerance can be zero.
        (self.start - self.end).squared_length() <= tolerance_squared
            && (self.start - self.ctrl1).squared_length() <= tolerance_squared
            && (self.end - self.ctrl2).squared_length() <= tolerance_squared
    }

    /// Compute the real roots of the cubic bezier function with
    /// parameters of the form a*t^3 + b*t^2 + c*t + d for each dimension
    /// using cardano's algorithm (code adapted from github.com/nical/lyon)
    /// returns an ArrayVec of the present roots (max 3)
    #[allow(clippy::many_single_char_names)] // this is math, get over it
    pub(crate) fn real_roots(
        &self,
        a: NativeFloat,
        b: NativeFloat,
        c: NativeFloat,
        d: NativeFloat,
    ) -> ArrayVec<[NativeFloat; 3]> {
        let mut result = ArrayVec::new();

        // check if can be handled below cubic order
        if a.abs() < EPSILON {
            if b.abs() < EPSILON {
                if c.abs() < EPSILON {
                    // no solutions
                    return result;
                }
                // is linear equation
                result.push(-d / c);
                return result;
            }
            // is quadratic equation
            let delta = c * c - b * d * 4.0;
            if delta > 0.0 {
                let sqrt_delta = sqrt(delta);
                result.push((-c - sqrt_delta) / (b * 2.0));
                result.push((-c + sqrt_delta) / (b * 2.0));
            } else if delta.abs() < EPSILON {
                result.push(-c / (b * 2.0));
            }
            return result;
        }

        // is cubic equation -> use cardano's algorithm
        let frac_1_3 = NativeFloat::from(1.0 / 3.0);

        let bn = b / a;
        let cn = c / a;
        let dn = d / a;

        let delta0: NativeFloat = (cn * 3.0 - bn * bn) / 9.0;
        let delta1: NativeFloat = (bn * cn * 9.0 - dn * 27.0 - bn * bn * bn * 2.0) / 54.0;
        let delta_01: NativeFloat = delta0 * delta0 * delta0 + delta1 * delta1;

        if delta_01 >= NativeFloat::from(0.0) {
            let delta_p_sqrt: NativeFloat = delta1 + sqrt(delta_01);
            let delta_m_sqrt: NativeFloat = delta1 - sqrt(delta_01);

            let s = pow(delta_p_sqrt.signum() * delta_p_sqrt.abs(), frac_1_3);
            let t = pow(delta_m_sqrt.signum() * delta_m_sqrt.abs(), frac_1_3);

            result.push(-bn * frac_1_3 + (s + t));

            // Don't add the repeated root when s + t == 0.
            if (s - t).abs() < EPSILON && (s + t).abs() >= EPSILON {
                result.push(-bn * frac_1_3 - (s + t) / 2.0);
            }
        } else {
            // TODO(lucasw) exact same code is in root.rs?
            let theta = acos(delta1 / sqrt(-delta0 * delta0 * delta0));
            let two_sqrt_delta0 = sqrt(-delta0) * 2.0;
            result.push(two_sqrt_delta0 * cos(theta * frac_1_3) - bn * frac_1_3);
            result
                .push(two_sqrt_delta0 * cos((theta + 2.0 * PI) * frac_1_3) - bn * frac_1_3);
            result
                .push(two_sqrt_delta0 * cos((theta + 4.0 * PI) * frac_1_3) - bn * frac_1_3);
        }

        result
    }

    /// Solves the cubic bezier function given a particular coordinate axis value
    /// by solving the roots for the axis functions
    /// Parameters:
    /// value: the coordinate value on the particular axis
    /// axis: the index of the axis
    /// Returns those roots of the function that are in the interval [0.0, 1.0].
    #[allow(dead_code)]
    fn solve_t_for_axis(&self, value: NativeFloat, axis: usize) -> ArrayVec<[NativeFloat; 3]> {
        let mut result = ArrayVec::new();
        // check if all points are the same or if the curve is really just a line
        if self.is_a_point(EPSILON)
            || (self.are_points_colinear(EPSILON)
                && (self.start - self.end).squared_length() < EPSILON)
        {
            return result;
        }
        let a = -self.start.axis(axis) + self.ctrl1.axis(axis) * 3.0 - self.ctrl2.axis(axis) * 3.0
            + self.end.axis(axis);
        let b =
            self.start.axis(axis) * 3.0 - self.ctrl1.axis(axis) * 6.0 + self.ctrl2.axis(axis) * 3.0;
        let c = -self.start.axis(axis) * 3.0 + self.ctrl1.axis(axis) * 3.0;
        let d = self.start.axis(axis) - value;

        let roots = self.real_roots(a, b, c, d);
        for root in roots {
            if root > 0.0 && root < 1.0 {
                result.push(root);
            }
        }

        result
    }

    /// Return the bounding box of the curve as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box(&self) -> [(NativeFloat, NativeFloat); PDIM] {
        // calculate coefficients for the derivative: at^2 + bt + c
        // from the expansion of the cubic bezier curve: sum_i=0_to_3( binomial(3, i) * t^i * (1-t)^(n-i) )
        // yields coeffcients
        // po: [1, -2,  1]
        // p1: [0,  2, -2]
        // p2: [0,  0,  1]
        //      c   b   a
        let mut bounds = [(0.0, 0.0); PDIM];
        let derivative = self.derivative();
        // calculate coefficients for derivative
        let a: P = derivative.start + derivative.ctrl * -2.0 + derivative.end;
        let b: P = derivative.start * -2.0 + derivative.ctrl * 2.0;
        let c: P = derivative.start;

        // calculate roots for t over x axis and plug them into the bezier function
        //  to get x,y values (make vec 2 bigger for t=0,t=1 values)
        // loop over any of the points dimensions (they're all the same)
        for (dim, _) in a.into_iter().enumerate() {
            let mut extrema: ArrayVec<[NativeFloat; 4]> = ArrayVec::new();
            extrema.extend(
                derivative
                    .real_roots(a.axis(dim), b.axis(dim), c.axis(dim))
                    .into_iter(),
            );
            // only retain roots for which t is in [0..1]
            extrema.retain(|root| -> bool { root > &mut 0.0 && root < &mut 1.0 });
            // evaluates roots in original function
            for t in extrema.iter_mut() {
                *t = self.eval_casteljau(*t).axis(dim);
            }
            // add y-values for start and end point as candidates
            extrema.push(self.start.axis(dim));
            extrema.push(self.end.axis(dim));
            // sort to get min and max values for bounding box
            extrema.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            // determine xmin, xmax, ymin, ymax, from the set {B(xroots), B(yroots), B(0), B(1)}
            // (Intermediate control points can't form a boundary)
            // .unwrap() is ok as it can never be empty as it always at least contains the endpoints
            bounds[dim] = (extrema[0], *extrema.last().unwrap());
        }
        bounds
    }
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::PointN;
    use super::*;
    use std::println;
    use core::f64::consts::PI;
    #[test]
    fn circle_approximation_error() {
        // define closure for unit circle
        let circle_error =
            |p: Position, radius: Length| -> Length { (p.x * p.x + p.y * p.y).sqrt() - radius };

        let radius = 1.0;
        // define control points for 4 bezier segments
        // control points are chosen for minimum radial distance error
        // according to: http://spencermortensen.com/articles/bezier-circle/
        // TODO don't hardcode values
        let c = 0.551915024494 * radius;
        let max_drift_perc = 0.019608; // radial drift percent
        let max_error = max_drift_perc * 0.011 * radius; // absolute max radial error

        let bezier_quadrants = [
            CubicBezier2::new(
                &Position::from_array_meter([0.0, radius]),
                &Position::from_array_meter([c, radius]),
                &Position::from_array_meter([radius, c]),
                &Position::from_array_meter([radius, 0.0]),
            ),
            CubicBezier2::new(
                &Position::from_array_meter([radius, 0.0]),
                &Position::from_array_meter([radius, -c]),
                &Position::from_array_meter([c, -radius]),
                &Position::from_array_meter([0.0, -radius]),
            ),
            CubicBezier2::new(
                &Position::from_array_meter([0.0, -radius]),
                &Position::from_array_meter([-c, -radius]),
                &Position::from_array_meter([-radius, -c]),
                &Position::from_array_meter([-radius, 0.0]),
            ),
            CubicBezier2::new(
                &Position::from_array_meter([-radius, 0.0]),
                &Position::from_array_meter([-radius, c]),
                &Position::from_array_meter([-c, radius]),
                &Position::from_array_meter([0.0, radius]),
            ),
        ];
        let radius = Length::new::<meter>(radius);
        let nsteps = 1000;
        for t in 0..=nsteps {
            let t = ParametricTFrac(t as NativeFloat / (nsteps as NativeFloat));

            for bezier_quadrant in &bezier_quadrants {
                let point = bezier_quadrant.eval(t);
                let error = circle_error(point, radius).get::<meter>();
                assert!(error.abs() <= max_error, "{} <= {max_error}", error.abs());
            }
        }
    }

    #[test]
    fn circle_circumference_approximation() {
        // define control points for 4 cubic bezier segments to best approximate a unit circle
        // control points are chosen for minimum radial distance error, see circle_approximation_error() in this file
        // given this, the circumference will also be close to 2*pi
        // (remember arclen also works by linear approximation, not the true integral, so we have to accept error)!
        // This approximation is unfeasable if desired accuracy is greater than ~2 decimal places (at 1000 steps)
        // TODO don't hardcode values, solve for them
        let c = 0.551915024494 as NativeFloat;
        let max_error = 1e-2 as NativeFloat;
        let nsteps = 1e3 as usize;
        let tau = 2. * PI as NativeFloat;

        let bezier_quadrant_1 = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 1.0]),
            ctrl1: PointN::new([c, 1.0]),
            ctrl2: PointN::new([1.0, c]),
            end: PointN::new([1.0, 0.0]),
        };
        let bezier_quadrant_2 = CubicBezier::<_, 2> {
            start: PointN::new([1.0, 0.0]),
            ctrl1: PointN::new([1.0, -c]),
            ctrl2: PointN::new([c, -1.0]),
            end: PointN::new([0.0, -1.0]),
        };
        let bezier_quadrant_3 = CubicBezier::<_, 2> {
            start: PointN::new([0.0, -1.0]),
            ctrl1: PointN::new([-c, -1.0]),
            ctrl2: PointN::new([-1.0, -c]),
            end: PointN::new([-1.0, 0.0]),
        };
        let bezier_quadrant_4 = CubicBezier::<_, 2> {
            start: PointN::new([-1.0, 0.0]),
            ctrl1: PointN::new([-1.0, c]),
            ctrl2: PointN::new([-c, 1.0]),
            end: PointN::new([0.0, 1.0]),
        };
        let circumference = (bezier_quadrant_1.arclen(nsteps)
            + bezier_quadrant_2.arclen(nsteps)
            + bezier_quadrant_3.arclen(nsteps)
            + bezier_quadrant_4.arclen(nsteps)) as NativeFloat;
        //dbg!(circumference);
        //dbg!(tau);
        assert!(((tau + max_error) > circumference) && ((tau - max_error) < circumference),
            "{tau} +/- {max_error} >/< {circumference}");

        let circumference = (bezier_quadrant_1.arclen_castlejau(None)
            + bezier_quadrant_2.arclen_castlejau(None)
            + bezier_quadrant_3.arclen_castlejau(None)
            + bezier_quadrant_4.arclen_castlejau(None)) as NativeFloat;
        let max_error = max_error * 0.1;
        assert!(((tau + max_error) > circumference) && ((tau - max_error) < circumference),
            "{tau} +/- {max_error} >/< {circumference}");
    }

    #[test]
    fn find_paramteric_t() {
        let c = 0.551915024494;
        let bezier0 = CubicBezier::<_, 2> {
            start: PointN::new([-1.0, 0.0]),
            ctrl1: PointN::new([-1.0, c]),
            ctrl2: PointN::new([-c, 1.0]),
            end: PointN::new([0.0, 1.0]),
        };

        // TODO(lucasw) for straight lines the handles need to be 1/3 the length for the parametric
        // t function to work
        let bezier_straight_thirds = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 0.0]),
            ctrl1: PointN::new([3.0, 0.0]),
            ctrl2: PointN::new([6.0, 0.0]),
            end: PointN::new([9.0, 9.0]),
        };

        let bezier_diag_neg = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 0.0]),
            ctrl1: PointN::new([-1.0, -1.0]),
            ctrl2: PointN::new([-9.0, -9.0]),
            end: PointN::new([-10.0, -10.0]),
        };

        let bezier_diag = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 0.0]),
            ctrl1: PointN::new([1.0, 1.0]),
            ctrl2: PointN::new([9.0, 9.0]),
            end: PointN::new([10.0, 10.0]),
        };

        let bezier_ud = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 0.0]),
            ctrl1: PointN::new([0.0, 1.0]),
            ctrl2: PointN::new([0.0, 9.0]),
            end: PointN::new([0.0, 10.0]),
        };

        let bezier_lr = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 0.0]),
            ctrl1: PointN::new([1.0, 0.0]),
            ctrl2: PointN::new([9.0, 0.0]),
            end: PointN::new([10.0, 0.0]),
        };

        for bezier in [bezier0, bezier_straight_thirds, bezier_diag_neg, bezier_diag, bezier_ud, bezier_lr] {
        // for bezier in [bezier0, bezier_straight_thirds, bezier_lr] {
            let b_len = bezier.arclen_castlejau(None);
            for sc in [0.0, 0.05, 0.1, 0.2, 0.5, 0.6, 0.9, 1.0] {
                let desired_len = sc * b_len;
                let (len, parametric_t) = bezier.desired_len_to_parametric_t(desired_len, None);
                let (left, _right) = bezier.split(parametric_t);
                let achieved_len = left.arclen_castlejau(None);
                assert!((desired_len - achieved_len).abs() < 0.1, "{bezier:?}\neuclidean t {sc} -> parametric t value {parametric_t}, desired_len {desired_len} -> {len} or split {achieved_len}");
            }

            for t in [0.0, 0.05, 0.1, 0.5, 0.8, 1.0] {
                let (left, right) = bezier.split(t);
                let left_len = left.arclen_castlejau(None);
                let (_len, t2) = bezier.desired_len_to_parametric_t(left_len, None);

                let right_len = right.arclen_castlejau(None);
                assert!((left_len + right_len - b_len).abs() < 0.001, "{left_len} + {right_len} = {} = {b_len}", left_len + right_len);
                assert!(
                    (t - t2).abs() < 0.01,
                    "parametric t {t} -> length {left_len} / {b_len} -> solved parametric t {t2}",
                );
            }
        }
    }

    #[test]
    fn eval_equivalence_casteljau() {
        // all eval methods should be approximately equivalent for well defined test cases
        // and not equivalent where numerical stability becomes an issue for normal eval
        let bezier = CubicBezier::<_, 2>::new(
            PointN::new([0.0, 1.77]),
            PointN::new([1.1, -1.0]),
            PointN::new([4.3, 3.0]),
            PointN::new([3.2, -4.0]),
        );

        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as NativeFloat * 1.0 / (nsteps as NativeFloat);
            let p1 = bezier.eval(t);
            let p2 = bezier.eval_casteljau(t);
            let err = p2 - p1;
            assert!(err.squared_length() < EPSILON);
        }
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 1.77]),
            ctrl1: PointN::new([2.9, 0.0]),
            ctrl2: PointN::new([4.3, 3.0]),
            end: PointN::new([3.2, -4.0]),
        };
        // split it at an arbitrary point
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // compare left and right subcurves with parent curve
        // take the difference of the two points which must not exceed the absolute error
        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as NativeFloat * 1.0 / (nsteps as NativeFloat);
            // left
            let mut err = bezier.eval(t / 2.0) - left.eval(t);
            assert!(err.squared_length() < EPSILON);
            // right
            err = bezier.eval((t * 0.5) + 0.5) - right.eval(t);
            assert!(err.squared_length() < EPSILON);
        }
    }

    #[test]
    fn bounding_box_contains() {
        // check if bounding box for a curve contains all points (with some approximation error)
        let bezier = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 1.77]),
            ctrl1: PointN::new([2.9, 0.0]),
            ctrl2: PointN::new([4.3, -3.0]),
            end: PointN::new([3.2, 4.0]),
        };

        let bounds = bezier.bounding_box();

        let max_err = 1e-2;

        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = t as NativeFloat * 1.0 / (nsteps as NativeFloat);
            let p = bezier.eval_casteljau(t);
            //dbg!(t);
            //dbg!(p);
            //dbg!(xmin-max_err, ymin-max_err, xmax+max_err, ymax+max_err);
            for (idx, axis) in p.into_iter().enumerate() {
                assert!((axis >= (bounds[idx].0 - max_err)) && (axis <= (bounds[idx].1 + max_err)))
            }
        }
    }

    #[test]
    fn closest_to_point_on_straight_line() {
        let bezier = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 0.0]),
            ctrl1: PointN::new([1.0, 0.0]),
            ctrl2: PointN::new([9.0, 0.0]),
            end: PointN::new([10.0, 0.0]),
        };

        for desired_len in [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0] {
            let (_found_len, parametric_t) = bezier.desired_len_to_parametric_t(desired_len, None);
            let expected_distance = 2.0;
            let point_off_line = PointN::new([desired_len, expected_distance]);
            let (closest_point_on_line, t, distance) = bezier.closest_to_point(point_off_line);
            println!("{desired_len} {parametric_t} -> {closest_point_on_line:?} {t} {distance}");
            assert!((t - parametric_t).abs() < 0.01, "t {t}, expected {parametric_t}");
            assert!((distance - expected_distance).abs() < 0.002, "distance {distance}, expected {expected_distance}");
        }
    }

    #[test]
    fn distance_to_point() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let curve = CubicBezier::<_, 2> {
            start: PointN::new([0.0, 1.77]),
            ctrl1: PointN::new([1.1, -1.0]),
            ctrl2: PointN::new([4.3, 3.0]),
            end: PointN::new([3.2, -4.0]),
        };
        assert!(
            curve.distance_to_point(PointN::new([-5.1, -5.6]))
                > curve.distance_to_point(PointN::new([5.1, 5.6]))
        );

        // make sure the default nsteps isn't that bad
        for pt_x in [-3.0, -2.0, 0.0, 4.0, 5.0, 6.0] {
            let pt = PointN::new([pt_x, 4.0]);
            let good_nsteps = 512;
            let good_distance = curve.distance_to_point_with_nsteps(pt, good_nsteps);
            let test_distance = curve.distance_to_point(pt);
            assert!((good_distance - test_distance).abs() < 0.15, "distance {test_distance:.4}, expected {good_distance:.4}");
        }
    }
}
