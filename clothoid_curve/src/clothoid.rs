/*
#[cfg(not(feature = "f64"))]
pub type Float = f32;
#[cfg(not(feature = "f64"))]
use core::f32::consts::{FRAC_2_SQRT_PI, FRAC_PI_2, PI};

#[cfg(feature = "f64")]
pub type Float = f64;
#[cfg(feature = "f64")]
use std::f64::consts::{FRAC_2_SQRT_PI, FRAC_PI_2, PI};
*/

use core::clone::Clone;
use core::cmp::PartialEq;
use core::iter::Iterator;
use core::default::Default;
use core::prelude::rust_2024::derive;
use core::ops::Neg;

use serde::{Deserialize, Serialize};

use typenum::{N2, Z0};
use uom::si::marker::AngleKind;
use uom::num_traits::Zero;
use uom::si::{
    ISQ, Quantity, SI,
    angle::radian,
    area::square_meter,
    length::meter,
};

/// put angle into -pi, pi range
pub fn angle_unwrap(angle: Angle) -> Angle {
    Angle::new::<radian>((angle.get::<radian>() + PI) % (2.0 * PI) - PI)
}

#[cfg(feature = "std")]
extern crate alloc;
#[cfg(feature = "std")]
use alloc::vec::Vec;

/*
dimension: ISQ<
        N2,     // length
        Z0,     // mass
        Z0,     // time
        Z0,     // electric current
        Z0,     // thermodynamic temperature
        Z0,     // amount of substance
        Z0>;    // luminous intensity
*/
pub type CurvaturePerLength = Quantity<ISQ<N2, Z0, Z0, Z0, Z0, Z0, Z0, dyn AngleKind>, SI<V>, V>;

// TODO(lucasw) need a custom reciprocal_square_meter unit
// impl CurvaturePerLength {
pub fn curvature_per_meter(val: Float) -> CurvaturePerLength {
    (1.0 / Area::new::<square_meter>(1.0 / val)).into()
}

/// turn CurvaturePerLength into a float
pub fn curvature_per_meter_float(cpl: CurvaturePerLength) -> Float {
    let area: Area = (1.0 / cpl).into();
    1.0 / area.get::<square_meter>()
}
//}

// TODO(lucasw) put Position in a different file
#[derive(Default, Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Position {
    pub x: Length,
    pub y: Length,
}

impl Neg for Position {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl Position {
    pub fn from_array_meter(p: [Float; 2]) -> Self {
        Position {
            x: Length::new::<meter>(p[0]),
            y: Length::new::<meter>(p[1]),
        }
    }

    pub fn as_array_meter(&self) -> [Float; 2] {
        [self.x.get::<meter>(), self.y.get::<meter>()]
    }
}

// TODO(lucasw) 'the method `clamp` exists for struct `Quantity<..., ..., f64>`, but its trait
// bounds were not satisfied'
pub fn curvature_clamp(curvature: Curvature, min: Curvature, max: Curvature) -> Curvature {
    if curvature > max {
        max
    } else if curvature < min {
        min
    } else {
        curvature
    }
}

#[derive(Clone, PartialEq)]
pub struct Clothoid {
    /// start point xy
    pub xy0: Position,
    /// start point theta/yaw/heading
    pub theta0: Angle,
    /// cached cos and sin of theta
    pub cos_theta0: Float,
    pub sin_theta0: Float,
    kappa0: Curvature,     // start point curvature 1/r
    dk: CurvaturePerLength, // curvature rate, how much curvature changes per unit length, end theta will be
    // theta(s) = theta + theta' * s + 1/2 * theta0'' * s^2
    // theta(s) = theta + kappa0 * s + 1/2 * dk * s^2
    // with s = length
    pub length: Length, // how long the curve is (end kappa will be length * dk)
}

impl Default for Clothoid {
    fn default() -> Self {
        Self {
            xy0: Position::default(),
            theta0: Angle::zero(),
            cos_theta0: 1.0,
            sin_theta0: 0.0,
            kappa0: Curvature::zero(),
            dk: CurvaturePerLength::zero(),
            length: Length::zero(),
        }
    }
}

/*
// This function calculates the fresnel cosine and sine integrals.
// Input:
// y = values for which fresnel integrals have to be evaluated
//
// Output:
// FresnelC = fresnel cosine integral of y
// FresnelS = fresnel sine integral of y
//
// Adapted from:
// Atlas for computing mathematical functions : an illustrated guide for
// practitioners, with programs in C and Mathematica / William J. Thompson.
// New York : Wiley, c1997.
//
// Author: Venkata Sivakanth Telasula
// email: sivakanth.telasula@gmail.com
// date: August 11, 2005
*/
#[allow(clippy::excessive_precision)]
const FRN: &[Float] = &[
    0.49999988085884732562,
    1.3511177791210715095,
    1.3175407836168659241,
    1.1861149300293854992,
    0.7709627298888346769,
    0.4173874338787963957,
    0.19044202705272903923,
    0.06655998896627697537,
    0.022789258616785717418,
    0.0040116689358507943804,
    0.0012192036851249883877,
];

#[allow(clippy::excessive_precision)]
const FRD: &[Float] = &[
    1.0,
    2.7022305772400260215,
    4.2059268151438492767,
    4.5221882840107715516,
    3.7240352281630359588,
    2.4589286254678152943,
    1.3125491629443702962,
    0.5997685720120932908,
    0.20907680750378849485,
    0.07159621634657901433,
    0.012602969513793714191,
    0.0038302423512931250065,
];

#[allow(clippy::excessive_precision)]
const GN: &[Float] = &[
    0.50000014392706344801,
    0.032346434925349128728,
    0.17619325157863254363,
    0.038606273170706486252,
    0.023693692309257725361,
    0.007092018516845033662,
    0.0012492123212412087428,
    0.00044023040894778468486,
    -8.80266827476172521e-6,
    -1.4033554916580018648e-8,
    2.3509221782155474353e-10,
];

#[allow(clippy::excessive_precision)]
const GD: &[Float] = &[
    1.0,
    2.0646987497019598937,
    2.9109311766948031235,
    2.6561936751333032911,
    2.0195563983177268073,
    1.1167891129189363902,
    0.57267874755973172715,
    0.19408481169593070798,
    0.07634808341431248904,
    0.011573247407207865977,
    0.0044099273693067311209,
    -0.00009070958410429993314,
];

// adapted from ebertolazzi/Clothoids.git Fresnel.cc
//
//

// Compute Fresnel integrals C(x) and S(x)
//
// \f[
//   S(x) = \int_0^x \sin t^2 \,\mathrm{d} t, \qquad
//   C(x) = \int_0^x \cos t^2 \,\mathrm{d} t
// \f]
//
// **Example:**
//
// | \f$ x \f$ | \f$ C(x) \f$  |  \f$ S(x) \f$ |
// | :-------: | :-----------: | :-----------: |
// | 0.0       | 0.00000000    |  0.00000000   |
// | 0.5       | 0.49234423    |  0.06473243   |
// | 1.0       | 0.77989340    |  0.43825915   |
// | 1.5       | 0.44526118    |  0.69750496   |
// | 2.0       | 0.48825341    |  0.34341568   |
// | 2.5       | 0.45741301    |  0.61918176   |
//
// **Adapted from:**
//
// - *William J. Thompson*, Atlas for computing mathematical functions :
//   an illustrated guide for practitioners, with programs in C and Mathematica,
//   Wiley, 1997.
//
// **Author:**
//
// - *Venkata Sivakanth Telasula*,
//   email: sivakanth.telasula@gmail.com,
//   date: August 11, 2005
//
// \param[in]  y Argument of \f$ C(y) \f$ and \f$ S(y) \f$
// \param[out] C \f$ C(x) \f$
// \param[out] S \f$ S(x) \f$

fn fresnel_cs(y: Float) -> (Float, Float) {
    let eps = 1E-15;
    let x = y.abs();

    let mut c_value: Float;
    let mut s_value: Float;

    if x < 1.0 {
        let s = FRAC_PI_2 * (x * x);
        let t = -s * s;

        // Cosine integral series
        {
            let mut twofn = 0.0;
            let mut fact = 1.0;
            let mut denterm = 1.0;
            let mut numterm = 1.0;
            let mut sum: Float = 1.0;
            loop {
                twofn += 2.0;
                fact *= twofn * (twofn - 1.0);
                denterm += 4.0;
                numterm *= t;
                let term = numterm / (fact * denterm);
                sum += term;
                if term.abs() <= eps * sum.abs() {
                    break;
                }
            }

            c_value = x * sum;
        }

        // Sine integral series
        {
            let mut twofn = 1.0;
            let mut fact = 1.0;
            let mut denterm = 3.0;
            let mut numterm = 1.0;
            let mut sum: Float = numterm / denterm;
            loop {
                twofn += 2.0;
                fact *= twofn * (twofn - 1.0);
                denterm += 4.0;
                numterm *= t;
                let term = numterm / (fact * denterm);
                sum += term;
                if term.abs() <= (eps * sum.abs()) {
                    break;
                }
            }

            s_value = FRAC_PI_2 * sum * (x * x * x);
        }
    } else if x < 6.0 {
        // Rational approximation for f
        let f: Float;
        {
            let mut sumn = 0.0;
            let mut sumd = FRD[11];
            for k in (0..=10).rev() {
                sumn = FRN[k] + x * sumn;
                // println!("    sumn = FRN[{}] {} + x {} * sumn {}", FRN[k], k, x, sumn);
                sumd = FRD[k] + x * sumd;
                // println!("    sumd = FRD[{}] {} + x {} * sumd {}", FRD[k], k, x, sumd);
            }
            f = sumn / sumd;
            // println!("  f = sumn {} / sumd {}", sumn, sumd);
        }

        // Rational approximation for g
        let g: Float;
        {
            let mut sumn = 0.0;
            let mut sumd = GD[11];
            for k in (0..=10).rev() {
                sumn = GN[k] + x * sumn;
                sumd = GD[k] + x * sumd;
            }
            g = sumn / sumd;
            // println!("  g = sumn {} / sumd {}", sumn, sumd);
        }

        let u_value = FRAC_PI_2 * (x * x);
        let sin_u = sin(u_value);
        let cos_u = cos(u_value);
        c_value = 0.5 + f * sin_u - g * cos_u;
        s_value = 0.5 - f * cos_u - g * sin_u;

        // println!("  x {}, U {}, f {}, g {}", x, U, f, g);
    } else {
        // x >= 6; asymptotic expansions for  f  and  g

        let s = PI * x * x;
        let t = -1.0 / (s * s);

        // Expansion for f
        let mut numterm = -1.0;
        let mut term = 1.0;
        let mut sum = 1.0;
        // let mut oldterm =  1.0;
        let eps10 = 0.1 * eps;

        loop {
            numterm += 4.0;
            term *= numterm * (numterm - 2.0) * t;
            sum += term;
            let absterm = term.abs();
            /*
            UTILS_ASSERT(
                oldterm >= absterm,
                "In fresnel_cs f not converged to eps, x = {} oldterm = {} absterm = {}\n",
                x, oldterm, absterm
                );
            oldterm = absterm;
            */
            if absterm <= eps10 * sum.abs() {
                break;
            }
        }

        let f = sum / (PI * x);

        //  Expansion for  g
        numterm = -1.0;
        term = 1.0;
        sum = 1.0;
        // oldterm =  1.0;

        loop {
            numterm += 4.0;
            term *= numterm * (numterm + 2.0) * t;
            sum += term;
            let absterm = term.abs();
            /*
            UTILS_ASSERT(
                oldterm >= absterm,
                "In fresnel_cs g not converged to eps, x = {} oldterm = {} absterm = {}\n",
                x, oldterm, absterm
                );
            oldterm = absterm;
            */
            if absterm <= eps10 * sum.abs() {
                break;
            }
        }

        let g0 = PI * x;
        let g = sum / (g0 * g0 * x);

        let u_value = FRAC_PI_2 * (x * x);
        let sin_u = sin(u_value);
        let cos_u = cos(u_value);
        c_value = 0.5 + f * sin_u - g * cos_u;
        s_value = 0.5 - f * cos_u - g * sin_u;
    }
    if y < 0.0 {
        c_value = -c_value;
        s_value = -s_value;
    }

    (c_value, s_value)
}

fn lommel_reduced(mu: Float, nu: Float, b: Float) -> Float {
    let mut tmp = 1.0 / ((mu + nu + 1.0) * (mu - nu + 1.0));
    let mut res = tmp;
    for n in 1..=100 {
        let nf = n as Float;
        tmp *= (-b / (2.0 * nf + mu - nu + 1.0)) * (b / (2.0 * nf + mu + nu + 1.0));
        res += tmp;
        if tmp.abs() < (res.abs() * 1e-50) {
            break;
        }
    }
    res
}

fn eval_xyazero(nk: usize, b: Float) -> ([Float; 43], [Float; 43]) {
    let mut x: [Float; 43] = [0.0; 43];
    let mut y: [Float; 43] = [0.0; 43];
    let sb = sin(b);
    let cb = cos(b);
    let b2 = b * b;
    let threshold = 1e-3;
    if b.abs() < threshold {
        x[0] = 1.0 - (b2 / 6.0) * (1.0 - (b2 / 20.0) * (1.0 - (b2 / 42.0)));
        y[0] = (b / 2.0) * (1.0 - (b2 / 12.0) * (1.0 - (b2 / 30.0)));
    } else {
        x[0] = sb / b;
        y[0] = (1.0 - cb) / b;
    }
    // use recurrence in the stable part
    let mut m = floor(2.0 * b) as usize;
    if m >= nk {
        m = nk - 1;
    }
    if m < 1 {
        m = 1;
    }
    for k in 1..m {
        let kf = k as Float;
        x[k] = (sb - kf * y[k - 1]) / b;
        y[k] = (kf * x[k - 1] - cb) / b;
    }
    //  use Lommel for the unstable part
    if m < nk {
        let a = b * sb;
        let d = sb - b * cb;
        let b = b * d;
        let c = -b2 * sb;
        let m_offset = m as Float + 0.5;
        let mut r_la = lommel_reduced(m_offset, 1.5, b);
        let mut r_ld = lommel_reduced(m_offset, 0.5, b);
        for k in m..nk {
            let kf = k as Float;
            let k_offset = kf + 1.5;
            let r_lb = lommel_reduced(k_offset, 0.5, b);
            let r_lc = lommel_reduced(k_offset, 1.5, b);
            x[k] = (kf * a * r_la + b * r_lb + cb) / (1.0 + kf);
            y[k] = (c * r_lc + sb) / (2.0 + kf) + d * r_ld;
            r_la = r_lc;
            r_ld = r_lb;
        }
    }

    (x, y)
}

fn eval_xy_a_small(a: Float, b: Float, p: usize) -> (Float, Float) {
    // UTILS_ASSERT(
    //   p < 11 && p > 0, "In eval_xy_a_small p = {} must be in 1..10\n", p
    // );

    let nkk = 4 * p + 3; // max 43
    let (x0, y0) = eval_xyazero(nkk, b);

    let mut x = x0[0] - (a / 2.0) * y0[2];
    let mut y = y0[0] + (a / 2.0) * x0[2];

    let mut t = 1.0;
    let aa = -a * a / 4.0; // controllare!
    for n in 1..=p {
        t *= aa / ((2 * n * (2 * n - 1)) as Float);
        let bf = a / ((4 * n + 2) as Float);
        let jj = 4 * n;
        x += t * (x0[jj] - bf * y0[jj + 2]);
        y += t * (y0[jj] + bf * x0[jj + 2]);
    }
    (x, y)
}

fn eval_xy_a_large(a: Float, b: Float) -> (Float, Float) {
    let s = a.signum();
    let absa = a.abs();
    let m_1_sqrt_pi = FRAC_2_SQRT_PI * 0.5;
    let z = m_1_sqrt_pi * sqrt(absa);
    let ell = s * b * m_1_sqrt_pi / sqrt(absa);
    let g = -0.5 * s * (b * b) / absa;
    let cg = cos(g) / z;
    let sg = sin(g) / z;

    // println!("ell {}, z {}", ell, z);
    let (cl, sl) = fresnel_cs(ell);
    let (cz, sz) = fresnel_cs(ell + z);
    // println!("cl {}, sl {}, cz {}, sz {}", cl, sl, cz, sz);

    let d_c0 = cz - cl;
    let d_s0 = sz - sl;

    let x = cg * d_c0 - s * sg * d_s0;
    let y = sg * d_c0 + s * cg * d_s0;

    (x, y)
}

// let (f_c, f_s) = fresnel_cs3(self.curvature_rate() * s * s, self.curvature() * s, self.theta0);
fn fresnel_cs3(a: Float, b: Float, c: Angle) -> (Float, Float) {
    let threshold = 0.01;
    let a_series_size = 3;
    let xx: Float;
    let yy: Float;
    if a.abs() < threshold {
        // println!("small");
        (xx, yy) = eval_xy_a_small(a, b, a_series_size);
    } else {
        // println!("large");
        (xx, yy) = eval_xy_a_large(a, b);
    };

    let cosc = cos(c.get::<radian>());
    let sinc = sin(c.get::<radian>());

    let int_c = xx * cosc - yy * sinc;
    let int_s = xx * sinc + yy * cosc;
    // println!("xx {:0.2}, yy {:0.2}, int_c {:0.2}, int_s {:0.2}", xx, yy, int_c, int_s);

    (int_c, int_s)
}

impl Clothoid {
    pub fn create(
        x0: Length,
        y0: Length,
        theta0: Angle,
        curvature0: Curvature,
        curvature_rate: CurvaturePerLength,
        length: Length,
    ) -> Self {
        Self {
            xy0: Position { x: x0, y: y0 },
            theta0: angle_unwrap(theta0),
            cos_theta0: cos(theta0.get::<radian>()),
            sin_theta0: sin(theta0.get::<radian>()),
            kappa0: curvature0,
            dk: curvature_rate,
            length,
        }
    }

    pub fn get_start_theta(&self) -> Angle {
        self.theta0
    }

    pub fn curvature(&self) -> Curvature {
        self.kappa0
    }

    pub fn set_curvature(&mut self, curvature: Curvature) {
        self.kappa0 = curvature;
    }

    pub fn curvature_rate(&self) -> CurvaturePerLength {
        self.dk
    }

    pub fn set_curvature_rate(&mut self, curvature_rate: CurvaturePerLength) {
        self.dk = curvature_rate;
    }

    pub fn zero_curvature(&self) -> Self {
        let mut straight = self.clone();
        straight.set_curvature(Curvature::zero());
        straight.set_curvature_rate(CurvaturePerLength::zero());
        straight
    }

    // s is length along the curve, x and y will be in same units
    pub fn get_xy(&self, s: Length) -> Position {
        let (f_c, f_s) = fresnel_cs3(
            (self.curvature_rate() * s * s).into(),
            (self.curvature() * s).into(),
            self.theta0,
        );
        // println!("{:0.2} {:0.2} {:0.2}", s, f_c, f_s);
        let x = self.xy0.x + s * f_c;
        let y = self.xy0.y + s * f_s;
        Position { x, y }
    }

    /// get a new Clothoid at this location along the current one
    pub fn get_clothoid(&self, s: Length) -> Self {
        let xy_s = self.get_xy(s);
        // https://github.com/ebertolazzi/Clothoids/blob/master/src/Clothoids/Fresnel.hxx#L142
        // theta(s) = theta + theta' * s + 1/2 * theta0'' * s^2
        let delta_curvature_s: Curvature = (s * self.curvature_rate()).into();
        let theta_s = self.theta0 + Angle::new::<radian>(
            (s * (self.curvature() + 0.5 * delta_curvature_s)).into()
        );
        let kappa_s = self.curvature() + delta_curvature_s; // curvature changes linearly with curvature_rate

        // TODO(lucasw) just use the same length as starting clothoid, or set to self.length - s
        // but only if s < self.length?
        let length = self.length;
        Self {
            xy0: xy_s,
            theta0: angle_unwrap(theta_s),
            cos_theta0: cos(theta_s.get::<radian>()),
            sin_theta0: sin(theta_s.get::<radian>()),
            kappa0: kappa_s,
            dk: self.dk, // curvature rate is constant through the clothoid segment
            length,
        }
    }

    pub fn get_end_clothoid(&self) -> Self {
        self.get_clothoid(self.length)
    }

    pub fn get_points<const NUM: usize>(&self) -> [[Float; 2]; NUM] {
        let mut xys = [[0.0; 2]; NUM];

        let step = self.length / ((NUM - 1) as Float);

        for (i, xys_i) in xys.iter_mut().take(NUM).enumerate() {
            let s: Length = i as Float * step;
            *xys_i = self.get_xy(s).as_array_meter();
        }

        xys
    }

    #[cfg(feature = "std")]
    pub fn get_points_num(&self, mut num: usize) -> Vec<[Float; 2]> {
        let mut xys = Vec::<[Float; 2]>::new();

        if num == 0 {
            num = 1;
        }

        let step = self.length / ((num - 1) as Float);

        for i in 0..num {
            let s: Length = i as Float * step;
            xys.push(self.get_xy(s).as_array_meter());
        }

        xys
    }

    fn dist_sq(p0: &Position, p1: &Position) -> Area {
        let dx = p1.x - p0.x;
        let dy = p1.y - p0.y;
        dx * dx + dy * dy
    }

    /// find the closest/nearest point on the clothoid to the provided point
    /// s0 initial position on curve to start searching from
    /// return the xy location, s distance along curve, distance to point, and number of iterations
    #[cfg(feature = "std")]
    pub fn get_nearest(&self, pt: &Position, s0: Length) -> (Position, Length, Length, usize) {
        // argmin isn't no_std so trying brent_search
        // use rustamath_mnmz::brent_search;
        use rustamath_mnmz::golden_section_search;

        // TODO(lucasw) f32 brent_search?
        let cost_fn = |s: f64| -> f64 {
            let p = self.get_xy(Length::new::<meter>(s as Float));
            Self::dist_sq(pt, &p).get::<square_meter>().into()
        };
        let s0: f64 = s0.get::<meter>().into();
        let max_len: f64 = self.length.get::<meter>().into();
        // TODO(lucasw) the underlying speed of the phenomena
        // causing s to change from the previous value (if that's what is being fed
        // in as an initial value governs what these should be- they could be tighter
        // if the initial value can be a better guess (e.g. a known rate of change
        // is added to the previous solution)
        let bracket_a: f64 = (s0 - 1.0).clamp(0.0, max_len);
        let bracket_b: f64 = (s0 + 1.0).clamp(0.0, max_len);
        let tolerance = 0.001;
        let max_iterations = 30;
        // TODO(lucasw) longer curves that exceed pi radians
        // with more curvature rate don't work as well,  it'll get stuck in a spiral
        // let (s_min, _cost_at_min, iterations) = brent_search(cost_fn, bracket_a, bracket_b, tolerance, max_iterations);
        let (s_min, distance_min, iterations) = golden_section_search(cost_fn, bracket_a, bracket_b, tolerance, max_iterations);

        let s_min = s_min.clamp(0.0, max_len);
        let s = Length::new::<meter>(s_min as Float);
        let distance = Length::new::<meter>(distance_min as Float);
        let p = self.get_xy(s);
        (p, s, distance, iterations)
    }

    /// linear search
    pub fn get_nearest_with_step(&self, pt: &Position, s0: Length, fr: Length, num: usize) -> (Position, Length) {
        let p0 = self.get_xy(s0);
        let dist_sq0 = Self::dist_sq(pt, &p0);
        let mut dist_sq_min = dist_sq0;
        let mut p_min = p0;
        let mut s_min = s0;

        // search a little ways in both directions
        for step in [-fr, fr] {
            let mut s = s0;
            for _ in 0..num {
                s += step;
                if s < Length::zero() || s > self.length {
                    break;
                }
                let p = self.get_xy(s);
                let dist_sq = Self::dist_sq(pt, &p);
                if dist_sq < dist_sq_min {
                    dist_sq_min = dist_sq;
                    p_min = p;
                    s_min = s;
                    // TODO(lucasw) give up if less than tolerance improvement
                } else {
                    break;
                }
            }
        }

        (p_min, s_min)
    }
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::*;
    use std::format;
    use uom::si::curvature::radian_per_meter;

    #[test]
    fn get_points() {
        let curvature = Curvature::new::<radian_per_meter>(1.0);
        let length: Length = (PI / (curvature * 2.0)).into();
        let c0 = Clothoid::create(Length::zero(), Length::zero(), Angle::zero(), curvature, CurvaturePerLength::zero(), length);

        const NUM: usize = 32;
        let pts0 = c0.get_points::<NUM>();
        assert_eq!(pts0.len(), NUM);

        // redundant, std is set for tests
        #[cfg(feature = "std")]
        {
            let pts1 = c0.get_points_num(NUM);
            assert_eq!(pts1.len(), NUM);

            for i in 0..NUM {
                let p0 = pts0[i];
                let p1 = pts1[i];
                assert_eq!(p0[0], p1[0]);
                assert_eq!(p0[1], p1[1]);
            }
        }

        let pt = Position { x: Length::new::<meter>(2.0), y: Length::new::<meter>(1.0) };
        let (pos, s, distance, iterations) = c0.get_nearest(&pt, Length::zero());
        // TODO(lucasw) use assert float eq
        assert!((pos.x.get::<meter>() - 1.0).abs() < 0.001);
        assert!((pos.y.get::<meter>() - 1.0).abs() < 0.001);
        assert!((distance.get::<meter>() - 1.0).abs() < 0.001);
        println!("{pos:?}");
        println!("{s:?}");
        println!("{distance:?}");
        println!("{iterations:?}");
    }

    #[test]
    fn curvatures() {
        // radius = 2.0
        let curvature = Curvature::new::<radian_per_meter>(0.5);
        let length = Length::new::<meter>(1.0);
        let clothoid0 = Clothoid::create(Length::zero(), Length::zero(), Angle::zero(), curvature, CurvaturePerLength::zero(), length);

        // sample the clothoid at the end
        let clothoid1 = clothoid0.get_end_clothoid();
        // curvature_rate should not change
        assert!(clothoid0.curvature_rate() == clothoid1.curvature_rate());

        // with no curvature rate, output clothoid curvature should always match input
        let msg = format!(
            "{:?} -> {:?} at rate {:?}",
            clothoid0.curvature(),
            clothoid1.curvature(),
            clothoid0.curvature_rate()
        );
        assert!(clothoid0.curvature() == clothoid1.curvature(), "{}", msg);
    }
}
