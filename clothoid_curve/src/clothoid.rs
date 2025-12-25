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

/// put angle into -pi, pi range
pub fn angle_unwrap(angle_radians: Float) -> Float {
    (angle_radians + PI) % (2.0 * PI) - PI
}

#[cfg(feature = "std")]
extern crate alloc;
#[cfg(feature = "std")]
use alloc::vec::Vec;

#[derive(Clone, PartialEq)]
pub struct Clothoid {
    pub x0: Float,     // start point x
    pub y0: Float,     // start point y
    pub theta0: Float, // start point theta/yaw/heading
    kappa0: Float,     // start point curvature 1/r
    dk: Float, // curvature rate, how much curvature changes per unit length, end theta will be
    // theta(s) = theta + theta' * s + 1/2 * theta0'' * s^2
    // theta(s) = theta + kappa0 * s + 1/2 * dk * s^2
    // with s = length
    pub length: Float, // how long the curve is (end kappa will be length * dk)
}

impl Default for Clothoid {
    fn default() -> Self {
        Self {
            x0: 0.0,
            y0: 0.0,
            theta0: 0.0,
            kappa0: 0.0,
            dk: 0.0,
            length: 1.0, // TODO(lucasw) should this be 0.0?
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

fn fresnel_cs3(a: Float, b: Float, c: Float) -> (Float, Float) {
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

    let cosc = cos(c);
    let sinc = sin(c);

    let int_c = xx * cosc - yy * sinc;
    let int_s = xx * sinc + yy * cosc;
    // println!("xx {:0.2}, yy {:0.2}, int_c {:0.2}, int_s {:0.2}", xx, yy, int_c, int_s);

    (int_c, int_s)
}

impl Clothoid {
    pub fn create(
        x0: Float,
        y0: Float,
        theta0: Float,
        kappa0: Float,
        dk: Float,
        length: Float,
    ) -> Self {
        Self {
            x0,
            y0,
            theta0: angle_unwrap(theta0),
            kappa0,
            dk,
            length,
        }
    }

    pub fn get_start_theta(&self) -> Float {
        self.theta0
    }

    pub fn curvature(&self) -> Float {
        self.kappa0
    }

    pub fn curvature_rate(&self) -> Float {
        self.dk
    }

    // s is length along the curve, x and y will be in same units
    fn get_xy(&self, s: Float) -> (Float, Float) {
        let (f_c, f_s) = fresnel_cs3(self.dk * s * s, self.kappa0 * s, self.theta0);
        // println!("{:0.2} {:0.2} {:0.2}", s, f_c, f_s);
        let x = self.x0 + s * f_c;
        let y = self.y0 + s * f_s;
        (x, y)
    }

    pub fn get_xy_array(&self, s: Float) -> [Float; 2] {
        let (x, y) = self.get_xy(s);
        [x, y]
    }

    /// get a new Clothoid at this location along the current one
    pub fn get_clothoid(&self, s: Float) -> Self {
        let (x_s, y_s) = self.get_xy(s);
        // https://github.com/ebertolazzi/Clothoids/blob/master/src/Clothoids/Fresnel.hxx#L142
        // theta(s) = theta + theta' * s + 1/2 * theta0'' * s^2
        let theta_s = self.theta0 + s * (self.kappa0 + 0.5 * s * self.dk);
        let kappa_s = self.kappa0 + s * self.dk; // curvature changes linearly with curvature_rate

        // TODO(lucasw) just use the same length as starting clothoid, or set to self.length - s
        // but only if s < self.length?
        let length = self.length;
        Self {
            x0: x_s,
            y0: y_s,
            theta0: angle_unwrap(theta_s),
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
        let mut s = 0.0;

        for xys_i in xys.iter_mut().take(NUM) {
            let xy = self.get_xy(s);
            // println!("{:0.2} {:?}", s, xy);
            *xys_i = [xy.0, xy.1];
            s += step;
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
        let mut s = 0.0;

        for xys_i in xys.iter_mut().take(num) {
            let xy = self.get_xy(s);
            // println!("{:0.2} {:?}", s, xy);
            *xys_i = [xy.0, xy.1];
            s += step;
        }

        xys
    }
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::*;
    use std::format;

    #[test]
    fn curvatures() {
        // radius = 2.0
        let curvature = 0.5;
        let length = 1.0;
        let clothoid0 = Clothoid::create(0.0, 0.0, 0.0, curvature, 0.0, length);

        // sample the clothoid at the end
        let clothoid1 = clothoid0.get_end_clothoid();
        // curvature_rate should not change
        assert!(clothoid0.curvature_rate() == clothoid1.curvature_rate());

        // with no curvature rate, output clothoid curvature should always match input
        let msg = format!(
            "{} -> {} at rate {}",
            clothoid0.curvature(),
            clothoid1.curvature(),
            clothoid0.curvature_rate()
        );
        assert!(clothoid0.curvature() == clothoid1.curvature(), "{}", msg);
    }
}
