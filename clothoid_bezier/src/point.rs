use serde::{Deserialize, Serialize};
use core::iter::IntoIterator;
use core::slice;

use super::*;
// use super::Point;

use core::{
    ops::{Add, Mul, Sub},
};

/// The Point trait defines constituents Trait requirements and is the only interface on which the library relies.
///
pub trait Point:
    Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    // + Mul<Self::Scalar, Output = Self>
    + Mul<NativeFloat, Output = Self>
    + Copy
    + PartialEq
    //+ PartialOrd
    + Default
    + IntoIterator
{
    /*
    type Scalar: Float
        + Default
        + PartialEq
        + From<NativeFloat>
        + Into<NativeFloat>
        + Add<NativeFloat, Output = Self::Scalar>
        + Sub<NativeFloat, Output = Self::Scalar>
        + Mul<NativeFloat, Output = Self::Scalar>
        + Div<NativeFloat, Output = Self::Scalar>
        + Sum<NativeFloat>;
    */
    // Returns the component of the Point on its axis corresponding to index e.g. [0, 1, 2] -> [x, y, z]
    /// Panics if index is greater than implementors dimension
    // TODO maybe remove in favour of iterator (?)
    fn axis(&self, index: usize) -> NativeFloat;

    // Returns the squared L2-Norm of the Point interpreted as a Vector
    // TODO this could be moved into the library because computability is ensured by its existing trait bounds
    fn squared_length(&self) -> NativeFloat;

    fn distance(&self, b: &Self) -> NativeFloat;
}

/// Point with dimensions of constant generic size N and of type NatiiveFloat
///
/// (Implemented as Newtype Pattern on an array
/// see book or https://www.worthe-it.co.za/blog/2020-10-31-newtype-pattern-in-rust.html)
/// This type only interacts with the library through
/// the point trait, so you are free to use your own
/// Point/Coord/Vec structures instead by implementing the (small) trait
/// The trait `Serialize` is not implemented for `[f64; N]`
/// so just hard-coding to Point2
#[derive(Debug, Copy, Clone, Deserialize, Serialize)]
pub struct Point2([NativeFloat; 2]);

impl Point2 {
    pub fn new(array: [NativeFloat; 2]) -> Self {
        Point2(array)
    }
}

/// Initialize with the Default value for the underlying type
impl Default for Point2 {
    fn default() -> Self {
        Point2([NativeFloat::default(); 2])
    }
}

impl PartialEq for Point2
{
    fn eq(&self, other: &Self) -> bool {
        for i in 0..2 {
            if self.0[i] != other.0[i] {
                return false;
            }
        }
        true
    }
}

impl Add for Point2
{
    type Output = Self;

    fn add(self, other: Point2) -> Point2 {
        let mut res = self;
        for i in 0..2 {
            res.0[i] = self.0[i] + other.0[i];
        }
        res
    }
}

/// This is not required by the Point trait or library but
/// convenient if you want to use the type externally
impl Add<NativeFloat> for Point2
{
    type Output = Self;

    fn add(self, _rhs: NativeFloat) -> Point2 {
        let mut res = self;
        for i in 0..2 {
            res.0[i] = self.0[i] + _rhs;
        }
        res
    }
}

impl Sub for Point2
{
    type Output = Self;

    fn sub(self, other: Point2) -> Point2 {
        let mut res = self;
        for i in 0..2 {
            res.0[i] = self.0[i] - other.0[i];
        }
        res
    }
}

/// This is not required by the Point trait or library but
/// convenient if you want to use the type externally
impl Sub<NativeFloat> for Point2
{
    type Output = Self;

    fn sub(self, _rhs: NativeFloat) -> Point2 {
        let mut res = self;
        for i in 0..2 {
            res.0[i] = self.0[i] - _rhs;
        }
        res
    }
}

impl Mul<NativeFloat> for Point2
where
    // The mulitplication is done by mulitpling T * U => T, this
    // trait bound for T will specify this requirement as the mul operator is
    // translated to using the first operand as self and the second as rhs.
{
    type Output = Point2;

    fn mul(self, rhs: NativeFloat) -> Point2 {
        let mut res = self;
        for i in 0..res.0.len() {
            res.0[i] *= rhs;
        }
        res
    }
}

impl IntoIterator for Point2 {
    type Item = NativeFloat;
    type IntoIter = core::array::IntoIter<Self::Item, 2>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIterator::into_iter(self.0)
    }
}

impl<'a> IntoIterator for &'a mut Point2 {
    type Item = &'a mut NativeFloat;
    type IntoIter = slice::IterMut<'a, NativeFloat>;

    fn into_iter(self) -> slice::IterMut<'a, NativeFloat> {
        self.0.iter_mut()
    }
}

impl Point for Point2
where
{
    fn axis(&self, index: usize) -> NativeFloat {
        assert!(index <= 2);
        self.0[index]
    }

    fn squared_length(&self) -> NativeFloat {
        let mut sqr_dist: NativeFloat = 0.0;
        for i in 0..2 {
            sqr_dist += self.0[i] * self.0[i];
        }
        sqr_dist
    }

    fn distance(&self, b: &Point2) -> NativeFloat {
        sqrt((*self - *b).squared_length())
    }
}

#[cfg(test)]
mod tests {
    use super::Point2;
    use super::*;
    #[test]
    fn distance_test() {
        let a = Point2::new([0.0, 0.0]);
        let b = Point2::new([1.0, 1.0]);
        let dist = a.distance(&b);
        // assert!((dist - sqrt(2.0)).abs() < 0.00001);
        assert_eq!(dist, sqrt(2.0));
    }
}
