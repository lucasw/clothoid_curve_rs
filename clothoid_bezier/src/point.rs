// use serde::{Deserialize, Serialize};
use core::iter::IntoIterator;
use core::slice;

use super::*;
// use super::Point;

use core::{
    ops::{Add, Mul, Sub},
};

// const PDIM: usize = 2;

/// The Point trait defines constituents Trait requirements and is the only interface on which the library relies.
///
/// The associated constant DIM is necessary so that the memory layout of
/// its implementing type can be made known to the library, whenever new instances are returned.
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
    const DIM: usize;
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
/// TODO(lucasw) the trait `Serialize` is not implemented for `[f64; N]`
#[derive(Debug, Copy, Clone)] // , Deserialize, Serialize)]
pub struct PointN<const N: usize>([NativeFloat; N]);

impl<const N: usize> PointN<N> {
    pub fn new(array: [NativeFloat; N]) -> Self {
        PointN(array)
    }
}

/// Initialize with the Default value for the underlying type
impl<const N: usize> Default for PointN<N> {
    fn default() -> Self {
        PointN([NativeFloat::default(); N])
    }
}

impl<const N: usize> PartialEq for PointN<N>
{
    fn eq(&self, other: &Self) -> bool {
        for i in 0..N {
            if self.0[i] != other.0[i] {
                return false;
            }
        }
        true
    }
}

impl<const N: usize> Add for PointN<N>
{
    type Output = Self;

    fn add(self, other: PointN<N>) -> PointN<N> {
        let mut res = self;
        for i in 0..N {
            res.0[i] = self.0[i] + other.0[i];
        }
        res
    }
}

/// This is not required by the Point trait or library but
/// convenient if you want to use the type externally
impl<const N: usize> Add<NativeFloat> for PointN<N>
{
    type Output = Self;

    fn add(self, _rhs: NativeFloat) -> PointN<N> {
        let mut res = self;
        for i in 0..N {
            res.0[i] = self.0[i] + _rhs;
        }
        res
    }
}

impl<const N: usize> Sub for PointN<N>
{
    type Output = Self;

    fn sub(self, other: PointN<N>) -> PointN<N> {
        let mut res = self;
        for i in 0..N {
            res.0[i] = self.0[i] - other.0[i];
        }
        res
    }
}

/// This is not required by the Point trait or library but
/// convenient if you want to use the type externally
impl<const N: usize> Sub<NativeFloat> for PointN<N>
{
    type Output = Self;

    fn sub(self, _rhs: NativeFloat) -> PointN<N> {
        let mut res = self;
        for i in 0..N {
            res.0[i] = self.0[i] - _rhs;
        }
        res
    }
}

impl<const N: usize> Mul<NativeFloat> for PointN<N>
where
    // The mulitplication is done by mulitpling T * U => T, this
    // trait bound for T will specify this requirement as the mul operator is
    // translated to using the first operand as self and the second as rhs.
{
    type Output = PointN<N>;

    fn mul(self, rhs: NativeFloat) -> PointN<N> {
        let mut res = self;
        for i in 0..res.0.len() {
            res.0[i] *= rhs;
        }
        res
    }
}

impl<const N: usize> IntoIterator for PointN<N> {
    type Item = NativeFloat;
    type IntoIter = core::array::IntoIter<Self::Item, N>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIterator::into_iter(self.0)
    }
}

impl<'a, const N: usize> IntoIterator for &'a mut PointN<N> {
    type Item = &'a mut NativeFloat;
    type IntoIter = slice::IterMut<'a, NativeFloat>;

    fn into_iter(self) -> slice::IterMut<'a, NativeFloat> {
        self.0.iter_mut()
    }
}

impl<const N: usize> Point for PointN<N>
where
{
    const DIM: usize = { N };

    fn axis(&self, index: usize) -> NativeFloat {
        assert!(index <= N);
        self.0[index]
    }

    fn squared_length(&self) -> NativeFloat {
        let mut sqr_dist: NativeFloat = 0.0;
        for i in 0..N {
            sqr_dist += self.0[i] * self.0[i];
        }
        sqr_dist
    }

    fn distance(&self, b: &PointN<N>) -> NativeFloat {
        sqrt((*self - *b).squared_length())
    }
}

#[cfg(test)]
mod tests {
    use super::PointN;
    use super::*;
    #[test]
    fn distance_test() {
        let a = PointN::new([0.0, 0.0]);
        let b = PointN::new([1.0, 1.0]);
        let dist = a.distance(&b);
        // assert!((dist - sqrt(2.0)).abs() < 0.00001);
        assert_eq!(dist, sqrt(2.0));
    }
}
