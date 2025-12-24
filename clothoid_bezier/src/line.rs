use super::Point;
use super::*;

/// LineSegment defined by a start and an endpoint, evaluable anywhere inbetween using interpolation parameter t: [0,1] in eval()
///
/// A LineSegment is equal to a linear Bezier curve, which is why there is no
/// specialized type for that case.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LineSegment<P, const PDIM: usize> {
    pub(crate) start: P,
    pub(crate) end: P,
}

impl<P, const PDIM: usize> LineSegment<P, PDIM>
where
    P: Point,
{
    pub fn new(start: P, end: P) -> Self {
        LineSegment { start, end }
    }

    pub fn eval(&self, t: NativeFloat) -> P {
        self.start + (self.end - self.start) * t
    }

    pub fn split(&self, t: NativeFloat) -> (Self, Self) {
        // compute the split point by interpolation
        let ctrl_ab = self.start + (self.start - self.end) * t;

        (
            LineSegment {
                start: self.start,
                end: ctrl_ab,
            },
            LineSegment {
                start: ctrl_ab,
                end: self.end,
            },
        )
    }

    // DEPRECATED
    // pub fn to_line(&self) -> Line<P> {
    //     Line {
    //         origin: self.start,
    //         vector: self.end - self.start,
    //     }
    // }

    /// Return the distance from the LineSegment to Point p by calculating the projection
    pub fn distance_to_point(&self, p: P) -> NativeFloat {
        let l2 = (self.end - self.start).squared_length();
        // if start and endpoint are approx the same, return the distance to either
        if l2 < NativeFloat::from(EPSILON) {
            sqrt((self.start - p).squared_length())
        } else {
            let v1 = p - self.start;
            let v2 = self.end - self.start;
            let mut dot = NativeFloat::from(0.0);
            for (i, _) in v1.into_iter().enumerate() {
                dot += v1.axis(i) * v2.axis(i);
            }
            // v1 and v2 will by definition always have the same number of axes and produce a value for each Item
            // dot = v1.into_iter()
            //         .zip(v2.into_iter())
            //         .map(|(x1, x2)| x1 * x2)
            //         .sum::<NativeFloat>();
            let mut t = NativeFloat::from(0.0);
            if dot / l2 < NativeFloat::from(1.0) {
                t = dot / l2;
            }
            if t < NativeFloat::from(0.0) {
                t = NativeFloat::from(0.0);
            }
            let projection = self.start + (self.end - self.start) * t; // Projection falls on the segment

            sqrt((p - projection).squared_length())
        }
    }

    /// Sample the coordinate axis of the segment at t (expecting t between 0 and 1).
    pub fn axis(&self, t: NativeFloat, axis: usize) -> NativeFloat {
        self.start.axis(axis) + (self.end.axis(axis) - self.start.axis(axis)) * t
    }

    /// Return the derivative function.
    /// The derivative is also a bezier curve but of degree n-1 - In the case of a line the derivative vector.
    pub fn derivative(&self) -> P {
        self.end - self.start
    }

    pub(crate) fn root(&self, a: NativeFloat, b: NativeFloat) -> ArrayVec<[NativeFloat; 1]> {
        let mut r = ArrayVec::new();
        if a.abs() < EPSILON {
            return r;
        }
        r.push(-b / a);
        r
    }

    /// Return the bounding box of the line as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box(&self) -> [(NativeFloat, NativeFloat); PDIM] {
        let mut bounds = [(NativeFloat::default(), NativeFloat::default()); PDIM];

        // find min/max for that particular axis
        // TODO shoul be rewritten once 'Iterator' is implemented on P to get rid of .axis() method
        for (i, _) in self.start.into_iter().enumerate() {
            if self.start.axis(i) < self.end.axis(i) {
                bounds[i] = (self.start.axis(i), self.end.axis(i));
            } else {
                bounds[i] = (self.end.axis(i), self.start.axis(i));
            }
        }

        bounds
    }
}

#[cfg(test)]
mod tests {
    use super::point::PointN;
    use super::*;
    /// Check whether a line segment interpolation p + t*(q-p) at t=0.5
    /// yields equal distance to the start (p)/end (q) points (up to machine accuracy).
    #[test]
    fn line_segment_interpolation() {
        let start = PointN::new([0.0, 1.77]);
        let end = PointN::new([4.3, 3.0]);
        let line = LineSegment::<_, 2> {
            start,
            end,
        };

        let mid = line.eval(0.5);
        assert!((mid - line.start).squared_length() - (mid - line.end).squared_length() < EPSILON)
    }

    /// Check whether classic pythagorean equality holds for sides 3, 4 with hypothenuse 5
    #[test]
    fn line_segment_distance_to_point() {
        // 3D cause why not
        let line = LineSegment::<_, 2> {
            start: PointN::new([0.0, 1.0, 0.0]),
            end: PointN::new([3.0, 1.0, 0.0]),
        };
        // dist to start should be 4; dist to end should be 5
        let p1 = PointN::new([0.0, 5.0, 0.0]);
        assert!((line.distance_to_point(p1) - 4.0).abs() < EPSILON);
        assert!(((p1 - line.start).squared_length().sqrt() - 4.0).abs() < EPSILON);
        assert!(((p1 - line.end).squared_length().sqrt() - 5.0).abs() < EPSILON);
        // dist to midpoint (t=0.5) should be 1
        let p2 = PointN::new([1.5, 2.0, 0.0]);
        assert!(((p2 - line.eval(0.5)).squared_length().sqrt() - 1.0).abs() < EPSILON);
    }
}
