//! Mean orbital elements in the units SGP4 works in.

use crate::time::JulianDate;

/// The inputs to SGP4 initialization, in radians and minutes.
///
/// Build these from a TLE with [`Tle::elements`](crate::Tle::elements), or
/// fill them in directly from another source such as an OMM message.
/// The elements are Kozai mean elements fit with SGP4; osculating elements
/// from other theories will not give correct results.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elements {
    /// Time the elements are valid at (UTC).
    pub epoch: JulianDate,
    /// Drag term, in inverse Earth radii.
    pub bstar: f64,
    /// First time derivative of mean motion, rad/min². SGP4 ignores it.
    pub mean_motion_dot: f64,
    /// Second time derivative of mean motion, rad/min³. SGP4 ignores it.
    pub mean_motion_ddot: f64,
    pub eccentricity: f64,
    /// Argument of perigee, radians.
    pub arg_of_perigee: f64,
    /// Inclination, radians.
    pub inclination: f64,
    /// Mean anomaly, radians.
    pub mean_anomaly: f64,
    /// Kozai mean motion, radians per minute.
    pub mean_motion: f64,
    /// Right ascension of the ascending node, radians.
    pub raan: f64,
}
