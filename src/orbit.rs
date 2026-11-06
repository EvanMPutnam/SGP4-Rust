//! Classical (osculating) orbital elements from a position and velocity.
//!
//! Ports Vallado's `rv2coe` and `newtonnu`. Values that are undefined for an
//! orbit's shape, such as the argument of perigee of a circular orbit, are
//! `None` instead of Vallado's 999999.1 sentinel.

use std::f64::consts::{PI, TAU};

/// Below this, eccentricity and inclination count as zero.
const SMALL: f64 = 1.0e-8;

/// Osculating orbital elements. Angles are in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassicalElements {
    /// Semi-latus rectum, km.
    pub semi_latus_rectum: f64,
    /// Semi-major axis, km. Negative for hyperbolic orbits, `None` for parabolic.
    pub semi_major_axis: Option<f64>,
    pub eccentricity: f64,
    pub inclination: f64,
    /// Right ascension of the ascending node. `None` for equatorial orbits.
    pub raan: Option<f64>,
    /// `None` for circular or equatorial orbits.
    pub arg_of_perigee: Option<f64>,
    /// `None` for circular orbits.
    pub true_anomaly: Option<f64>,
    /// For circular orbits, the argument of latitude or true longitude stands in.
    pub mean_anomaly: Option<f64>,
    /// Circular inclined orbits only.
    pub arg_of_latitude: Option<f64>,
    /// Circular equatorial orbits only.
    pub true_longitude: Option<f64>,
    /// Elliptical equatorial orbits only.
    pub longitude_of_periapsis: Option<f64>,
}

impl ClassicalElements {
    /// Computes elements from a position (km) and velocity (km/s) with
    /// gravitational parameter `mu` (km³/s², see [`GravityConstants::mu`](crate::GravityConstants::mu)).
    ///
    /// Returns `None` when the angular momentum is near zero, as for
    /// straight-line motion through the center of the Earth.
    pub fn from_state(
        position: [f64; 3],
        velocity: [f64; 3],
        mu: f64,
    ) -> Option<ClassicalElements> {
        let r = position;
        let v = velocity;
        let magr = mag(r);
        let magv = mag(v);
        let hbar = cross(r, v);
        let magh = mag(hbar);
        if magh <= SMALL {
            return None;
        }

        let nbar = [-hbar[1], hbar[0], 0.0];
        let magn = mag(nbar);
        let c1 = magv * magv - mu / magr;
        let rdotv = dot(r, v);
        let ebar = [0, 1, 2].map(|i| (c1 * r[i] - rdotv * v[i]) / mu);
        let ecc = mag(ebar);

        let sme = 0.5 * magv * magv - mu / magr;
        let semi_major_axis = (sme.abs() > SMALL).then(|| -mu / (2.0 * sme));
        let p = magh * magh / mu;
        let incl = (hbar[2] / magh).acos();

        let equatorial = incl < SMALL || (incl - PI).abs() < SMALL;
        let circular = ecc < SMALL;

        let raan = (magn > SMALL).then(|| {
            let omega = clamp_unit(nbar[0] / magn).acos();
            if nbar[1] < 0.0 { TAU - omega } else { omega }
        });

        let arg_of_perigee = (!circular && !equatorial).then(|| {
            let argp = angle(nbar, ebar);
            if ebar[2] < 0.0 { TAU - argp } else { argp }
        });

        let true_anomaly = (!circular).then(|| {
            let nu = angle(ebar, r);
            if rdotv < 0.0 { TAU - nu } else { nu }
        });

        let arg_of_latitude = (circular && !equatorial).then(|| {
            let arglat = angle(nbar, r);
            if r[2] < 0.0 { TAU - arglat } else { arglat }
        });

        let longitude_of_periapsis = (!circular && equatorial).then(|| {
            let lonper = clamp_unit(ebar[0] / ecc).acos();
            let lonper = if ebar[1] < 0.0 { TAU - lonper } else { lonper };
            if incl > PI / 2.0 {
                TAU - lonper
            } else {
                lonper
            }
        });

        let true_longitude = (circular && equatorial && magr > SMALL).then(|| {
            let truelon = clamp_unit(r[0] / magr).acos();
            let truelon = if r[1] < 0.0 { TAU - truelon } else { truelon };
            if incl > PI / 2.0 {
                TAU - truelon
            } else {
                truelon
            }
        });

        let mean_anomaly = match true_anomaly {
            Some(nu) => anomalies_from_true(ecc, nu).map(|a| a.mean),
            None => arg_of_latitude.or(true_longitude),
        };

        Some(ClassicalElements {
            semi_latus_rectum: p,
            semi_major_axis,
            eccentricity: ecc,
            inclination: incl,
            raan,
            arg_of_perigee,
            true_anomaly,
            mean_anomaly,
            arg_of_latitude,
            true_longitude,
            longitude_of_periapsis,
        })
    }
}

/// The anomalies that correspond to a true anomaly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Anomalies {
    /// Eccentric anomaly; hyperbolic or parabolic anomaly for open orbits.
    pub eccentric: f64,
    pub mean: f64,
}

/// Converts a true anomaly to eccentric and mean anomaly (Vallado's `newtonnu`).
///
/// Returns `None` for a hyperbolic orbit past its asymptote, or a parabolic
/// orbit beyond 168 degrees, where Vallado's method has no answer.
/// For closed orbits the mean anomaly is in `[0, 2π)`.
pub fn anomalies_from_true(eccentricity: f64, true_anomaly: f64) -> Option<Anomalies> {
    let ecc = eccentricity;
    let nu = true_anomaly;
    let (e0, m) = if ecc.abs() < SMALL {
        (nu, nu)
    } else if ecc < 1.0 - SMALL {
        let sine = ((1.0 - ecc * ecc).sqrt() * nu.sin()) / (1.0 + ecc * nu.cos());
        let cose = (ecc + nu.cos()) / (1.0 + ecc * nu.cos());
        let e0 = sine.atan2(cose);
        (e0, e0 - ecc * e0.sin())
    } else if ecc > 1.0 + SMALL {
        if nu.abs() + 0.00001 >= PI - (1.0 / ecc).acos() {
            return None;
        }
        let sine = ((ecc * ecc - 1.0).sqrt() * nu.sin()) / (1.0 + ecc * nu.cos());
        let e0 = sine.asinh();
        (e0, ecc * e0.sinh() - e0)
    } else if nu.abs() < 168.0_f64.to_radians() {
        let e0 = (nu * 0.5).tan();
        (e0, e0 + (e0 * e0 * e0) / 3.0)
    } else {
        return None;
    };

    if ecc < 1.0 {
        Some(Anomalies {
            eccentric: e0 % TAU,
            mean: m.rem_euclid(TAU),
        })
    } else {
        Some(Anomalies {
            eccentric: e0,
            mean: m,
        })
    }
}

fn mag(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Angle between two non-zero vectors, in `[0, π]`.
fn angle(a: [f64; 3], b: [f64; 3]) -> f64 {
    clamp_unit(dot(a, b) / (mag(a) * mag(b))).acos()
}

/// Clamps a cosine that rounding pushed just outside `[-1, 1]`.
fn clamp_unit(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn elliptic_inclined_orbit() {
        // a = 26560 km, e = 0.1, i = 55°, RAAN = 40°, argp = 30°, nu = 10°.
        let r = [8374.048172997831, 18547.27500988975, 12603.838242542852];
        let v = [-3.290244896520, -0.275171108848, 2.719387109636];
        let coe = ClassicalElements::from_state(r, v, 398_600.441_8).unwrap();

        approx(coe.semi_latus_rectum, 26_294.4, 1e-2);
        approx(coe.semi_major_axis.unwrap(), 26_560.0, 1e-2);
        approx(coe.eccentricity, 0.1, 1e-8);
        approx(coe.inclination, 55f64.to_radians(), 1e-8);
        approx(coe.raan.unwrap(), 40f64.to_radians(), 1e-8);
        approx(coe.arg_of_perigee.unwrap(), 30f64.to_radians(), 1e-8);
        approx(coe.true_anomaly.unwrap(), 10f64.to_radians(), 1e-8);
        approx(coe.mean_anomaly.unwrap(), 0.142215, 1e-6);
        assert_eq!(coe.arg_of_latitude, None);
        assert_eq!(coe.true_longitude, None);
        assert_eq!(coe.longitude_of_periapsis, None);
    }

    #[test]
    fn circular_equatorial_orbit_uses_true_longitude() {
        let mu = 398_600.8;
        let r = [0.0, -7000.0, 0.0];
        let v = [(mu / 7000.0_f64).sqrt(), 0.0, 0.0];
        let coe = ClassicalElements::from_state(r, v, mu).unwrap();
        assert!(coe.eccentricity < 1e-8);
        assert_eq!(coe.raan, None);
        assert_eq!(coe.arg_of_perigee, None);
        assert_eq!(coe.true_anomaly, None);
        approx(coe.true_longitude.unwrap(), 1.5 * PI, 1e-12);
        assert_eq!(coe.mean_anomaly, coe.true_longitude);
    }

    #[test]
    fn rectilinear_motion_has_no_elements() {
        assert_eq!(
            ClassicalElements::from_state([7000.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 398_600.8),
            None
        );
    }

    #[test]
    fn anomalies_match_python_sgp4() {
        let a = anomalies_from_true(0.1, 1.0).unwrap();
        approx(a.eccentric, 0.917_912_038_666_774_1, 1e-12);
        approx(a.mean, 0.838_478_542_901_907_3, 1e-12);

        let a = anomalies_from_true(0.1, 40f64.to_radians()).unwrap();
        approx(a.eccentric, 0.636_094_245_730_360_9, 1e-12);
        approx(a.mean, 0.576_688_435_198_464_7, 1e-12);

        let parabolic = anomalies_from_true(1.0, 2.9).unwrap();
        approx(parabolic.eccentric, 8.238_092_752_965_605, 1e-12);
        approx(parabolic.mean, 194.600_699_894_828_98, 1e-12);

        let hyperbolic = anomalies_from_true(1.1, 2.7).unwrap();
        approx(hyperbolic.eccentric, 4.262_200_676_156_417, 1e-12);
        approx(hyperbolic.mean, 34.761_340_820_283_72, 1e-12);
        let hyperbolic = anomalies_from_true(1.5, 0.5).unwrap();
        approx(hyperbolic.eccentric, 0.229_385_302_037_439_15, 1e-12);
        approx(hyperbolic.mean, 0.117_718_026_442_174_41, 1e-12);
    }

    #[test]
    fn anomalies_undefined_past_the_asymptote() {
        assert_eq!(anomalies_from_true(1.5, 2.5), None);
        assert_eq!(anomalies_from_true(1.0, 3.0), None);
        assert_eq!(
            anomalies_from_true(0.0, 1.2345).unwrap(),
            Anomalies {
                eccentric: 1.2345,
                mean: 1.2345
            }
        );
    }
}
