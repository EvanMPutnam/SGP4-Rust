//! Utility routines translated from Vallado's `sgp4ext.cpp`.

use std::f64::consts::PI;

use crate::functions::days2mdhms;

/// Sentinel used by Vallado to indicate undefined values.
pub const UNDEFINED: f64 = 999_999.1;
pub const INFINITE: f64 = 999_999.9;

/// Magnitude of a 3-vector.
pub fn mag(x: &[f64; 3]) -> f64 {
    (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt()
}

/// Cross product `vec1 x vec2`.
pub fn cross(vec1: &[f64; 3], vec2: &[f64; 3]) -> [f64; 3] {
    [
        vec1[1] * vec2[2] - vec1[2] * vec2[1],
        vec1[2] * vec2[0] - vec1[0] * vec2[2],
        vec1[0] * vec2[1] - vec1[1] * vec2[0],
    ]
}

/// Dot product of two 3-vectors.
pub fn dot(x: &[f64; 3], y: &[f64; 3]) -> f64 {
    x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
}

/// Angle between two vectors, or `UNDEFINED` if either is too small.
///
/// Returns an angle in radians in the range `[0, π]` (Python code
/// comments said `-pi to pi`, but the implementation uses `acos`).
pub fn angle(vec1: &[f64; 3], vec2: &[f64; 3]) -> f64 {
    let small = 1.0e-8;

    let magv1 = mag(vec1);
    let magv2 = mag(vec2);

    if magv1 * magv2 > small * small {
        let mut temp = dot(vec1, vec2) / (magv1 * magv2);
        if temp.abs() > 1.0 {
            // clamp to [-1, 1] like Python's copysign logic
            temp = temp.signum().copysign(1.0);
        }
        temp.acos()
    } else {
        UNDEFINED
    }
}

/// Solve Kepler’s equation when true anomaly is known (Vallado `newtonnu`).
///
/// Returns `(e0, m)`:
/// - `e0` – eccentric anomaly (or equivalent parameter)
/// - `m`  – mean anomaly
pub fn newtonnu(ecc: f64, nu: f64) -> (f64, f64) {
    let mut e0 = 999_999.9;
    let mut m = 999_999.9;
    let small = 1.0e-8;

    // Circular
    if ecc.abs() < small {
        m = nu;
        e0 = nu;
    } else if ecc < 1.0 - small {
        // Elliptic
        let sine = ((1.0 - ecc * ecc).sqrt() * nu.sin()) / (1.0 + ecc * nu.cos());
        let cose = (ecc + nu.cos()) / (1.0 + ecc * nu.cos());
        e0 = sine.atan2(cose);
        m = e0 - ecc * e0.sin();
    } else if ecc > 1.0 + small {
        // Hyperbolic
        if ecc > 1.0 && (nu.abs() + 0.00001) < PI - (1.0 / ecc).acos() {
            let sine = ((ecc * ecc - 1.0).sqrt() * nu.sin()) / (1.0 + ecc * nu.cos());
            e0 = sine.asinh();
            m = ecc * e0.sinh() - e0;
        }
    } else if nu.abs() < 168.0 * PI / 180.0 {
        // Parabolic
        e0 = (nu * 0.5).tan();
        m = e0 + (e0 * e0 * e0) / 3.0;
    }

    if ecc < 1.0 {
        // Normalize to [0, 2π)
        let twopi = 2.0 * PI;
        m = m % twopi;
        if m < 0.0 {
            m += twopi;
        }
        e0 = e0 % twopi;
    }

    (e0, m)
}

/// Type of orbit used internally in `rv2coe`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OrbitType {
    /// Elliptical / parabolic / hyperbolic inclined (`"ei"` in Python).
    Ei,
    /// Circular equatorial (`"ce"`).
    Ce,
    /// Circular inclined (`"ci"`).
    Ci,
    /// Elliptical equatorial (`"ee"`).
    Ee,
}

/// Convert position/velocity vectors to classical orbital elements (`rv2coe`).
///
/// Returns:
/// `(p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper)`.
pub fn rv2coe(
    r: &[f64; 3],
    v: &[f64; 3],
    mu: f64,
) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let mut hbar = [0.0_f64; 3];
    let mut nbar = [0.0_f64; 3];
    let mut ebar = [0.0_f64; 3];

    let twopi = 2.0 * PI;
    let halfpi = 0.5 * PI;
    let small = 1.0e-8;

    let magr = mag(r);
    let magv = mag(v);

    // cross(r, v, hbar)
    hbar = cross(r, v);
    let magh = mag(&hbar);

    if magh <= small {
        // Undefined orbit
        return (
            UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED,
            UNDEFINED, UNDEFINED, UNDEFINED,
        );
    }

    // nbar = [-h_y, h_x, 0]
    nbar[0] = -hbar[1];
    nbar[1] = hbar[0];
    nbar[2] = 0.0;

    let magn = mag(&nbar);
    let c1 = magv * magv - mu / magr;
    let rdotv = dot(r, v);

    for i in 0..3 {
        ebar[i] = (c1 * r[i] - rdotv * v[i]) / mu;
    }
    let ecc = mag(&ebar);

    // semi-major axis / semi-latus rectum
    let sme = 0.5 * magv * magv - mu / magr;
    let a = if sme.abs() > small {
        -mu / (2.0 * sme)
    } else {
        INFINITE
    };
    let p = magh * magh / mu;

    // inclination
    let hk = hbar[2] / magh;
    let incl = hk.acos();

    // Determine orbit type
    let mut typeorbit = OrbitType::Ei;
    if ecc < small {
        // circular
        if incl < small || (incl - PI).abs() < small {
            typeorbit = OrbitType::Ce;
        } else {
            typeorbit = OrbitType::Ci;
        }
    } else if incl < small || (incl - PI).abs() < small {
        typeorbit = OrbitType::Ee;
    }

    // Longitude of ascending node, omega
    let omega = if magn > small {
        let mut temp = nbar[0] / magn;
        if temp.abs() > 1.0 {
            temp = temp.signum().copysign(1.0);
        }
        let mut omega = temp.acos();
        if nbar[1] < 0.0 {
            omega = twopi - omega;
        }
        omega
    } else {
        UNDEFINED
    };

    // Argument of perigee
    let argp = if typeorbit == OrbitType::Ei {
        let mut argp = angle(&nbar, &ebar);
        if ebar[2] < 0.0 {
            argp = twopi - argp;
        }
        argp
    } else {
        UNDEFINED
    };

    // True anomaly
    let nu = if matches!(typeorbit, OrbitType::Ei | OrbitType::Ee) {
        let mut nu = angle(&ebar, r);
        if rdotv < 0.0 {
            nu = twopi - nu;
        }
        nu
    } else {
        UNDEFINED
    };

    // Argument of latitude (circular inclined)
    let mut arglat = UNDEFINED;
    let mut m = UNDEFINED;
    if typeorbit == OrbitType::Ci {
        arglat = angle(&nbar, r);
        if r[2] < 0.0 {
            arglat = twopi - arglat;
        }
        m = arglat;
    }

    // Longitude of periapsis (elliptical equatorial)
    let lonper = if ecc > small && typeorbit == OrbitType::Ee {
        let mut temp = ebar[0] / ecc;
        if temp.abs() > 1.0 {
            temp = temp.signum().copysign(1.0);
        }
        let mut lonper = temp.acos();
        if ebar[1] < 0.0 {
            lonper = twopi - lonper;
        }
        if incl > halfpi {
            lonper = twopi - lonper;
        }
        lonper
    } else {
        UNDEFINED
    };

    // True longitude (circular equatorial)
    let truelon = if magr > small && typeorbit == OrbitType::Ce {
        let mut temp = r[0] / magr;
        if temp.abs() > 1.0 {
            temp = temp.signum().copysign(1.0);
        }
        let mut truelon = temp.acos();
        if r[1] < 0.0 {
            truelon = twopi - truelon;
        }
        if incl > halfpi {
            truelon = twopi - truelon;
        }
        m = truelon;
        truelon
    } else {
        UNDEFINED
    };

    // Mean anomaly for elliptical cases
    if matches!(typeorbit, OrbitType::Ei | OrbitType::Ee) {
        let (_e, m_out) = newtonnu(ecc, nu);
        m = m_out;
    }

    (p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper)
}

/// Vallado-style Julian Date (single `f64`), not the split `(jd, fr)` version.
///
/// This is *different* from `functions::jday` in your port, which returns
/// `(jd, fr)`. This corresponds to the `jday()` at the bottom of `sgp4ext.cpp`.
pub fn jday(year: i32, mon: i32, day: i32, hr: i32, minute: i32, sec: f64) -> f64 {
    // faithful translation of the Python / Vallado formula
    let year_f = year as f64;
    let mon_f = mon as f64;
    let day_f = day as f64;
    let hr_f = hr as f64;
    let min_f = minute as f64;

    367.0 * year_f - (7.0 * (year_f + ((mon_f + 9.0) / 12.0).floor()) * 0.25).floor()
        + (275.0 * mon_f / 9.0).floor()
        + day_f
        + 1_721_013.5
        + (((sec / 60.0 + min_f) / 60.0) + hr_f) / 24.0
}

/// Inverse of `jday` — convert Julian Date to calendar date/time (`invjday`).
///
/// Returns `(year, month, day, hour, minute, second)`.
///
/// This uses `functions::days2mdhms` from your library, like the Python code.
pub fn invjday(jd: f64) -> (i32, i32, i32, i32, i32, f64) {
    // --------------- find year and days of the year ---------------
    let temp = jd - 2_415_019.5;
    let tu = temp / 365.25;
    let mut year = 1900 + tu.floor() as i32;
    let mut leapyrs = (((year - 1901) as f64 * 0.25).floor()) as i32;

    // optional nudge by 8.64e-7 sec to get even outputs
    let mut days = temp - ((year - 1900) as f64 * 365.0 + leapyrs as f64) + 0.00000000001;

    // check for case of beginning of a year
    if days < 1.0 {
        year -= 1;
        leapyrs = (((year - 1901) as f64 * 0.25).floor()) as i32;
        days = temp - ((year - 1900) as f64 * 365.0 + leapyrs as f64);
    }

    // use your existing days2mdhms from `functions`
    let (mon, day, hr, minute, mut sec) = days2mdhms(year, days, false);
    sec -= 0.00000086400;

    // Cast everything to i32 for a clean API
    (year, mon as i32, day as i32, hr as i32, minute as i32, sec)
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn approx(actual: f64, expected: f64, tol: f64) {
        let diff = (actual - expected).abs();
        assert!(
            diff <= tol,
            "expected {expected}, got {actual} (|Δ| = {diff})"
        );
    }

    #[test]
    fn mag_basic() {
        let v = [3.0_f64, 4.0, 12.0];
        let m = mag(&v);
        approx(m, 13.0, 1e-12);
    }

    #[test]
    fn dot_basic() {
        let a = [1.0_f64, 2.0, 3.0];
        let b = [4.0_f64, -5.0, 6.0];
        let d = dot(&a, &b);
        approx(d, 12.0, 1e-12); // 1*4 + 2*(-5) + 3*6 = 12
    }

    #[test]
    fn cross_basic() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0_f64, 1.0, 0.0];
        let out = cross(&a, &b);
        approx(out[0], 0.0, 1e-12);
        approx(out[1], 0.0, 1e-12);
        approx(out[2], 1.0, 1e-12);
    }

    #[test]
    fn angle_orthogonal() {
        let x = [1.0_f64, 0.0, 0.0];
        let y = [0.0_f64, 1.0, 0.0];
        let theta = angle(&x, &y);
        approx(theta, PI / 2.0, 1e-12);
    }

    #[test]
    fn newtonnu_circular_degenerates_to_nu() {
        let ecc = 0.0_f64;
        let nu = 1.2345_f64;
        let (e0, m) = newtonnu(ecc, nu);
        approx(e0, nu, 1e-12);
        approx(m, nu, 1e-12);
    }

    /// Elliptic example — values taken from a direct run of Vallado's algorithm.
    #[test]
    fn newtonnu_elliptic_example() {
        let ecc = 0.1_f64;
        let nu = 1.0_f64; // radians

        let (e0, m) = newtonnu(ecc, nu);

        // Reference values from the Python/Vallado implementation.
        let expected_e0 = 0.917_912_038_666_774_1_f64;
        let expected_m = 0.838_478_542_901_907_3_f64;

        approx(e0, expected_e0, 1e-12);
        approx(m, expected_m, 1e-12);
    }

    /// Hyperbolic example — values taken from a direct run of Vallado's algorithm.
    #[test]
    fn newtonnu_hyperbolic_example() {
        let ecc = 1.5_f64;
        let nu = 0.5_f64; // radians

        let (e0, m) = newtonnu(ecc, nu);

        // Reference values from the Python/Vallado implementation.
        let expected_e0 = 0.229_385_302_037_439_15_f64;
        let expected_m = 0.117_718_026_442_174_41_f64;

        approx(e0, expected_e0, 1e-12);
        approx(m, expected_m, 1e-12);
    }

    #[test]
    fn jday_and_invjday_round_trip() {
        let year = 2020;
        let mon = 2;
        let day = 29;
        let hr = 23;
        let minute = 59;
        let sec = 30.5;

        let jd = jday(year, mon, day, hr, minute, sec);
        let (yy, mm, dd, hh, min, s) = invjday(jd);

        assert_eq!(yy, year);
        assert_eq!(mm, mon);
        assert_eq!(dd, day);
        assert_eq!(hh, hr);
        assert_eq!(min, minute);
        approx(s, sec, 1e-4);
    }

    #[test]
    fn rv2coe_elliptic_inclined_example() {
        use crate::ext::rv2coe;

        // Elliptic, inclined, non-equatorial orbit ("ei" case)
        // a = 26560 km, e = 0.1, i = 55 deg, RAAN = 40 deg,
        // argp = 30 deg, nu = 10 deg, mu = 398600.4418 km^3/s^2.
        // r, v computed from those COEs in perifocal frame then rotated to ECI.
        let r = [
            8374.048172997831_f64,
            18547.275009889749_f64,
            12603.838242542852_f64,
        ];
        let v = [-3.290244896520_f64, -0.275171108848_f64, 2.719387109636_f64];
        let mu = 398_600.4418_f64;

        let (p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper) = rv2coe(&r, &v, mu);

        // Use slightly relaxed tolerances because we don't care about
        // sub-nanometer precision here, just that the implementation matches
        // the Vallado-style math.
        let tol_km = 1.0e-2;
        let tol_ecc = 1.0e-8;
        let tol_ang = 1.0e-8;

        // Expected values from the COEs we started with
        approx(p, 26294.4, tol_km);
        approx(a, 26560.0, tol_km);
        approx(ecc, 0.1, tol_ecc);

        // Angles in radians: 55°, 40°, 30°, 10°, and mean anomaly ~0.142215 rad
        let deg_to_rad = std::f64::consts::PI / 180.0;
        approx(incl, 55.0 * deg_to_rad, tol_ang);
        approx(omega, 40.0 * deg_to_rad, tol_ang);
        approx(argp, 30.0 * deg_to_rad, tol_ang);
        approx(nu, 10.0 * deg_to_rad, tol_ang);
        approx(m, 0.142215, 1.0e-6);

        // We *expect* these to be undefined for an "ei" orbit; you can either
        // ignore them or assert they hit the sentinel.
        assert!((arglat - 999_999.1).abs() < 1.0e-6);
        assert!((truelon - 999_999.1).abs() < 1.0e-6);
        assert!((lonper - 999_999.1).abs() < 1.0e-6);
    }
}
