use crate::propagation::get_grav_const;

/// Earth gravity model parameters used by SGP4.
///
/// This mirrors python-sgp4's `EarthGravity` namedtuple:
/// (tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EarthGravity {
    pub tumin: f64,
    pub mu: f64,
    pub radiusearthkm: f64,
    pub xke: f64,
    pub j2: f64,
    pub j3: f64,
    pub j4: f64,
    pub j3oj2: f64,
}

impl EarthGravity {
    /// Build an `EarthGravity` from the underlying `getgravconst` helper.
    pub fn from_model(which: &str) -> Self {
        let grav_const = get_grav_const(which);
        EarthGravity {
            tumin: grav_const.tumin,
            mu: grav_const.mu,
            radiusearthkm: grav_const.radiusearthkm,
            xke: grav_const.xke,
            j2: grav_const.j2,
            j3: grav_const.j3,
            j4: grav_const.j4,
            j3oj2: grav_const.j3oj2,
        }
    }
}

/// WGS-72 (old) gravity model parameters.
pub fn wgs72old() -> EarthGravity {
    EarthGravity::from_model("wgs72old")
}

/// WGS-72 gravity model parameters.
pub fn wgs72() -> EarthGravity {
    EarthGravity::from_model("wgs72")
}

/// WGS-84 gravity model parameters.
pub fn wgs84() -> EarthGravity {
    EarthGravity::from_model("wgs84")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "expected {b}, got {a}");
    }

    #[test]
    fn earth_gravity_wrappers_match_getgravconst() {
        for (which, f) in [
            ("wgs72old", wgs72old as fn() -> EarthGravity),
            ("wgs72", wgs72),
            ("wgs84", wgs84),
        ] {
            let eg = f();
            let grav_const = get_grav_const(which);
            approx_eq(eg.tumin, grav_const.tumin, 1e-12);
            approx_eq(eg.mu, grav_const.mu, 1e-9);
            approx_eq(eg.radiusearthkm, grav_const.radiusearthkm, 1e-9);
            approx_eq(eg.xke, grav_const.xke, 1e-12);
            approx_eq(eg.j2, grav_const.j2, 1e-12);
            approx_eq(eg.j3, grav_const.j3, 1e-14);
            approx_eq(eg.j4, grav_const.j4, 1e-14);
            approx_eq(eg.j3oj2, grav_const.j3oj2, 1e-14);
        }
    }
}
