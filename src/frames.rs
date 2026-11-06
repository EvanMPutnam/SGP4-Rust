//! Earth-fixed and geodetic coordinates for SGP4 output.

use crate::time::JulianDate;

/// WGS-84 equatorial radius, km.
const WGS84_A: f64 = 6378.137;
/// WGS-84 flattening.
const WGS84_F: f64 = 1.0 / 298.257_223_563;

/// Rotates a TEME vector into the Earth-fixed frame at `time`.
///
/// This is a rotation by Greenwich mean sidereal time. It ignores polar
/// motion, which moves the result by at most about 15 m.
pub fn teme_to_ecef(teme: [f64; 3], time: JulianDate) -> [f64; 3] {
    let (sin, cos) = time.gmst().sin_cos();
    let [x, y, z] = teme;
    [cos * x + sin * y, -sin * x + cos * y, z]
}

/// A point relative to the WGS-84 ellipsoid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geodetic {
    /// Geodetic latitude, radians, positive north.
    pub latitude: f64,
    /// Longitude, radians in `(-π, π]`, positive east.
    pub longitude: f64,
    /// Height above the ellipsoid, km.
    pub altitude_km: f64,
}

impl Geodetic {
    /// Converts an Earth-fixed position in km.
    pub fn from_ecef(ecef: [f64; 3]) -> Geodetic {
        let [x, y, z] = ecef;
        let e2 = WGS84_F * (2.0 - WGS84_F);
        let p = x.hypot(y);
        let longitude = y.atan2(x);

        // Fixed-point iteration on latitude; converges to 1e-12 rad within a
        // few steps for any point outside the Earth's core.
        let mut latitude = z.atan2(p * (1.0 - e2));
        let mut n = WGS84_A;
        for _ in 0..10 {
            let sin = latitude.sin();
            n = WGS84_A / (1.0 - e2 * sin * sin).sqrt();
            let next = (z + e2 * n * sin).atan2(p);
            let converged = (next - latitude).abs() < 1e-12;
            latitude = next;
            if converged {
                break;
            }
        }
        let (sin, cos) = latitude.sin_cos();
        // Stable at the poles, unlike p / cos(latitude) - n.
        let altitude_km = p * cos + (z + e2 * n * sin) * sin - n;
        Geodetic {
            latitude,
            longitude,
            altitude_km,
        }
    }

    /// Converts a TEME position in km, such as [`State::position`](crate::State::position).
    pub fn from_teme(teme: [f64; 3], time: JulianDate) -> Geodetic {
        Geodetic::from_ecef(teme_to_ecef(teme, time))
    }

    /// The Earth-fixed position of this point, in km.
    pub fn to_ecef(self) -> [f64; 3] {
        let e2 = WGS84_F * (2.0 - WGS84_F);
        let (sin_lat, cos_lat) = self.latitude.sin_cos();
        let (sin_lon, cos_lon) = self.longitude.sin_cos();
        let n = WGS84_A / (1.0 - e2 * sin_lat * sin_lat).sqrt();
        [
            (n + self.altitude_km) * cos_lat * cos_lon,
            (n + self.altitude_km) * cos_lat * sin_lon,
            (n * (1.0 - e2) + self.altitude_km) * sin_lat,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn axis_points() {
        let g = Geodetic::from_ecef([WGS84_A + 400.0, 0.0, 0.0]);
        assert_eq!((g.latitude, g.longitude), (0.0, 0.0));
        assert!((g.altitude_km - 400.0).abs() < 1e-9);

        let polar_radius = WGS84_A * (1.0 - WGS84_F);
        let g = Geodetic::from_ecef([0.0, 0.0, -(polar_radius + 35.0)]);
        assert!((g.latitude + FRAC_PI_2).abs() < 1e-12);
        assert!((g.altitude_km - 35.0).abs() < 1e-9);
    }

    #[test]
    fn round_trips_through_ecef() {
        for (lat, lon, alt) in [
            (51.5, -0.12, 0.05),
            (-33.9, 151.2, 550.0),
            (89.9, 10.0, 35_786.0),
            (0.3, 179.9, 20_200.0),
        ] {
            let point = Geodetic {
                latitude: f64::to_radians(lat),
                longitude: f64::to_radians(lon),
                altitude_km: alt,
            };
            let back = Geodetic::from_ecef(point.to_ecef());
            assert!((back.latitude - point.latitude).abs() < 1e-11, "{lat}");
            assert!((back.longitude - point.longitude).abs() < 1e-12, "{lon}");
            assert!((back.altitude_km - alt).abs() < 1e-6, "{alt}");
        }
    }

    #[test]
    fn teme_to_ecef_rotates_by_gmst() {
        let time = JulianDate::from_calendar(1992, 8, 20, 12, 14, 0.0);
        let gmst = time.gmst();
        // A vector pointing at the Greenwich meridian in TEME lands on +x.
        let teme = [7000.0 * gmst.cos(), 7000.0 * gmst.sin(), 1000.0];
        let ecef = teme_to_ecef(teme, time);
        assert!((ecef[0] - 7000.0).abs() < 1e-9 && ecef[1].abs() < 1e-9);
        assert_eq!(ecef[2], 1000.0);
    }
}
