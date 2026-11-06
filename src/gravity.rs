//! Earth gravity models used by SGP4.

use std::fmt;
use std::str::FromStr;

/// The Earth gravity model SGP4 runs with.
///
/// TLEs are fit with WGS-72, so [`Gravity::Wgs72`] is the default and the
/// right choice for TLE data. The others exist to reproduce published results.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Gravity {
    /// WGS-72 with the low-precision `xke` used by early SGP4 code.
    Wgs72Old,
    #[default]
    Wgs72,
    Wgs84,
}

/// Numeric constants for one [`Gravity`] model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GravityConstants {
    /// Gravitational parameter, km³/s².
    pub mu: f64,
    /// Equatorial radius of the Earth, km.
    pub radius_km: f64,
    /// sqrt(mu) in Earth radii^1.5 per minute.
    pub xke: f64,
    /// Minutes per time unit (1 / `xke`).
    pub tumin: f64,
    pub j2: f64,
    pub j3: f64,
    pub j4: f64,
    pub j3oj2: f64,
}

impl Gravity {
    pub fn constants(self) -> GravityConstants {
        let (mu, radius_km, j2, j3, j4): (f64, f64, f64, f64, f64) = match self {
            Gravity::Wgs72Old | Gravity::Wgs72 => {
                let mu = if self == Gravity::Wgs72Old {
                    398_600.799_64
                } else {
                    398_600.8
                };
                (
                    mu,
                    6378.135,
                    0.001_082_616,
                    -0.000_002_538_81,
                    -0.000_001_655_97,
                )
            }
            Gravity::Wgs84 => (
                398_600.5,
                6378.137,
                0.001_082_629_989_05,
                -0.000_002_532_153_06,
                -0.000_001_610_987_61,
            ),
        };
        let xke = match self {
            Gravity::Wgs72Old => 0.074_366_916_1,
            _ => 60.0 / (radius_km * radius_km * radius_km / mu).sqrt(),
        };
        GravityConstants {
            mu,
            radius_km,
            xke,
            tumin: 1.0 / xke,
            j2,
            j3,
            j4,
            j3oj2: j3 / j2,
        }
    }
}

impl fmt::Display for Gravity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Gravity::Wgs72Old => "wgs72old",
            Gravity::Wgs72 => "wgs72",
            Gravity::Wgs84 => "wgs84",
        })
    }
}

/// Returned when parsing a [`Gravity`] name fails.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnknownGravityModel(pub String);

impl fmt::Display for UnknownGravityModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "unknown gravity model {:?} (expected wgs72old, wgs72, or wgs84)",
            self.0
        )
    }
}

impl std::error::Error for UnknownGravityModel {}

impl FromStr for Gravity {
    type Err = UnknownGravityModel;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "wgs72old" => Ok(Gravity::Wgs72Old),
            "wgs72" => Ok(Gravity::Wgs72),
            "wgs84" => Ok(Gravity::Wgs84),
            _ => Err(UnknownGravityModel(s.to_string())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constants_match_vallado() {
        let old = Gravity::Wgs72Old.constants();
        assert_eq!(old.xke, 0.074_366_916_1);
        assert_eq!(old.mu, 398_600.799_64);

        let wgs72 = Gravity::Wgs72.constants();
        // python-sgp4 2.27 values.
        assert!((wgs72.xke - 0.074_366_916_133_173_42).abs() < 1e-16);
        assert_eq!(wgs72.radius_km, 6378.135);

        let wgs84 = Gravity::Wgs84.constants();
        assert_eq!(wgs84.radius_km, 6378.137);
        assert!((wgs84.xke - 0.074_366_853_168_713_85).abs() < 1e-16);
        assert!((wgs84.tumin * wgs84.xke - 1.0).abs() < 1e-15);
    }

    #[test]
    fn parses_names_case_insensitively() {
        assert_eq!("WGS84".parse::<Gravity>(), Ok(Gravity::Wgs84));
        assert_eq!("wgs72old".parse::<Gravity>(), Ok(Gravity::Wgs72Old));
        assert!("egm96".parse::<Gravity>().is_err());
        for g in [Gravity::Wgs72Old, Gravity::Wgs72, Gravity::Wgs84] {
            assert_eq!(g.to_string().parse::<Gravity>(), Ok(g));
        }
    }
}
