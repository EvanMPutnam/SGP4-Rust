//! An initialized SGP4 propagator.

use std::f64::consts::TAU;

use crate::elements::Elements;
use crate::error::PropagationError;
use crate::gravity::{Gravity, GravityConstants};
use crate::model::Model;
use crate::time::JulianDate;
use crate::tle::Tle;

/// Which of two historical variants of SGP4 to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum OpsMode {
    /// Reproduces the U.S. Air Force Space Command code, including its
    /// approximate sidereal time and node wrapping.
    Afspc,
    /// Vallado's improved mode. Use this unless you need to match AFSPC output.
    #[default]
    Improved,
}

/// Position and velocity in the TEME (true equator, mean equinox) frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct State {
    /// Kilometers.
    pub position: [f64; 3],
    /// Kilometers per second.
    pub velocity: [f64; 3],
}

/// A satellite ready to propagate.
///
/// Propagation takes `&self`, so a `Satellite` can be shared across threads
/// and propagated to any number of times in any order.
///
/// ```
/// use sgp4_rust::{Satellite, Tle};
///
/// let tle = Tle::parse(
///     "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991",
///     "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482",
/// )?;
/// let iss = Satellite::from_tle(&tle)?;
/// let state = iss.propagate(90.0)?; // 90 minutes after the TLE epoch
/// let [x, y, z] = state.position;
/// assert!(((x * x + y * y + z * z).sqrt() - 6790.0).abs() < 20.0);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct Satellite {
    elements: Elements,
    model: Model,
}

impl Satellite {
    /// Initializes SGP4 for `elements`.
    ///
    /// Fails if the elements are out of range or cannot be propagated even
    /// at their own epoch.
    pub fn new(
        elements: Elements,
        gravity: Gravity,
        mode: OpsMode,
    ) -> Result<Self, PropagationError> {
        if !(0.0..1.0).contains(&elements.eccentricity) {
            return Err(PropagationError::Eccentricity(elements.eccentricity));
        }
        if elements.mean_motion.is_nan() || elements.mean_motion <= 0.0 {
            return Err(PropagationError::MeanMotion(elements.mean_motion));
        }
        let model = Model::new(&elements, gravity.constants(), mode);
        model.propagate(0.0)?;
        Ok(Self { elements, model })
    }

    /// Initializes SGP4 for a TLE with the WGS-72 constants TLEs are fit with.
    pub fn from_tle(tle: &Tle) -> Result<Self, PropagationError> {
        Self::new(tle.elements(), Gravity::Wgs72, OpsMode::Improved)
    }

    /// The state `minutes_since_epoch` minutes after the element epoch.
    /// Negative values propagate backward.
    pub fn propagate(&self, minutes_since_epoch: f64) -> Result<State, PropagationError> {
        self.model.propagate(minutes_since_epoch)
    }

    /// The state at a UTC instant, given as a [`JulianDate`] or a
    /// `chrono::DateTime<Utc>`.
    pub fn propagate_to(&self, time: impl Into<JulianDate>) -> Result<State, PropagationError> {
        self.propagate(time.into().minutes_since(self.elements.epoch))
    }

    pub fn elements(&self) -> &Elements {
        &self.elements
    }

    pub fn epoch(&self) -> JulianDate {
        self.elements.epoch
    }

    pub fn gravity(&self) -> GravityConstants {
        self.model.gravity
    }

    /// Orbital period at epoch, in minutes, from the Brouwer mean motion.
    pub fn period_minutes(&self) -> f64 {
        TAU / self.model.no_unkozai
    }

    /// Mean semi-major axis at epoch, in km.
    pub fn semi_major_axis_km(&self) -> f64 {
        self.model.a * self.model.gravity.radius_km
    }

    /// Mean perigee altitude above the equatorial radius at epoch, in km.
    pub fn perigee_altitude_km(&self) -> f64 {
        self.model.altp * self.model.gravity.radius_km
    }

    /// Mean apogee altitude above the equatorial radius at epoch, in km.
    pub fn apogee_altitude_km(&self) -> f64 {
        self.model.alta * self.model.gravity.radius_km
    }

    /// Whether the period is 225 minutes or more, which switches on the
    /// lunar-solar and resonance terms (the "SDP4" part of SGP4).
    pub fn is_deep_space(&self) -> bool {
        self.period_minutes() >= 225.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{TimeZone, Utc};

    const ISS_1: &str = "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991";
    const ISS_2: &str = "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482";
    const VANGUARD_1: &str =
        "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    const VANGUARD_2: &str =
        "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

    fn iss() -> Satellite {
        Satellite::from_tle(&Tle::parse(ISS_1, ISS_2).unwrap()).unwrap()
    }

    fn assert_close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() <= tolerance, "{actual:?} != {expected:?}");
        }
    }

    #[test]
    fn matches_python_sgp4_for_iss() {
        // python-sgp4 2.27: Satrec.twoline2rv(...).sgp4(2458826.5, 0.8625)
        let state = iss()
            .propagate_to(JulianDate::new(2_458_826.5, 0.8625))
            .unwrap();
        assert_close(
            state.position,
            [-6088.922358183078, -936.1367482888298, -2866.4408776941464],
            1e-8,
        );
        assert_close(
            state.velocity,
            [-1.5255342771668234, -5.5384966257541395, 5.068164791538606],
            1e-11,
        );
    }

    #[test]
    fn datetime_and_julian_date_agree() {
        let sat = iss();
        let time = Utc.with_ymd_and_hms(2019, 12, 9, 20, 42, 0).unwrap();
        let from_datetime = sat.propagate_to(time).unwrap();
        let from_jd = sat
            .propagate_to(JulianDate::new(2_458_826.5, 0.8625))
            .unwrap();
        assert_close(from_datetime.position, from_jd.position, 1e-6);
    }

    #[test]
    fn gravity_models_match_python_sgp4() {
        // python-sgp4 test_all_three_gravity_models_with_sgp4init.
        let elements = Tle::parse(VANGUARD_1, VANGUARD_2).unwrap().elements();
        let t = 309.671_107_200_015_3;
        let cases = [
            (
                Gravity::Wgs72Old,
                [-3754.251473242793, 7876.346815095482, 4719.220855042922],
            ),
            (
                Gravity::Wgs72,
                [-3754.2514743216166, 7876.346817439062, 4719.220856478582],
            ),
            (
                Gravity::Wgs84,
                [-3754.2437675772426, 7876.3549956188945, 4719.227897029576],
            ),
        ];
        for (gravity, expected) in cases {
            let sat = Satellite::new(elements, gravity, OpsMode::Improved).unwrap();
            assert_close(sat.propagate(t).unwrap().position, expected, 1e-6);
        }
    }

    #[test]
    fn rejects_out_of_range_elements() {
        let elements = Tle::parse(ISS_1, ISS_2).unwrap().elements();
        let bad_ecc = Elements {
            eccentricity: 1.2,
            ..elements
        };
        assert_eq!(
            Satellite::new(bad_ecc, Gravity::Wgs72, OpsMode::Improved).unwrap_err(),
            PropagationError::Eccentricity(1.2)
        );
        for n in [0.0, -0.05, f64::NAN] {
            let bad_n = Elements {
                mean_motion: n,
                ..elements
            };
            let err = Satellite::new(bad_n, Gravity::Wgs72, OpsMode::Improved).unwrap_err();
            assert_eq!(err.code(), 2);
        }
    }

    #[test]
    fn orbit_summary_for_iss() {
        let sat = iss();
        assert!((sat.period_minutes() - 92.9).abs() < 0.1);
        assert!((sat.perigee_altitude_km() - 412.0).abs() < 10.0);
        assert!((sat.apogee_altitude_km() - 422.0).abs() < 10.0);
        assert!(!sat.is_deep_space());
    }

    #[test]
    fn satellites_are_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<Satellite>();
    }
}
