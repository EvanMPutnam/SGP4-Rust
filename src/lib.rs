//! SGP4/SDP4 satellite propagation.
//!
//! A Rust port of David Vallado's reference SGP4 code (2020 revision), by way
//! of Brandon Rhodes's python-sgp4. It reproduces Vallado's published test
//! vectors to within 2e-7 km.
//!
//! ```
//! use sgp4_rust::{Geodetic, JulianDate, Satellite, Tle};
//!
//! let tle = Tle::parse(
//!     "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991",
//!     "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482",
//! )?;
//! let iss = Satellite::from_tle(&tle)?;
//!
//! // Position and velocity in the TEME frame, km and km/s.
//! let time = JulianDate::from_calendar(2019, 12, 9, 20, 42, 0.0);
//! let state = iss.propagate_to(time)?;
//!
//! // Where that is over the Earth.
//! let ground = Geodetic::from_teme(state.position, time);
//! println!(
//!     "{:.2}°, {:.2}° at {:.0} km",
//!     ground.latitude.to_degrees(),
//!     ground.longitude.to_degrees(),
//!     ground.altitude_km
//! );
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! [`Satellite::propagate_to`] also accepts a `chrono::DateTime<Utc>`.
//! Read whole TLE files with [`tle::parse_all`].

pub mod alpha5;
mod elements;
mod error;
pub mod frames;
mod gravity;
mod model;
pub mod orbit;
mod satellite;
mod time;
pub mod tle;

pub use elements::Elements;
pub use error::{PropagationError, TleError, TleLine};
pub use frames::Geodetic;
pub use gravity::{Gravity, GravityConstants, UnknownGravityModel};
pub use satellite::{OpsMode, Satellite, State};
pub use time::JulianDate;
pub use tle::Tle;
