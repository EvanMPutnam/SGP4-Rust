//! Julian dates and sidereal time.

use chrono::{DateTime, TimeDelta, Utc};
use std::f64::consts::TAU;

const MINUTES_PER_DAY: f64 = 1440.0;
const SECONDS_PER_DAY: f64 = 86_400.0;
/// Julian date of 1970-01-01 00:00 UTC.
const UNIX_EPOCH_JD: f64 = 2_440_587.5;
/// Julian date of J2000.0 (2000-01-01 12:00 TT).
const J2000_JD: f64 = 2_451_545.0;

/// A UTC instant as a Julian date held in two parts.
///
/// A single `f64` near 2.46 million only resolves about 40 µs, so the date is
/// kept as a `day` (usually the midnight `.5` value) plus a small `fraction`.
/// This is the same split Vallado uses for `jdsatepoch` / `jdsatepochF`.
/// The parts are not normalized: `(2458826.5, 1.25)` is a valid date.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JulianDate {
    pub day: f64,
    pub fraction: f64,
}

impl JulianDate {
    pub const fn new(day: f64, fraction: f64) -> Self {
        Self { day, fraction }
    }

    /// Builds a date from UTC calendar fields with Vallado's `jday` formula.
    ///
    /// `day` holds the Julian date of the calendar day's midnight and
    /// `fraction` the time of day. Valid for years 1900 through 2100.
    pub fn from_calendar(
        year: i32,
        month: u32,
        day: u32,
        hour: u32,
        minute: u32,
        second: f64,
    ) -> Self {
        let (year, month) = (f64::from(year), f64::from(month));
        let jd = 367.0 * year - (7.0 * (year + ((month + 9.0) / 12.0).floor()) * 0.25).floor()
            + (275.0 * month / 9.0).floor()
            + f64::from(day)
            + 1_721_013.5;
        let fraction =
            (second + f64::from(minute) * 60.0 + f64::from(hour) * 3600.0) / SECONDS_PER_DAY;
        Self::new(jd, fraction)
    }

    /// The date as a single `f64`, at reduced precision.
    pub fn value(self) -> f64 {
        self.day + self.fraction
    }

    /// Minutes from `earlier` to `self`, computed without losing precision.
    pub fn minutes_since(self, earlier: JulianDate) -> f64 {
        ((self.day - earlier.day) + (self.fraction - earlier.fraction)) * MINUTES_PER_DAY
    }

    /// The date `minutes` later (or earlier, if negative).
    pub fn add_minutes(self, minutes: f64) -> Self {
        Self::new(self.day, self.fraction + minutes / MINUTES_PER_DAY)
    }

    /// Greenwich mean sidereal time in radians, in `[0, 2π)`.
    ///
    /// Uses the IAU-82 expression from Vallado's `gstime`, treating UTC as UT1.
    pub fn gmst(self) -> f64 {
        let tut1 = ((self.day - J2000_JD) + self.fraction) / 36_525.0;
        let seconds = -6.2e-6 * tut1 * tut1 * tut1
            + 0.093_104 * tut1 * tut1
            + (876_600.0 * 3600.0 + 8_640_184.812_866) * tut1
            + 67_310.548_41;
        // 1 second of sidereal time is 1/240 degree.
        (seconds.to_radians() / 240.0).rem_euclid(TAU)
    }

    /// Converts to a chrono UTC timestamp, rounded to the nanosecond.
    ///
    /// Returns `None` if the date is not finite or is outside chrono's range.
    pub fn to_datetime(self) -> Option<DateTime<Utc>> {
        let offset = self.day - UNIX_EPOCH_JD;
        let whole_days = offset.floor();
        let nanos = ((offset - whole_days) + self.fraction) * SECONDS_PER_DAY * 1e9;
        if !nanos.is_finite() || whole_days.abs() > 1e8 {
            return None;
        }
        DateTime::UNIX_EPOCH
            .checked_add_signed(TimeDelta::try_days(whole_days as i64)?)?
            .checked_add_signed(TimeDelta::nanoseconds(nanos.round() as i64))
    }
}

impl From<DateTime<Utc>> for JulianDate {
    fn from(time: DateTime<Utc>) -> Self {
        let seconds = time.timestamp();
        let day = seconds.div_euclid(86_400) as f64 + UNIX_EPOCH_JD;
        let second_of_day =
            seconds.rem_euclid(86_400) as f64 + f64::from(time.timestamp_subsec_nanos()) * 1e-9;
        Self::new(day, second_of_day / SECONDS_PER_DAY)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    #[test]
    fn from_calendar_matches_python_sgp4() {
        // python-sgp4: jday(2019, 10, 9, 16, 57, 15) == (2458765.5, 0.7064236111111111)
        let jd = JulianDate::from_calendar(2019, 10, 9, 16, 57, 15.0);
        assert_eq!(jd.day, 2_458_765.5);
        assert!((jd.fraction - 0.706_423_611_111_111_1).abs() < 1e-15);

        let jd = JulianDate::from_calendar(2020, 2, 11, 13, 57, 0.0);
        assert_eq!(jd.day, 2_458_890.5);
        assert!((jd.fraction - 0.58125).abs() < 1e-15);
    }

    #[test]
    fn datetime_round_trip_keeps_microseconds() {
        let time = Utc.with_ymd_and_hms(2020, 2, 29, 23, 59, 30).unwrap()
            + TimeDelta::microseconds(500_123);
        let jd = JulianDate::from(time);
        assert_eq!(
            jd,
            JulianDate::from_calendar(2020, 2, 29, 23, 59, 30.500_123)
        );
        let back = jd.to_datetime().unwrap();
        assert!((back - time).num_nanoseconds().unwrap().abs() <= 1);
    }

    #[test]
    fn datetime_before_unix_epoch() {
        let time = Utc.with_ymd_and_hms(1958, 3, 17, 6, 0, 0).unwrap();
        let jd = JulianDate::from(time);
        assert_eq!(jd, JulianDate::from_calendar(1958, 3, 17, 6, 0, 0.0));
        assert_eq!(jd.to_datetime(), Some(time));
    }

    #[test]
    fn minutes_since_is_exact_across_days() {
        let a = JulianDate::new(2_458_826.5, 0.69339541);
        let b = JulianDate::new(2_458_827.5, 0.19339541);
        assert!((b.minutes_since(a) - 720.0).abs() < 1e-9);
        assert!((a.add_minutes(720.0).minutes_since(b)).abs() < 1e-9);
    }

    #[test]
    fn gmst_matches_vallado_example() {
        // Vallado, Fundamentals of Astrodynamics, example 3-5:
        // 1992-08-20 12:14 UT1 gives GMST = 152.578 787 886 degrees.
        let jd = JulianDate::from_calendar(1992, 8, 20, 12, 14, 0.0);
        assert!((jd.gmst().to_degrees() - 152.578_787_886).abs() < 1e-6);
    }

    #[test]
    fn non_finite_dates_have_no_datetime() {
        assert_eq!(JulianDate::new(f64::NAN, 0.0).to_datetime(), None);
    }
}
