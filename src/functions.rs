//! General-purpose date/time routines used by SGP4.
//!
//! Ported from Brandon Rhodes' python-sgp4 `sgp4.ext` helpers.

/// Compute the Julian date as (integer-ish part, fractional part).
///
/// This mirrors Vallado’s `jday()` and the python-sgp4 behavior:
/// the first value is the JD at midnight starting the calendar day
/// (always `.5` fractional part), and the second is the fractional
/// offset for the given time of day.
///
/// Example from python-sgp4:
/// ```text
/// jday(2020, 2, 11, 13, 57, 0) -> (2458890.5, 0.58125)
/// ```
pub fn jday(
    year: i32,
    mon: i32,
    day: i32,
    hr: i32,
    minute: i32,
    sec: f64,
) -> (f64, f64) {
    // Careful port of the weird floor/float math from Python.
    let year_f = year as f64;
    let mon_f = mon as f64;
    let day_f = day as f64;

    let term1 = 367.0 * year_f;

    let term2 = {
        // ((mon + 9) // 12.0)
        let m = ((mon as f64 + 9.0) / 12.0).floor();
        // 7 * (year + m) * 0.25 // 1.0
        (7.0 * (year_f + m) * 0.25).floor()
    };

    let term3 = {
        // 275 * mon / 9.0 // 1.0
        (275.0 * mon_f / 9.0).floor()
    };

    let jd = term1 - term2 + term3 + day_f + 1721013.5;

    let fr = (sec + (minute as f64) * 60.0 + (hr as f64) * 3600.0) / 86400.0;

    (jd, fr)
}

/// Convert a fractional "day of year" into (month, day, hour, minute, seconds).
///
/// `days` is 1.0 for the start of Jan 1, 2.0 for the start of Jan 2, etc.
/// Seconds are optionally rounded to the nearest microsecond if
/// `round_to_microsecond` is true (mirroring the original intent).
pub fn days2mdhms(
    year: i32,
    days: f64,
    round_to_microsecond: bool,
) -> (i32, i32, i32, i32, f64) {
    // seconds since start of the year
    let mut second = days * 86400.0;

    if round_to_microsecond {
        second = round_to_n_decimals(second, 6);
    }

    // minute, second = divmod(second, 60.0)
    let minute_f = (second / 60.0).floor();
    second -= minute_f * 60.0;

    if round_to_microsecond {
        second = round_to_n_decimals(second, 6);
    }

    let mut minute = minute_f as i32;

    // hour, minute = divmod(minute, 60)
    let mut hour = minute / 60;
    minute %= 60;

    // day_of_year, hour = divmod(hour, 24)
    let mut day_of_year = hour / 24;
    hour %= 24;

    // Python day_of_year is 0-based after that divmod, but the helper
    // assumes a 1-based count inside its formula, so this matches.
    let is_leap = year % 400 == 0 || (year % 4 == 0 && year % 100 != 0);

    let (mut month, mut day) = day_of_year_to_month_day(day_of_year, is_leap);

    // Behave like original in case of overflow
    if month == 13 {
        month = 12;
        day += 31;
    }

    (month, day, hour, minute, second)
}

/// Internal helper: turn "day of year" into (month, day).
///
/// Port of `_day_of_year_to_month_day()` from python-sgp4.
pub fn day_of_year_to_month_day(day_of_year: i32, is_leap: bool) -> (i32, i32) {
    // In the Python code, is_leap and comparisons are used as ints (0/1).
    let leap_i = if is_leap { 1 } else { 0 };

    let february_bump = {
        let cond = if day_of_year >= 60 + leap_i { 1 } else { 0 };
        (2 - leap_i) * cond
    };

    let august = if day_of_year >= 215 { 1 } else { 0 };

    // month, day = divmod(2 * (day_of_year - 1 + 30 * august + february_bump), 61)
    let num = 2 * (day_of_year - 1 + 30 * august + february_bump);
    let mut month = num / 61;
    let mut day = num % 61;

    month += 1 - august;
    day /= 2;
    day += 1;

    (month, day)
}

/// Round a f64 to `n` decimal places.
fn round_to_n_decimals(x: f64, n: i32) -> f64 {
    if n <= 0 {
        return x.round();
    }
    let factor = 10_f64.powi(n);
    (x * factor).round() / factor
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) {
        let diff = (a - b).abs();
        assert!(
            diff <= tol,
            "expected {b:.15e}, got {a:.15e}, diff={diff:.3e}"
        );
    }

    #[test]
    fn jday_matches_python_example() {
        // From the python docstring:
        // >>> jd, fr = jday(2020, 2, 11, 13, 57, 0)
        // jd == 2458890.5, fr == 0.58125
        let (jd, fr) = jday(2020, 2, 11, 13, 57, 0.0);

        approx_eq(jd, 2458890.5, 1.0e-9);
        approx_eq(fr, 0.58125, 1.0e-12);
    }

    #[test]
    fn days2mdhms_basic_cases() {
        // days2mdhms(2000, 1.0) -> (1, 1, 0, 0, 0.0)
        let (m1, d1, h1, min1, s1) = days2mdhms(2000, 1.0, true);
        assert_eq!((m1, d1, h1, min1), (1, 1, 0, 0));
        approx_eq(s1, 0.0, 1.0e-9);

        // days2mdhms(2000, 32.0) -> (2, 1, 0, 0, 0.0)  # Feb 1, leap year
        let (m2, d2, h2, min2, s2) = days2mdhms(2000, 32.0, true);
        assert_eq!((m2, d2, h2, min2), (2, 1, 0, 0));
        approx_eq(s2, 0.0, 1.0e-9);

        // days2mdhms(2000, 366.0) -> (12, 31, 0, 0, 0.0)
        let (m3, d3, h3, min3, s3) = days2mdhms(2000, 366.0, true);
        assert_eq!((m3, d3, h3, min3), (12, 31, 0, 0));
        approx_eq(s3, 0.0, 1.0e-9);
    }

    #[test]
    fn day_of_year_to_month_day_leap_and_nonleap() {
        // Non-leap: 2001-03-01 is day 60 (1-based); after the divmods it's 59 here
        let (m, d) = day_of_year_to_month_day(59, false);
        // That's Feb 28 in a non-leap year
        assert_eq!((m, d), (2, 28));

        let (m2, d2) = day_of_year_to_month_day(60, false);
        // Mar 1
        assert_eq!((m2, d2), (3, 1));

        // Leap year behavior: for a leap year, day_of_year=60 should be Feb 29.
        let (m3, d3) = day_of_year_to_month_day(60, true);
        assert_eq!((m3, d3), (2, 29));
    }
}
