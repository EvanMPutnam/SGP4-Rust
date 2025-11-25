mod alpha5;
mod propagation;
mod functions;
mod earth_gravity;
mod ext;
mod io;


// tests/sgp4_tests.rs  (or put this in a #[cfg(test)] mod in lib.rs)
//
// NOTE: You may need to adjust the `use` paths below to match
// how your crate is organized (module names, SatRec location, etc.)

use std::f64::consts::PI;

use crate::ext::{invjday, newtonnu, rv2coe};
use crate::functions::{days2mdhms, jday, day_of_year_to_month_day};
use crate::io::{compute_checksum, fix_checksum, verify_checksum};
use crate::io::twoline2satrec;
use crate::propagation::SatRec;

// ---------------------------------------------------------------------
// Shared helpers / constants
// ---------------------------------------------------------------------

const LINE1: &str =
    "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
const LINE2: &str =
    "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
const BAD2: &str =
    "2 00007  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413669";

// Handy approx helper (same style as in your existing ext.rs tests)
fn approx(expected: f64, got: f64, eps: f64) {
    let diff = (expected - got).abs();
    assert!(
        diff <= eps,
        "expected {expected}, got {got} (|Δ| = {diff}, eps = {eps})"
    );
}

// ---------------------------------------------------------------------
// 1. Core “SatRec from TLE” / gravity-model tests
// ---------------------------------------------------------------------

#[test]
fn satrec_built_with_twoline2rv_has_basic_identity_fields() {
    // Uses your Rust SatRec::twoline2rv wrapper.
    // Adjust whichconst string if yours differ.
    let sat = SatRec::twoline2rv(LINE1, LINE2, "wgs72");

    // Basic identity from TLE
    assert_eq!(sat.satnum_str, "00005");
    assert_eq!(sat.classification, 'U');
    assert_eq!(sat.operationmode, 'i');

    // Basic orbital parameters – values taken from Vallado / Python sgp4
    // and should match to reasonable precision.
    approx(sat.bstar, 2.8098e-05, 1.0e-10);
    approx(sat.ecco, 0.1859667, 1.0e-10);
    approx(sat.inclo, 34.2682_f64.to_radians(), 1.0e-10);
    approx(sat.nodeo, 348.7242_f64.to_radians(), 1.0e-10);
    approx(sat.argpo, 331.7664_f64.to_radians(), 1.0e-10);
    approx(sat.mo, 19.3264_f64.to_radians(), 1.0e-10);

    // These units should be radians/minute
    approx(sat.no_kozai, 0.04722944544077857, 1.0e-12);
}

#[test]
fn gravity_models_change_position_as_expected() {
    // This mirrors the Python tests:
    // - build the same TLE under three gravity models
    // - propagate to tsince = 309.67110720001529 min
    // - compare position vector

    let tsince = 309.67110720001529;

    // ---- WGS72OLD ----
    let mut sat_old =
        SatRec::twoline2rv(LINE1, LINE2, "wgs72old");
    let (e_old, r_old, _v_old) = sat_old.sgp4_tsince(tsince);
    assert_eq!(e_old, 0);
    approx(r_old[0], -3754.251473242793, 1.0e-6);
    approx(r_old[1], 7876.346815095482, 1.0e-6);
    approx(r_old[2], 4719.220855042922, 1.0e-6);

    // ---- WGS72 ----
    let mut sat_72 =
        SatRec::twoline2rv(LINE1, LINE2, "wgs72");
    let (e_72, r_72, _v_72) = sat_72.sgp4_tsince(tsince);
    assert_eq!(e_72, 0);
    approx(r_72[0], -3754.2514743216166, 1.0e-6);
    approx(r_72[1], 7876.346817439062, 1.0e-6);
    approx(r_72[2], 4719.220856478582, 1.0e-6);

    // ---- WGS84 ----
    let mut sat_84 =
        SatRec::twoline2rv(LINE1, LINE2, "wgs84");
    let (e_84, r_84, _v_84) = sat_84.sgp4_tsince(tsince);
    assert_eq!(e_84, 0);
    approx(r_84[0], -3754.2437675772426, 1.0e-6);
    approx(r_84[1], 7876.3549956188945, 1.0e-6);
    approx(r_84[2], 4719.227897029576, 1.0e-6);
}

// ---------------------------------------------------------------------
// 2. days2mdhms / jday / invjday / calendar helpers
// ---------------------------------------------------------------------

#[test]
fn days2mdhms_example_from_python_suite() {
    // Python test: days2mdhms(2020, 133.35625) == (5, 12, 8, 33, 0.0)
    let (mon, day, hr, min, sec) = days2mdhms(2020, 133.35625, true);
    assert_eq!((mon, day, hr, min), (5, 12, 8, 33));
    approx(sec, 0.0, 1.0e-9);
}

#[test]
fn jday_split_matches_expected() {
    // Python test: jday(2019, 10, 9, 16, 57, 15)
    // returns (2458765.5, 0.7064236111111111)
    let (jd, fr) = jday(2019, 10, 9, 16, 57, 15.0);
    approx(jd, 2458765.5, 1.0e-9);
    approx(fr, 0.706_423_611_111_111_1, 1.0e-12);
}

#[test]
fn jday_and_invjday_round_trip() {
    // Simple round-trip test: (JD -> calendar -> JD)
    let (jd, fr) = jday(2020, 2, 29, 23, 59, 30.5);
    let jd_full = jd + fr;

    let (year, mon, day, hr, min, sec) = invjday(jd_full);

    assert_eq!(year, 2020);
    assert_eq!(mon, 2);
    assert_eq!(day, 29);
    assert_eq!(hr, 23);
    assert_eq!(min, 59);
    approx(sec, 30.5, 2.0e-5); // allow tiny numerical drift

    // And back again:
    let (jd2, fr2) = jday(year, mon, day, hr, min, sec);
    let jd_full2 = jd2 + fr2;
    approx(jd_full, jd_full2, 2.0e-5);
}

#[test]
fn months_and_days_helper_matches_python_logic() {
    // Port of test_months_and_days() from the Python suite.
    // Non-leap year
    let month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    let mut day_of_year = 1;
    for (month, length) in month_lengths.iter().enumerate() {
        for day in 1..=*length {
            let (m, d) = day_of_year_to_month_day(day_of_year, false);
            assert_eq!((m, d), ((month + 1) as i32, day));
            day_of_year += 1;
        }
    }

    // Leap year: February has 29
    let mut month_lengths_leap = month_lengths;
    month_lengths_leap[1] = 29;
    day_of_year = 1;
    for (month, length) in month_lengths_leap.iter().enumerate() {
        for day in 1..=*length {
            let (m, d) = day_of_year_to_month_day(day_of_year, true);
            assert_eq!((m, d), ((month + 1) as i32, day));
            day_of_year += 1;
        }
    }
}

// ---------------------------------------------------------------------
// 3. newtonnu / rv2coe “special case” tests
// ---------------------------------------------------------------------

#[test]
fn newtonnu_elliptic_example_matches_python_values() {
    // Elliptic case from Python tests
    let ecc = 0.1;
    let nu = 40.0_f64.to_radians();

    let (e0, m) = newtonnu(ecc, nu);

    approx(e0, 0.703_368_172_594_7, 1.0e-9);
    approx(m, 0.635_116_435_676_593, 1.0e-9);
}

#[test]
fn newtonnu_parabolic_and_hyperbolic_examples() {
    // Parabolic-ish (ecc = 1.0) and hyperbolic (ecc > 1) tests
    // from the Python file test_hyperbolic_orbit().
    let (e0_p, m_p) = newtonnu(1.0, 2.9); // parabolic branch
    approx(e0_p, 8.238_092_752_965_605, 1.0e-12);
    approx(m_p, 194.600_699_894_828_98, 1.0e-12);

    let (e0_h, m_h) = newtonnu(1.1, 2.7); // hyperbolic branch
    approx(e0_h, 4.262_200_676_156_417, 1.0e-12);
    approx(m_h, 34.761_340_820_283_72, 1.0e-12);
}

#[test]
fn rv2coe_round_trip_simple_orbit() {
    // Simple sanity check of rv2coe:
    // Use a circular-ish LEO and see that we get reasonable elements.
    let mu = 398_600.8_f64; // km^3 / s^2 (WGS-72 mu)
    let r = [7000.0, 0.0, 0.0]; // km
    let v = [0.0, 7.546049108166282, 1.0]; // km/s, slightly inclined

    let (p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper) =
        rv2coe(&r, &v, mu);

    // Semi-major axis and p near radius (circular-ish)
    approx(a, 7000.0, 1.0e2); // big tolerance, just a rough check
    approx(p, a * (1.0 - ecc * ecc), 1.0e2);

    // Eccentricity small, inclination not crazy
    assert!(ecc < 0.1);
    assert!(incl > 0.0 && incl < PI);

    // These should be defined (not “undefined” sentinel)
    assert!(!p.is_nan());
    assert!(!a.is_nan());
    assert!(!incl.is_nan());
    assert!(!omega.is_nan());
    assert!(!argp.is_nan());
    assert!(!nu.is_nan());
    assert!(!m.is_nan());
    assert!(!arglat.is_nan() || arglat == 999_999.1);
    assert!(!truelon.is_nan() || truelon == 999_999.1);
    assert!(!lonper.is_nan() || lonper == 999_999.1);
}

// ---------------------------------------------------------------------
// 4. TLE checksum helpers
// ---------------------------------------------------------------------

#[test]
fn good_tle_checksum_accepts_line_and_round_trips() {
    for line in [LINE1, LINE2] {
        let checksum_char = line.chars().last().unwrap();
        let checksum = checksum_char.to_digit(10).unwrap();

        assert_eq!(compute_checksum(line), checksum);
        assert_eq!(fix_checksum(&line[..68]), line);
        verify_checksum(LINE1, LINE2)
    }
}

// ---------------------------------------------------------------------
// 5. Basic bad-format / mismatched-line tests
// ---------------------------------------------------------------------

#[test]
fn mismatched_object_numbers_between_lines_is_error() {
    let result = twoline2satrec(LINE1, BAD2, "wgs72", 'i');
    assert!(result.is_err());
}

#[test]
fn non_ascii_tle_lines_are_rejected() {
    // Inject a non-ASCII character into line 1
    let mut bad1 = LINE1.replace("23 ", "23 "); // NBSP (U+00A0)
    // Keep line 2 as-is
    let res1 = twoline2satrec(&bad1, LINE2, "wgs72", 'i');
    assert!(res1.is_err());

    // Inject non-ASCII into line 2
    let mut bad2 = LINE2.replace(" 34", " 34"); // NBSP again
    let res2 = twoline2satrec(LINE1, &bad2, "wgs72", 'i');
    assert!(res2.is_err());
}


#[test]
fn usage_example() {
    let s = "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991";
    let t = "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482";
    let mut satrec = SatRec::twoline2rv(s, t, "wgs72");

    let jd = 2458826.5;
    let fr = 0.8625;
    let response_values = satrec.sgp4(jd, fr);
    dbg!(response_values);
}