use chrono::{DateTime, TimeZone, Utc};
use std::f64::consts::PI;

use crate::ext::{invjday, jday};
use crate::functions::days2mdhms;
use crate::propagation::{sgp4init, SatRec};

const LINE1_FMT: &str =
    "1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN";
const LINE2_FMT: &str =
    "2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN";

fn tle_error_message(line_no: u8, given: &str) -> String {
    format!(
        "TLE format error\n\n...\n\n{0}\n{1}",
        if line_no == 1 { LINE1_FMT } else { LINE2_FMT },
        given
    )
}

pub fn twoline2satrec(
    longstr1: &str,
    longstr2: &str,
    whichconst: &str,   // <-- now &str, matching sgp4init
    opsmode: char,
) -> Result<SatRec, String> {
    let deg2rad = PI / 180.0;
    let xpdotp = 1440.0 / (2.0 * PI); // 229.1831180523293

    let line1 = longstr1.trim_end_matches(&['\r', '\n'][..]).to_string();
    let line2 = longstr2.trim_end_matches(&['\r', '\n'][..]).to_string();

    if !line1.is_ascii() || !line2.is_ascii() {
        return Err(format!(
            "your TLE lines are broken because they contain non-ASCII characters:\n\n{}\n{}",
            line1, line2
        ));
    }

    // --- parse line 1 (same as before) ---
    if !(line1.len() >= 64
        && line1.starts_with("1 ")
        && line1.as_bytes().get(8) == Some(&b' ')
        && line1.as_bytes().get(23) == Some(&b'.')
        && line1.as_bytes().get(32) == Some(&b' ')
        && line1.as_bytes().get(34) == Some(&b'.')
        && line1.as_bytes().get(43) == Some(&b' ')
        && line1.as_bytes().get(52) == Some(&b' ')
        && line1.as_bytes().get(61) == Some(&b' ')
        && line1.as_bytes().get(63) == Some(&b' '))
    {
        return Err(tle_error_message(1, &line1));
    }

    let satnum_str = line1[2..7].trim().to_string();
    let two_digit_year: i32 = line1[18..20].trim().parse()
        .map_err(|e| format!("failed to parse epoch year from line 1: {e}"))?;
    let epochdays: f64 = line1[20..32].trim().parse()
        .map_err(|e| format!("failed to parse epochdays from line 1: {e}"))?;
    let ndot_raw: f64 = line1[33..43].trim().parse()
        .map_err(|e| format!("failed to parse ndot from line 1: {e}"))?;
    let nddot_str = format!("{}.{}", &line1[44..45], &line1[45..50]);
    let mut nddot: f64 = nddot_str.trim().parse()
        .map_err(|e| format!("failed to parse nddot from line 1: {e}"))?;
    let nexp: i32 = line1[50..52].trim().parse()
        .map_err(|e| format!("failed to parse nddot exponent from line 1: {e}"))?;
    let bstar_str = format!("{}.{}", &line1[53..54], &line1[54..59]);
    let mut bstar: f64 = bstar_str.trim().parse()
        .map_err(|e| format!("failed to parse bstar from line 1: {e}"))?;
    let ibexp: i32 = line1[59..61].trim().parse()
        .map_err(|e| format!("failed to parse bstar exponent from line 1: {e}"))?;
    let classification = line1.chars().nth(7).unwrap_or('U');
    let intldesg = line1[9..17].trim_end().to_string();
    let ephtype = line1.chars().nth(62).unwrap_or(' ');
    let elnum: i32 = line1[64..68].trim().parse()
        .map_err(|e| format!("failed to parse elnum from line 1: {e}"))?;

    // --- parse line 2 (same as before) ---
    if !(line2.len() >= 68
        && line2.starts_with("2 ")
        && line2.as_bytes().get(7) == Some(&b' ')
        && line2.as_bytes().get(11) == Some(&b'.')
        && line2.as_bytes().get(16) == Some(&b' ')
        && line2.as_bytes().get(20) == Some(&b'.')
        && line2.as_bytes().get(25) == Some(&b' ')
        && line2.as_bytes().get(33) == Some(&b' ')
        && line2.as_bytes().get(37) == Some(&b'.')
        && line2.as_bytes().get(42) == Some(&b' ')
        && line2.as_bytes().get(46) == Some(&b'.')
        && line2.as_bytes().get(51) == Some(&b' '))
    {
        return Err(tle_error_message(2, &line2));
    }

    if satnum_str != line2[2..7].trim() {
        return Err("Object numbers in lines 1 and 2 do not match".to_string());
    }

    let inclo_deg: f64 = line2[8..16].trim().parse()
        .map_err(|e| format!("failed to parse inclination from line 2: {e}"))?;
    let nodeo_deg: f64 = line2[17..25].trim().parse()
        .map_err(|e| format!("failed to parse RAAN from line 2: {e}"))?;
    let ecco_str = format!("0.{}", line2[26..33].replace(' ', "0"));
    let ecco: f64 = ecco_str.parse()
        .map_err(|e| format!("failed to parse eccentricity from line 2: {e}"))?;
    let argpo_deg: f64 = line2[34..42].trim().parse()
        .map_err(|e| format!("failed to parse argpo from line 2: {e}"))?;
    let mo_deg: f64 = line2[43..51].trim().parse()
        .map_err(|e| format!("failed to parse mo from line 2: {e}"))?;
    let no_kozai_rev_per_day: f64 = line2[52..63].trim().parse()
        .map_err(|e| format!("failed to parse mean motion from line 2: {e}"))?;
    let revnum = line2[63..68].trim().to_string();

    // --- unit conversions (same as Python) ---

    // mean motion rev/day -> rad/min
    let no_kozai = no_kozai_rev_per_day / xpdotp;

    // apply exponents to nddot, bstar
    nddot *= 10f64.powi(nexp);
    bstar *= 10f64.powi(ibexp);

    // convert to sgp4 units
    let ndot_sgp4 = ndot_raw / (xpdotp * 1440.0);
    let nddot_sgp4 = nddot / (xpdotp * 1440.0 * 1440.0);

    // angles to radians
    let inclo = inclo_deg * deg2rad;
    let nodeo = nodeo_deg * deg2rad;
    let argpo = argpo_deg * deg2rad;
    let mo = mo_deg * deg2rad;

    // --- epoch handling (Vallado 57–2056 rule) ---
    let year = if two_digit_year < 57 {
        two_digit_year + 2000
    } else {
        two_digit_year + 1900
    };

    let (mon, day, hr, minute, sec) = days2mdhms(year, epochdays, true);

    let jdsatepoch = jday(
        year,
        mon as i32,
        day as i32,
        hr as i32,
        minute as i32,
        sec,
    );

    // ----------------- build SatRec -----------------
    let mut satrec = SatRec {
        error: 0,
        error_message: None,
        operationmode: opsmode,
        satnum_str: satnum_str.clone(),
        classification,
        bstar,
        ndot: ndot_sgp4,
        nddot: nddot_sgp4,
        ecco,
        argpo,
        inclo,
        mo,
        no_kozai,
        nodeo,
        // everything else default for now; sgp4init will fill most fields
        ..Default::default()
    };

    // ------------- call sgp4init with correct args -------------
    //
    // epoch here is in *days since 1950-01-01 00:00*:
    //   epoch = jdsatepoch - 2433281.5
    //
    let epoch_days_1950 = jdsatepoch - 2433281.5;

    let _ok = sgp4init(
        whichconst,          // &str, e.g. "wgs72", "wgs84"
        opsmode,             // 'i' or 'a'
        &satnum_str,         // satn: &str
        epoch_days_1950,     // epoch
        bstar,               // xbstar
        ndot_sgp4,           // xndot
        nddot_sgp4,          // xnddot
        ecco,                // xecco
        argpo,               // xargpo
        inclo,               // xinclo
        mo,                  // xmo
        no_kozai,            // xno_kozai
        nodeo,               // xnodeo
        &mut satrec,         // satrec
    );

    // If you want to enforce sgp4init success:
    // if !_ok {
    //     return Err(satrec.error_message.unwrap_or_else(|| "sgp4init failed".to_string()));
    // }

    Ok(satrec)
}

/// Compute the TLE checksum for the given line (first 68 columns).
pub fn compute_checksum(line: &str) -> u32 {
    line.chars()
        .take(68)
        .map(|c| {
            if c.is_ascii_digit() {
                (c as u32 - b'0' as u32)
            } else if c == '-' {
                1
            } else {
                0
            }
        })
        .sum::<u32>()
        % 10
}

/// Return a new copy of the TLE `line`, with the correct checksum appended.
///
/// This discards any existing checksum at the end of the line, if a checksum
/// is already present.
pub fn fix_checksum(line: &str) -> String {
    // Keep only the first 68 chars and pad if needed.
    let body: String = line.chars().take(68).collect();
    let body_padded = format!("{:<68}", body);
    let checksum = compute_checksum(&body_padded);
    format!("{body_padded}{checksum}")
}

/// Verify the checksum of *one or two* TLE lines.
///
/// Panics if any of the lines fails its checksum, matching the Python
/// behavior (which raises ValueError). This is convenient for use with
/// `#[should_panic]` tests.
pub fn verify_checksum(line1: &str, line2: &str) {
    for line in [line1, line2] {
        // Column 69 (index 68) is the checksum character if present
        let checksum_char = line.chars().nth(68);

        // If there is no checksum digit, skip (same behavior as Python: `if not checksum.isdigit(): continue`)
        let Some(c) = checksum_char else {
            continue;
        };

        if !c.is_ascii_digit() {
            continue;
        }

        let given = c as u32 - b'0' as u32;
        let computed = compute_checksum(line);

        if given != computed {
            // Match the spirit of the Python error message
            panic!(
                "TLE line gives its checksum as {} but in fact tallies to {}:\n{}",
                given, computed, line
            );
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Simple helper for float comparisons.
    fn approx(expected: f64, got: f64, tol: f64) {
        let diff = (expected - got).abs();
        assert!(
            diff <= tol,
            "expected {}, got {} (|Δ| = {})",
            expected,
            got,
            diff
        );
    }

    #[test]
    fn twoline2satrec_parses_iss_tle_basic_fields() {
        // ISS (ZARYA) example from python-sgp4 docs (WGS72 default).
        // See: https://pypi.org/project/sgp4/ (Epoch section)
        let line1 = "1 25544U 98067A   19343.69339541  .00001264  00000-0  29621-4 0  9990";
        let line2 = "2 25544  51.6436 353.2489 0007413  85.7752  45.1475 15.50112970199556";

        // Adjust function name if yours is `twoline2rv` instead.
        let mut satrec = twoline2satrec(line1, line2, "wgs72", 'i')
            .expect("failed to parse ISS TLE");

        // ID / bookkeeping bits
        assert_eq!(satrec.satnum_str, "25544");
        assert_eq!(satrec.classification, 'U');

        // Orbital elements in *radians* / dimensionless, matching TLE values.
        approx(
            satrec.inclo,
            51.6436_f64.to_radians(),
            1.0e-12,
        );
        approx(
            satrec.nodeo,
            353.2489_f64.to_radians(),
            1.0e-12,
        );
        approx(satrec.ecco, 0.0007413, 1.0e-10);
        approx(
            satrec.argpo,
            85.7752_f64.to_radians(),
            1.0e-12,
        );
        approx(
            satrec.mo,
            45.1475_f64.to_radians(),
            1.0e-12,
        );

        // Mean motion: revs/day → rad/min like Vallado / python-sgp4.
        let mean_motion_revs_per_day = 15.50112970199556_f64;
        let xpdotp = 1440.0 / (2.0 * PI); // revs/day to rad/min factor
        let expected_no_kozai = mean_motion_revs_per_day / xpdotp;
        approx(satrec.no_kozai, expected_no_kozai, 1.0e-10);

        // sgp4init should have completed without errors.
        assert_eq!(satrec.error, 0, "sgp4init set error = {}", satrec.error);
    }

    #[test]
    fn verify_checksum_accepts_valid_tle() {
        // ISS TLE with *corrected* checksum on line1
        let line1 = "1 25544U 98067A   19343.69339541  .00001264  00000-0  29621-4 0  9997";
        let line2 = "2 25544  51.6436 353.2489 0007413  85.7752  45.1475 15.50112970199556";

        // If verify_checksum panics on error:
        verify_checksum(line1, line2);

    }

    #[test]
    #[should_panic(expected = "TLE line gives its checksum as")]
    fn verify_checksum_rejects_corrupted_tle() {
        // Deliberately corrupt one digit in line 1 so the checksum fails.
        let bad_line1 =
            "1 25544U 98067A   19343.69339541  .90001264  00000-0  29621-4 0  9990";
        let line2 =
            "2 25544  51.6436 353.2489 0007413  85.7752  45.1475 15.50112970199556";

        // Expect verify_checksum() to panic with the standard message
        verify_checksum(bad_line1, line2);
    }

    #[test]
    fn fix_checksum_round_trips_to_valid_line() {
        // Line with a bogus checksum at the end — last digit is intentionally wrong.
        let bad_line1 =
            "1 25544U 98067A   19343.69339541  .00001264  00000-0  29621-4 0  9999";
        let line2 =
            "2 25544  51.6436 353.2489 0007413  85.7752  45.1475 15.50112970199556";

        // Compute a corrected line and make sure verify_checksum accepts it.
        let fixed_line1 = fix_checksum(bad_line1);
        verify_checksum(&fixed_line1, &line2);
    }
}
