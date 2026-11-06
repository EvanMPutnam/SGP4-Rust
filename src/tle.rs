//! Two-line element sets.

use std::f64::consts::TAU;
use std::str::FromStr;

use crate::alpha5;
use crate::elements::Elements;
use crate::error::{TleError, TleLine};
use crate::time::JulianDate;

/// Revolutions per day per radian per minute (Vallado's `xpdotp`).
const XPDOTP: f64 = 1440.0 / TAU;

/// One parsed TLE, with values in the units the format prints them in.
///
/// Convert to SGP4 units with [`Tle::elements`].
#[derive(Debug, Clone, PartialEq)]
pub struct Tle {
    /// The title line of a three-line TLE, if there was one.
    pub name: Option<String>,
    /// NORAD catalog number, decoded from Alpha-5 if needed.
    pub catalog_number: u32,
    /// `U` (unclassified), `C` (classified), or `S` (secret).
    pub classification: char,
    /// COSPAR ID such as `98067A`. Empty for some old element sets.
    pub international_designator: String,
    pub epoch: JulianDate,
    /// First derivative of mean motion divided by two, rev/day².
    pub mean_motion_dot: f64,
    /// Second derivative of mean motion divided by six, rev/day³.
    pub mean_motion_ddot: f64,
    /// Drag term, in inverse Earth radii.
    pub bstar: f64,
    pub ephemeris_type: u8,
    pub element_set_number: u32,
    pub inclination_deg: f64,
    /// Right ascension of the ascending node, degrees.
    pub raan_deg: f64,
    pub eccentricity: f64,
    pub arg_of_perigee_deg: f64,
    pub mean_anomaly_deg: f64,
    /// Mean motion, revolutions per day.
    pub mean_motion: f64,
    pub revolution_number: u32,
}

impl Tle {
    /// Parses the two data lines of a TLE.
    ///
    /// Trailing whitespace is ignored, as is anything past column 69 (some
    /// files append extra data there). Checksums are not verified; call
    /// [`verify_checksum`] for that.
    pub fn parse(line1: &str, line2: &str) -> Result<Tle, TleError> {
        let line1 = Columns::new(line1, TleLine::One)?;
        let line2 = Columns::new(line2, TleLine::Two)?;

        line1.expect_layout(
            b'1',
            64,
            &[
                (8, b' '),
                (23, b'.'),
                (32, b' '),
                (34, b'.'),
                (43, b' '),
                (52, b' '),
                (61, b' '),
                (63, b' '),
            ],
        )?;
        line2.expect_layout(
            b'2',
            63,
            &[
                (7, b' '),
                (11, b'.'),
                (16, b' '),
                (20, b'.'),
                (25, b' '),
                (33, b' '),
                (37, b'.'),
                (42, b' '),
                (46, b'.'),
                (51, b' '),
            ],
        )?;

        let catalog_number = line1.catalog_number()?;
        if line2.catalog_number()? != catalog_number {
            return Err(TleError::CatalogMismatch);
        }

        let two_digit_year: i32 = line1.parse(18..20, "epoch year")?;
        let epoch_year = if two_digit_year < 57 {
            2000 + two_digit_year
        } else {
            1900 + two_digit_year
        };
        let epoch_day: f64 = line1.parse(20..32, "epoch day")?;

        Ok(Tle {
            name: None,
            catalog_number,
            classification: line1.text[7] as char,
            international_designator: line1.str(9..17).trim().to_string(),
            epoch: tle_epoch(epoch_year, epoch_day),
            mean_motion_dot: line1.parse(33..43, "mean motion derivative")?,
            mean_motion_ddot: line1.exponent_field(44, "mean motion second derivative")?,
            bstar: line1.exponent_field(53, "B*")?,
            ephemeris_type: line1.parse_or_zero(62..63, "ephemeris type")?,
            element_set_number: line1.parse_or_zero(64..68, "element set number")?,
            inclination_deg: line2.parse(8..16, "inclination")?,
            raan_deg: line2.parse(17..25, "right ascension")?,
            eccentricity: format!("0.{}", line2.str(26..33).replace(' ', "0"))
                .parse()
                .map_err(|_| line2.field_error(26..33, "eccentricity"))?,
            arg_of_perigee_deg: line2.parse(34..42, "argument of perigee")?,
            mean_anomaly_deg: line2.parse(43..51, "mean anomaly")?,
            mean_motion: line2.parse(52..63, "mean motion")?,
            revolution_number: line2.parse_or_zero(63..68, "revolution number")?,
        })
    }

    /// The elements in SGP4 units, ready for [`Satellite::new`](crate::Satellite::new).
    pub fn elements(&self) -> Elements {
        Elements {
            epoch: self.epoch,
            bstar: self.bstar,
            // Vallado's unit conversions, which leave out the factors of 2 and 6
            // folded into the printed values. SGP4 does not use either term.
            mean_motion_dot: self.mean_motion_dot / (XPDOTP * 1440.0),
            mean_motion_ddot: self.mean_motion_ddot / (XPDOTP * 1440.0 * 1440.0),
            eccentricity: self.eccentricity,
            arg_of_perigee: self.arg_of_perigee_deg.to_radians(),
            inclination: self.inclination_deg.to_radians(),
            mean_anomaly: self.mean_anomaly_deg.to_radians(),
            mean_motion: self.mean_motion / XPDOTP,
            raan: self.raan_deg.to_radians(),
        }
    }
}

/// Parses a single TLE, with or without a title line.
impl FromStr for Tle {
    type Err = TleError;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let mut all = parse_all(text);
        let tle = all.next().ok_or(TleError::Unpaired(TleLine::One))??;
        match all.next() {
            None => Ok(tle),
            Some(_) => Err(TleError::Layout(TleLine::One)),
        }
    }
}

/// Parses every TLE in a file's contents.
///
/// Accepts two-line and three-line (titled) sets, mixed freely. A title line
/// may start with `0 `, as in CelesTrak's 3LE format. Blank lines are skipped.
/// Each set yields its own `Result`, so one bad entry does not stop the rest.
pub fn parse_all(text: &str) -> impl Iterator<Item = Result<Tle, TleError>> + '_ {
    let mut lines = text
        .lines()
        .map(str::trim_end)
        .filter(|line| !line.is_empty())
        .peekable();
    std::iter::from_fn(move || {
        let mut name = None;
        loop {
            let line = lines.next()?;
            if line.starts_with("1 ") {
                let Some(line2) = lines.next_if(|next| next.starts_with("2 ")) else {
                    return Some(Err(TleError::Unpaired(TleLine::One)));
                };
                return Some(Tle::parse(line, line2).map(|tle| Tle { name, ..tle }));
            }
            if line.starts_with("2 ") {
                return Some(Err(TleError::Unpaired(TleLine::Two)));
            }
            name = Some(line.strip_prefix("0 ").unwrap_or(line).trim().to_string());
        }
    })
}

/// The TLE checksum of a line: the sum of its digits, counting `-` as 1, mod 10.
/// Only the first 68 columns count.
pub fn checksum(line: &str) -> u8 {
    let sum: u32 = line
        .bytes()
        .take(68)
        .map(|b| match b {
            b'0'..=b'9' => u32::from(b - b'0'),
            b'-' => 1,
            _ => 0,
        })
        .sum();
    (sum % 10) as u8
}

/// Checks the checksum digit in column 69. A line without one passes.
pub fn verify_checksum(line: &str) -> Result<(), TleError> {
    let Some(given) = line.as_bytes().get(68).filter(|b| b.is_ascii_digit()) else {
        return Ok(());
    };
    let expected = given - b'0';
    let computed = checksum(line);
    if expected == computed {
        Ok(())
    } else {
        let line = if line.starts_with('2') {
            TleLine::Two
        } else {
            TleLine::One
        };
        Err(TleError::Checksum {
            line,
            expected,
            computed,
        })
    }
}

/// Returns the line cut or padded to 68 columns, with a correct checksum appended.
pub fn fix_checksum(line: &str) -> String {
    let body: String = line.chars().take(68).collect();
    let body = format!("{body:<68}");
    let sum = checksum(&body);
    format!("{body}{sum}")
}

/// The epoch as a split Julian date, rebuilt from the TLE's year and day so
/// the fraction keeps exactly the eight digits the TLE prints.
fn tle_epoch(year: i32, day_of_year: f64) -> JulianDate {
    let whole_days = day_of_year.floor();
    let fraction = ((day_of_year - whole_days) * 1e8).round() / 1e8;
    // Julian date of the midnight starting January 0; valid 1901-2099.
    let jan0 = f64::from(year * 365 + (year - 1) / 4) + 1_721_044.5;
    JulianDate::new(jan0 + whole_days, fraction)
}

/// Fixed-column access to one TLE line.
struct Columns<'a> {
    text: &'a [u8],
    line: TleLine,
}

impl<'a> Columns<'a> {
    fn new(text: &'a str, line: TleLine) -> Result<Self, TleError> {
        if !text.is_ascii() {
            return Err(TleError::NonAscii(line));
        }
        Ok(Self {
            text: text.trim_end().as_bytes(),
            line,
        })
    }

    fn expect_layout(
        &self,
        number: u8,
        min_len: usize,
        fixed: &[(usize, u8)],
    ) -> Result<(), TleError> {
        let ok = self.text.len() >= min_len
            && self.text[0] == number
            && self.text[1] == b' '
            && fixed.iter().all(|&(i, c)| self.text[i] == c);
        if ok {
            Ok(())
        } else {
            Err(TleError::Layout(self.line))
        }
    }

    /// Columns `range`, clipped to the line's length.
    fn str(&self, range: std::ops::Range<usize>) -> &'a str {
        let end = range.end.min(self.text.len());
        let start = range.start.min(end);
        // The line was checked to be ASCII, so any byte range is valid UTF-8.
        std::str::from_utf8(&self.text[start..end]).expect("ASCII is UTF-8")
    }

    fn field_error(&self, range: std::ops::Range<usize>, field: &'static str) -> TleError {
        TleError::Field {
            line: self.line,
            field,
            text: self.str(range).to_string(),
        }
    }

    fn parse<T: FromStr>(
        &self,
        range: std::ops::Range<usize>,
        field: &'static str,
    ) -> Result<T, TleError> {
        self.str(range.clone())
            .trim()
            .parse()
            .map_err(|_| self.field_error(range, field))
    }

    /// Like `parse`, but a blank or missing field reads as zero.
    fn parse_or_zero<T: FromStr + Default>(
        &self,
        range: std::ops::Range<usize>,
        field: &'static str,
    ) -> Result<T, TleError> {
        if self.str(range.clone()).trim().is_empty() {
            Ok(T::default())
        } else {
            self.parse(range, field)
        }
    }

    /// A field like ` 12345-4`: sign, five mantissa digits with an implied
    /// leading decimal point, then a signed exponent. Reads as 0.12345e-4.
    fn exponent_field(&self, start: usize, field: &'static str) -> Result<f64, TleError> {
        let range = start..start + 8;
        let mantissa = format!(
            "{}.{}",
            self.str(start..start + 1),
            self.str(start + 1..start + 6)
        );
        let mantissa: f64 = mantissa
            .trim()
            .parse()
            .map_err(|_| self.field_error(range.clone(), field))?;
        let exponent: i32 = self.parse(start + 6..start + 8, field)?;
        Ok(mantissa * 10f64.powi(exponent))
    }

    fn catalog_number(&self) -> Result<u32, TleError> {
        alpha5::decode(self.str(2..7)).map_err(|_| self.field_error(2..7, "catalog number"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ISS_1: &str = "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991";
    const ISS_2: &str = "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482";

    #[test]
    fn parses_every_field() {
        let tle = Tle::parse(ISS_1, ISS_2).unwrap();
        assert_eq!(tle.name, None);
        assert_eq!(tle.catalog_number, 25544);
        assert_eq!(tle.classification, 'U');
        assert_eq!(tle.international_designator, "98067A");
        // python-sgp4 gives jdsatepoch 2458826.5, jdsatepochF 0.69339541.
        assert_eq!(tle.epoch, JulianDate::new(2_458_826.5, 0.693_395_41));
        assert_eq!(tle.mean_motion_dot, 0.000_017_64);
        assert_eq!(tle.mean_motion_ddot, 0.0);
        assert!((tle.bstar - 0.387_92e-4).abs() < 1e-20);
        assert_eq!(tle.ephemeris_type, 0);
        assert_eq!(tle.element_set_number, 999);
        assert_eq!(tle.inclination_deg, 51.6439);
        assert_eq!(tle.raan_deg, 211.2001);
        assert_eq!(tle.eccentricity, 0.000_741_7);
        assert_eq!(tle.arg_of_perigee_deg, 17.6667);
        assert_eq!(tle.mean_anomaly_deg, 85.6398);
        assert_eq!(tle.mean_motion, 15.501_034_72);
        assert_eq!(tle.revolution_number, 20248);
    }

    #[test]
    fn keeps_classification() {
        let line1 = fix_checksum(&ISS_1.replacen("25544U", "25544S", 1));
        assert_eq!(Tle::parse(&line1, ISS_2).unwrap().classification, 'S');
    }

    #[test]
    fn decodes_alpha5_catalog_numbers() {
        let line1 = ISS_1.replacen("25544", "A0042", 1);
        let line2 = ISS_2.replacen("25544", "A0042", 1);
        assert_eq!(Tle::parse(&line1, &line2).unwrap().catalog_number, 100_042);
    }

    #[test]
    fn signed_exponent_fields() {
        let line1 = "1 33333U 05037B   05333.02012661  .25992681  00000-0  24476-3 0  1534";
        let tle = Tle::parse(line1, &ISS_2.replacen("25544", "33333", 1)).unwrap();
        assert!((tle.bstar - 0.244_76e-3).abs() < 1e-18);
        assert_eq!(
            tle.epoch.day,
            JulianDate::from_calendar(2005, 11, 29, 0, 0, 0.0).day
        );

        let line1 = "1 16925U 86065D   06151.67415771  .02550794 -30915-6  18784-3 0  4486";
        let tle = Tle::parse(line1, &ISS_2.replacen("25544", "16925", 1)).unwrap();
        assert!((tle.mean_motion_ddot + 0.309_15e-6).abs() < 1e-20);
    }

    #[test]
    fn elements_use_sgp4_units() {
        let e = Tle::parse(ISS_1, ISS_2).unwrap().elements();
        assert!((e.inclination - 51.6439_f64.to_radians()).abs() < 1e-15);
        assert!((e.mean_motion - 15.501_034_72 * TAU / 1440.0).abs() < 1e-15);
    }

    #[test]
    fn rejects_malformed_lines() {
        assert_eq!(
            Tle::parse(&ISS_1[..60], ISS_2),
            Err(TleError::Layout(TleLine::One))
        );
        assert_eq!(
            Tle::parse(ISS_2, ISS_1),
            Err(TleError::Layout(TleLine::One))
        );
        assert_eq!(
            Tle::parse(ISS_1, &ISS_2.replace("51.6439", "51,6439")),
            Err(TleError::Layout(TleLine::Two))
        );
        assert_eq!(
            Tle::parse(&ISS_1.replace("98067A", "98067Ä"), ISS_2),
            Err(TleError::NonAscii(TleLine::One))
        );
        assert!(matches!(
            Tle::parse(&ISS_1.replace("19343", "19x43"), ISS_2),
            Err(TleError::Field {
                field: "epoch day",
                ..
            })
        ));
    }

    #[test]
    fn rejects_mismatched_catalog_numbers() {
        let line2 = ISS_2.replacen("25544", "25545", 1);
        assert_eq!(Tle::parse(ISS_1, &line2), Err(TleError::CatalogMismatch));
    }

    #[test]
    fn parse_all_handles_titles_and_errors() {
        let text = format!(
            "ISS (ZARYA)\n{ISS_1}\n{ISS_2}\n\n{ISS_1}\r\n{ISS_2}\r\n0 TITLED\n{ISS_1}\n{ISS_1}\n{ISS_2}\n"
        );
        let all: Vec<_> = parse_all(&text).collect();
        assert_eq!(all.len(), 4);
        assert_eq!(
            all[0].as_ref().unwrap().name.as_deref(),
            Some("ISS (ZARYA)")
        );
        assert_eq!(all[1].as_ref().unwrap().name, None);
        assert_eq!(all[2], Err(TleError::Unpaired(TleLine::One)));
        assert_eq!(all[3].as_ref().unwrap().name, None);
    }

    #[test]
    fn from_str_takes_exactly_one_tle() {
        let tle: Tle = format!("0 ISS\n{ISS_1}\n{ISS_2}").parse().unwrap();
        assert_eq!(tle.name.as_deref(), Some("ISS"));
        assert!(
            format!("{ISS_1}\n{ISS_2}\n{ISS_1}\n{ISS_2}")
                .parse::<Tle>()
                .is_err()
        );
        assert!("".parse::<Tle>().is_err());
    }

    #[test]
    fn checksums() {
        assert_eq!(checksum(ISS_1), 1);
        assert_eq!(checksum(ISS_2), 2);
        assert_eq!(verify_checksum(ISS_1), Ok(()));
        assert_eq!(verify_checksum(&ISS_1[..68]), Ok(()));
        let corrupted = ISS_1.replace(".00001764", ".90001764");
        assert_eq!(
            verify_checksum(&corrupted),
            Err(TleError::Checksum {
                line: TleLine::One,
                expected: 1,
                computed: 0
            })
        );
        assert_eq!(fix_checksum(&corrupted), format!("{}0", &corrupted[..68]));
        assert_eq!(fix_checksum(&ISS_2[..68]), ISS_2);
    }
}
