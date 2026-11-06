//! Alpha-5 catalog numbers.
//!
//! TLEs have five columns for the catalog number. Alpha-5 extends the range to
//! 339,999 by replacing the leading digit with a letter: `A` is 10, `B` is 11,
//! and so on, skipping `I` and `O` because they look like `1` and `0`.

use std::fmt;

/// The largest number Alpha-5 can encode, written `Z9999`.
pub const MAX: u32 = 339_999;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Alpha5Error {
    /// The number is above [`MAX`].
    TooLarge(u32),
    /// The text is not five digits, or a letter (other than `I` or `O`)
    /// followed by four digits.
    Invalid(String),
}

impl fmt::Display for Alpha5Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Alpha5Error::TooLarge(n) => {
                write!(f, "catalog number {n} is above the Alpha-5 limit of {MAX}")
            }
            Alpha5Error::Invalid(s) => write!(f, "{s:?} is not an Alpha-5 catalog number"),
        }
    }
}

impl std::error::Error for Alpha5Error {}

/// Encodes a catalog number in five characters.
///
/// ```
/// assert_eq!(sgp4_rust::alpha5::encode(25544).unwrap(), "25544");
/// assert_eq!(sgp4_rust::alpha5::encode(100_000).unwrap(), "A0000");
/// ```
pub fn encode(number: u32) -> Result<String, Alpha5Error> {
    if number > MAX {
        return Err(Alpha5Error::TooLarge(number));
    }
    if number < 100_000 {
        return Ok(format!("{number:05}"));
    }
    let mut letter = b'A' + (number / 10_000 - 10) as u8;
    if letter >= b'I' {
        letter += 1;
    }
    if letter >= b'O' {
        letter += 1;
    }
    Ok(format!("{}{:04}", letter as char, number % 10_000))
}

/// Decodes a catalog number. Leading and trailing spaces are ignored, and
/// purely numeric text of any length up to five digits is accepted.
///
/// ```
/// assert_eq!(sgp4_rust::alpha5::decode("Z9999").unwrap(), 339_999);
/// assert_eq!(sgp4_rust::alpha5::decode("    5").unwrap(), 5);
/// ```
pub fn decode(text: &str) -> Result<u32, Alpha5Error> {
    let invalid = || Alpha5Error::Invalid(text.to_string());
    let trimmed = text.trim();
    let mut chars = trimmed.chars();
    let first = chars.next().ok_or_else(invalid)?;
    let rest = chars.as_str();

    if first.is_ascii_digit() {
        if trimmed.len() > 5 || !rest.bytes().all(|b| b.is_ascii_digit()) {
            return Err(invalid());
        }
        return trimmed.parse().map_err(|_| invalid());
    }

    if !first.is_ascii_uppercase()
        || first == 'I'
        || first == 'O'
        || rest.len() != 4
        || !rest.bytes().all(|b| b.is_ascii_digit())
    {
        return Err(invalid());
    }
    let mut prefix = u32::from(first as u8 - b'A') + 10;
    if first > 'I' {
        prefix -= 1;
    }
    if first > 'O' {
        prefix -= 1;
    }
    let suffix: u32 = rest.parse().map_err(|_| invalid())?;
    Ok(prefix * 10_000 + suffix)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_range_is_zero_padded() {
        assert_eq!(encode(0).unwrap(), "00000");
        assert_eq!(encode(42).unwrap(), "00042");
        assert_eq!(encode(99_999).unwrap(), "99999");
    }

    #[test]
    fn letters_skip_i_and_o() {
        assert_eq!(encode(100_000).unwrap(), "A0000");
        assert_eq!(encode(170_000).unwrap(), "H0000");
        assert_eq!(encode(180_000).unwrap(), "J0000");
        assert_eq!(encode(230_000).unwrap(), "P0000");
        assert_eq!(encode(MAX).unwrap(), "Z9999");
        assert_eq!(encode(MAX + 1), Err(Alpha5Error::TooLarge(MAX + 1)));
    }

    #[test]
    fn every_encoded_number_decodes_back() {
        for n in (0..=MAX).step_by(997).chain([MAX]) {
            assert_eq!(decode(&encode(n).unwrap()).unwrap(), n);
        }
    }

    #[test]
    fn rejects_bad_text() {
        for bad in [
            "", "I0000", "O1234", "a0000", "A123", "A12345", "12a45", "123456", "-1234",
        ] {
            assert!(decode(bad).is_err(), "{bad:?} should be rejected");
        }
    }
}
