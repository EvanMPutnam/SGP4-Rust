/// Alpha 5 encoding of satellite numbers to fit in 5 characters.
///
/// This is a Rust translation of:
///   - to_alpha5(n)
///   - from_alpha5(s)
///
/// Semantics match the original Python as closely as possible.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Alpha5Error {
    /// Satellite number exceeds 339999, whose Alpha 5 encoding is "Z9999".
    SatelliteNumberTooLarge,
    /// Input string was empty.
    EmptyString,
    /// Failed to parse the numeric portion as an integer.
    ParseIntError,
}

impl std::fmt::Display for Alpha5Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Alpha5Error::SatelliteNumberTooLarge => write!(
                f,
                "satellite number cannot exceed 339999, whose Alpha 5 encoding is 'Z9999'"
            ),
            Alpha5Error::EmptyString => write!(f, "alpha5 string cannot be empty"),
            Alpha5Error::ParseIntError => write!(f, "failed to parse integer"),
        }
    }
}

impl std::error::Error for Alpha5Error {}

/// Encode a satellite number into Alpha 5 format.
///
/// * If `n < 100000`, returns a zero-padded 5-digit decimal string.
/// * If `n > 339999`, returns `Alpha5Error::SatelliteNumberTooLarge`.
/// * Otherwise, returns an Alpha-5 encoded string like `"A0000"`..`"Z9999"`.
pub fn to_alpha5(n: u32) -> Result<String, Alpha5Error> {
    if n < 100_000 {
        return Ok(format!("{:05}", n));
    }
    if n > 339_999 {
        return Err(Alpha5Error::SatelliteNumberTooLarge);
    }

    // Python:
    //   i, n = divmod(n, 10000)
    //   i += ord('A') - 10
    //   if i >= ord('I'): i += 1
    //   if i >= ord('O'): i += 1
    //   return '%c%04d' % (i, n)

    let i = n / 10_000;
    let remainder = n % 10_000;

    let mut code = (i as u8) + b'A' - 10;
    if code >= b'I' {
        code += 1;
    }
    if code >= b'O' {
        code += 1;
    }

    let c = code as char;
    Ok(format!("{c}{remainder:04}"))
}

/// Decode an Alpha 5 string back into a satellite number.
///
/// * If the first character is **not alphabetic**, it decodes the whole
///   string as a plain integer (like the Python `int(s)` path).
/// * If the first character **is alphabetic**, it decodes it as Alpha 5.
pub fn from_alpha5(s: &str) -> Result<u32, Alpha5Error> {
    if s.is_empty() {
        return Err(Alpha5Error::EmptyString);
    }

    let mut chars = s.chars();
    let first = chars.next().unwrap();

    // Python: if not s[0].isalpha(): return int(s)
    if !first.is_ascii_alphabetic() {
        return s
            .parse::<u32>()
            .map_err(|_| Alpha5Error::ParseIntError);
    }

    // Python:
    //   c, s = s[0], s[1:]
    //   n = ord(c) - ord('A') + 10
    //   n -= c > 'I'
    //   n -= c > 'O'
    //   return n * 10000 + int(s)
    let c = first;
    let rest: String = chars.collect();

    let mut n = (c as u8 - b'A') as u32 + 10;
    if c > 'I' {
        n -= 1;
    }
    if c > 'O' {
        n -= 1;
    }

    let suffix = rest
        .parse::<u32>()
        .map_err(|_| Alpha5Error::ParseIntError)?;

    Ok(n * 10_000 + suffix)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_alpha5_numeric() {
        assert_eq!(to_alpha5(0).unwrap(), "00000");
        assert_eq!(to_alpha5(42).unwrap(), "00042");
        assert_eq!(to_alpha5(99999).unwrap(), "99999");
    }

    #[test]
    fn test_to_alpha5_alpha_range() {
        assert_eq!(to_alpha5(100_000).unwrap(), "A0000");
        assert_eq!(to_alpha5(100_001).unwrap(), "A0001");
        assert_eq!(to_alpha5(339_999).unwrap(), "Z9999");
    }

    #[test]
    fn test_to_alpha5_too_large() {
        assert!(matches!(
            to_alpha5(340_000),
            Err(Alpha5Error::SatelliteNumberTooLarge)
        ));
    }

    #[test]
    fn test_from_alpha5_numeric() {
        assert_eq!(from_alpha5("00042").unwrap(), 42);
        assert_eq!(from_alpha5("99999").unwrap(), 99_999);
    }

    #[test]
    fn test_from_alpha5_alpha() {
        assert_eq!(from_alpha5("A0000").unwrap(), 100_000);
        assert_eq!(from_alpha5("Z9999").unwrap(), 339_999);
    }
}