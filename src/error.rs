//! Error types.

use std::fmt;

/// Which line of a TLE an error refers to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TleLine {
    One,
    Two,
}

impl fmt::Display for TleLine {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            TleLine::One => "line 1",
            TleLine::Two => "line 2",
        })
    }
}

/// A TLE could not be parsed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TleError {
    /// The line contains non-ASCII characters.
    NonAscii(TleLine),
    /// The line is too short, or its fixed columns (spaces, decimal points,
    /// line number) are not where the TLE format puts them.
    Layout(TleLine),
    /// A field holds text that is not a valid value.
    Field {
        line: TleLine,
        field: &'static str,
        text: String,
    },
    /// The two lines carry different catalog numbers.
    CatalogMismatch,
    /// The checksum digit does not match the line's contents.
    Checksum {
        line: TleLine,
        expected: u8,
        computed: u8,
    },
    /// A line 1 was not followed by a line 2, or vice versa.
    Unpaired(TleLine),
}

impl fmt::Display for TleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TleError::NonAscii(line) => write!(f, "TLE {line} contains non-ASCII characters"),
            TleError::Layout(line) => write!(f, "TLE {line} does not follow the TLE column layout"),
            TleError::Field { line, field, text } => {
                write!(f, "TLE {line} has an invalid {field}: {text:?}")
            }
            TleError::CatalogMismatch => {
                f.write_str("TLE lines 1 and 2 have different catalog numbers")
            }
            TleError::Checksum {
                line,
                expected,
                computed,
            } => write!(
                f,
                "TLE {line} gives checksum {expected} but its contents sum to {computed}"
            ),
            TleError::Unpaired(line) => write!(f, "TLE {line} is missing its other line"),
        }
    }
}

impl std::error::Error for TleError {}

/// SGP4 could not produce a state. Each variant carries the offending value.
///
/// [`PropagationError::code`] gives the matching error number from Vallado's
/// reference code.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PropagationError {
    /// Mean eccentricity left `[0, 1)`. Code 1.
    Eccentricity(f64),
    /// Mean motion is not positive. Code 2.
    MeanMotion(f64),
    /// Eccentricity after lunar-solar perturbations left `[0, 1]`. Code 3.
    PerturbedEccentricity(f64),
    /// Semi-latus rectum is negative. Code 4.
    SemiLatusRectum(f64),
    /// The orbit radius, in Earth radii, is below 1: the satellite has decayed. Code 6.
    Decayed(f64),
}

impl PropagationError {
    /// The error number Vallado's `sgp4` stores in `satrec.error`.
    pub fn code(self) -> u8 {
        match self {
            PropagationError::Eccentricity(_) => 1,
            PropagationError::MeanMotion(_) => 2,
            PropagationError::PerturbedEccentricity(_) => 3,
            PropagationError::SemiLatusRectum(_) => 4,
            PropagationError::Decayed(_) => 6,
        }
    }
}

impl fmt::Display for PropagationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PropagationError::Eccentricity(e) => {
                write!(
                    f,
                    "mean eccentricity {e} is not within range 0.0 <= e < 1.0"
                )
            }
            PropagationError::MeanMotion(n) => write!(f, "mean motion {n} is not positive"),
            PropagationError::PerturbedEccentricity(e) => {
                write!(
                    f,
                    "perturbed eccentricity {e} is not within range 0.0 <= e <= 1.0"
                )
            }
            PropagationError::SemiLatusRectum(p) => {
                write!(f, "semi-latus rectum {p} is less than zero")
            }
            PropagationError::Decayed(r) => write!(
                f,
                "orbit radius {r} Earth radii is below the surface: the satellite has decayed"
            ),
        }
    }
}

impl std::error::Error for PropagationError {}
