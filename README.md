# sgp4-rust

SGP4/SDP4 satellite propagation in Rust, ported from David Vallado's reference
code (2020 revision) by way of Brandon Rhodes's
[python-sgp4](https://pypi.org/project/sgp4/).

It reproduces all 666 rows of Vallado's verification output (`tcppver.out`)
to within 1.2e-7 km, and the same error codes at the same times. Across the
full public CelesTrak catalog (28,179 objects, five times each) it agrees with
python-sgp4 2.27 to within 2.6e-7 km and returns the same errors.

## Usage

```rust
use sgp4_rust::{Geodetic, Satellite, Tle};

let tle = Tle::parse(
    "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991",
    "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482",
)?;
let iss = Satellite::from_tle(&tle)?;

// Minutes after the TLE epoch...
let state = iss.propagate(90.0)?;
// ...or any UTC time, as a chrono::DateTime<Utc> or a JulianDate.
let now = chrono::Utc::now();
let state = iss.propagate_to(now)?;

// TEME position (km) and velocity (km/s), and the point below it.
let ground = Geodetic::from_teme(state.position, now.into());
println!("{:.1} km up", ground.altitude_km);
```

Whole files of two- or three-line TLEs:

```rust
let text = std::fs::read_to_string("catalog.tle")?;
let satellites: Vec<Satellite> = sgp4_rust::tle::parse_all(&text)
    .filter_map(Result::ok)
    .filter_map(|tle| Satellite::from_tle(&tle).ok())
    .collect();
```

`Satellite::propagate` takes `&self`, so one `Satellite` can be shared across
threads (with rayon, say) and propagated to times in any order.

## API

| Item | Purpose |
|---|---|
| `Tle` | A parsed TLE, in the format's own units (degrees, rev/day). `Tle::parse`, `str::parse`, or `tle::parse_all` for files. Alpha-5 catalog numbers are decoded. |
| `Elements` | Mean elements in SGP4 units (radians, rad/min). Build directly for non-TLE sources such as OMM. |
| `Satellite` | An initialized propagator. `new(elements, gravity, mode)` or `from_tle`. |
| `State` | TEME position and velocity. |
| `JulianDate` | A split Julian date that keeps full precision; converts to and from `chrono::DateTime<Utc>`; `gmst()`. |
| `Gravity`, `OpsMode` | WGS-72 (default), WGS-72 old, WGS-84; improved (default) or AFSPC mode. |
| `frames` | `teme_to_ecef` and `Geodetic` (WGS-84 latitude, longitude, altitude). |
| `orbit` | `ClassicalElements::from_state` (Vallado's `rv2coe`) with `Option` for undefined angles. |
| `tle::{checksum, verify_checksum, fix_checksum}`, `alpha5::{encode, decode}` | TLE housekeeping. |
| `TleError`, `PropagationError` | Errors. `PropagationError::code()` gives Vallado's error number. |

## Changes from 0.1

0.2 is a rewrite of the public API. The numerical code keeps Vallado's
structure and variable names so it can be checked line by line against the
reference.

- The `SatRec` struct with 120 public fields, `sgp4init`, `twoline2rv`, and
  `sgp4(jd, fr)` are replaced by `Tle`, `Elements`, and `Satellite`.
- Propagation takes `&self` instead of `&mut self`. The deep-space resonance
  integrator restarts from epoch on each call instead of caching its last
  step; it takes the same steps, so results are identical.
- Nothing panics on bad input. Unknown gravity models, malformed TLEs,
  checksum failures, and propagation failures return errors.
- Errors are enums with values, not integer codes and strings. A decayed
  satellite (error 6) is now an `Err`, as are codes 1–4.
- Angles that `rv2coe` cannot define are `None` instead of 999999.1.

### Bugs fixed

- A TLE's classification (`C`, `S`) was overwritten with `U`.
- The TLE epoch went through a single-`f64` Julian date, losing about 40 µs
  of precision. It is now rebuilt from the TLE's year and day, as python-sgp4
  does.
- `sgp4_array` panicked on any satellite that decayed during the run.
- `sgp4_tsince` returned `Ok` for a decayed satellite but `Err` for other
  failures.
- A TLE that failed at epoch (such as Vallado's case 33334) produced a record
  without an error.
- `from_alpha5` accepted `I`, `O`, and lowercase letters, decoding them to
  wrong numbers.
- After Kepler's equation, sine and cosine were recomputed from the final
  iterate; the reference uses the values from the last loop pass.
- Eccentricity outside `[0, 1)` or a non-positive or NaN mean motion produced
  NaN positions instead of an error.

## Tests

```sh
cargo test
```

`tests/verification.rs` runs Vallado's 33 verification satellites
(`tests/data/SGP4-VER.TLE`) and compares each row with `tests/data/tcppver.out`.
Both files come from Vallado's 2020 release as shipped with python-sgp4.
