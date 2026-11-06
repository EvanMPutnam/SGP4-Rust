//! Vallado's SGP4 verification suite.
//!
//! `SGP4-VER.TLE` lists 33 element sets, each followed by a start, stop, and
//! step time in minutes. `tcppver.out` is the output of Vallado's C++ code for
//! those times. Both files come from Vallado's 2020 release, as redistributed
//! with python-sgp4. python-sgp4 holds itself to 2e-7 against this file.

use sgp4_rust::{Gravity, OpsMode, Satellite, Tle};

const TOLERANCE: f64 = 2e-7;

/// One row of `tcppver.out`: minutes since epoch, then position and velocity.
type Row = [f64; 7];

/// Each satellite's rows, in file order. Satellite 20413 appears twice,
/// with different time ranges.
fn expected_blocks() -> Vec<(u32, Vec<Row>)> {
    let mut blocks: Vec<(u32, Vec<Row>)> = Vec::new();
    for line in include_str!("data/tcppver.out").lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.get(1) == Some(&"xx") {
            blocks.push((fields[0].parse().unwrap(), Vec::new()));
        } else if !fields.is_empty() {
            let values: Vec<f64> = fields[..7].iter().map(|f| f.parse().unwrap()).collect();
            blocks
                .last_mut()
                .unwrap()
                .1
                .push(values.try_into().unwrap());
        }
    }
    blocks
}

/// The times the C++ test driver visits for one satellite.
fn times(start: f64, stop: f64, step: f64) -> Vec<f64> {
    let mut times = vec![0.0];
    let mut t = start;
    while t <= stop {
        if !(t == 0.0 && start == 0.0) {
            times.push(t);
        }
        t += step;
    }
    if t - stop < step - 1e-6 {
        times.push(stop);
    }
    times
}

#[test]
fn matches_vallado_test_vectors() {
    let mut expected = expected_blocks().into_iter();
    let mut lines = include_str!("data/SGP4-VER.TLE").lines();
    let mut error_codes = Vec::new();
    let mut satellites = 0;
    let mut rows_checked = 0;
    let mut worst: f64 = 0.0;

    while let Some(line1) = lines.next() {
        if !line1.starts_with('1') {
            continue;
        }
        let line2 = lines.next().unwrap();
        let tle = Tle::parse(line1, &line2[..69]).unwrap();
        let schedule: Vec<f64> = line2[69..]
            .split_whitespace()
            .map(|f| f.parse().unwrap())
            .collect();
        let (number, want) = expected.next().unwrap();
        assert_eq!(number, tle.catalog_number);
        satellites += 1;

        let satellite = match Satellite::new(tle.elements(), Gravity::Wgs72, OpsMode::Improved) {
            Ok(satellite) => satellite,
            Err(error) => {
                // The C++ driver prints a placeholder row for a satellite that
                // fails at epoch.
                assert_eq!(
                    want.len(),
                    1,
                    "satellite {} failed at epoch: {error}",
                    tle.catalog_number
                );
                error_codes.push(error.code());
                continue;
            }
        };

        let mut got = Vec::new();
        for t in times(schedule[0], schedule[1], schedule[2]) {
            match satellite.propagate(t) {
                Ok(state) => {
                    let [x, y, z] = state.position;
                    let [vx, vy, vz] = state.velocity;
                    got.push([t, x, y, z, vx, vy, vz]);
                }
                Err(error) => {
                    error_codes.push(error.code());
                    break;
                }
            }
        }

        assert_eq!(
            got.len(),
            want.len(),
            "row count for satellite {}",
            tle.catalog_number
        );
        for (got, want) in got.iter().zip(&want) {
            // Times are printed to 8 decimal places.
            assert!(
                (got[0] - want[0]).abs() < 1e-7,
                "time for satellite {}",
                tle.catalog_number
            );
            for i in 1..7 {
                let diff = (got[i] - want[i]).abs();
                worst = worst.max(diff);
                assert!(
                    diff < TOLERANCE,
                    "satellite {} at t={}: column {i} is {} but tcppver.out has {}",
                    tle.catalog_number,
                    got[0],
                    got[i],
                    want[i]
                );
            }
            rows_checked += 1;
        }
    }

    assert_eq!(satellites, 33);
    assert_eq!(rows_checked, 666);
    // The C++ driver stops each satellite at its first error, in this order.
    assert_eq!(error_codes, [1, 1, 6, 6, 4, 3, 6]);
    eprintln!("checked {rows_checked} rows; largest difference {worst:.2e}");
}

#[test]
fn propagation_order_does_not_matter() {
    // Vallado caches the resonance integrator between calls; this crate does
    // not. Propagating backward, forward, and out of order must agree.
    let tle = Tle::parse(
        "1 09880U 77021A   06176.56157475  .00000421  00000-0  10000-3 0  9814",
        "2 09880  64.5968 349.3786 7069051 270.0229  16.3320  2.00813614112380",
    )
    .unwrap();
    let satellite = Satellite::from_tle(&tle).unwrap();
    let in_order: Vec<_> = (0..=20)
        .map(|i| satellite.propagate(f64::from(i) * 144.0).unwrap())
        .collect();
    for i in (0..=20).rev() {
        assert_eq!(
            satellite.propagate(f64::from(i) * 144.0).unwrap(),
            in_order[i as usize]
        );
    }
}
