# SGPT4 Rust Port

WARNING: Use at your own risk.

This project was an experiment using ChatGPT to convert over a python library for [SGP4](https://pypi.org/project/sgp4/) to rust.  The AI completely generated
all relevant source and test files, except a few manual fixes.

This repo may be broken in weird places, may not have idiomatic rust, and is not in an ideal state from a code quality perspective.  One day it may be though...

As of now this library should be mostly complete.  Feel free to make a PR or raise an issue if there is something missing here or if there are any errors.

## Usage
```rust
    let s = "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991";
    let t = "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482";
    let mut satrec = SatRec::twoline2rv(s, t, "wgs72");

    let jd = 2458826.5;
    let fr = 0.8625;
    let (e, r, v) = satrec.sgp4(jd, fr);
```

## Methodology
I started by breaking the Python SGP4 implementation down into small, independent pieces. For each function or file, I asked ChatGPT to produce a Rust version that stayed as close as possible to the original Python code. The goal wasn’t to get “clever” Rust, but a mostly 1-to-1 port so it would be easy to compare and reason about.

ChatGPT handled most of the translation well - it got the math formulas right the vast majority of the time. The main rough edges were around converting Python idioms and data types into idiomatic Rust, which took a bit of back-and-forth and some manual cleanup on my side.

As the code came together, I added tests that mirror known Python-sgp4 test cases and sanity checks (e.g., r/v against reference values, calendar conversions, etc.) so I could validate behavior as I went and keep the library backward-compatible.

I’m not an orbital mechanics expert (my knowledge is very much beginner level), so I leaned heavily on both the original Python implementation and these test cases to catch small math bugs and sign errors. With enough iteration, prompting, and testing, I was able to get the Rust port to line up closely with the Python library’s results.

## Future
- Demonstrations of functionality.  Visualizations?
- Performance tests.  I wonder how closely it performs to the Python app as well as the legacy C source code.
- Publish to Cargo?