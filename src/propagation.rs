// Direct-ish translation of python-sgp4's core SGP4 implementation.
// The math and flow closely follow Vallado's C++ code and Brandon Rhodes's
// Python port, with minimal "clever" Rust-ification.

use crate::ext::{invjday, jday};
use crate::io::twoline2satrec;
use log::debug;
use std::f64::consts::PI;

pub const DEG2RAD: f64 = PI / 180.0;
pub const TWOPI: f64 = 2.0 * PI;

/// Wrap a floating-point angle to [0, 2π).
#[inline]
fn wrap_to_2pi(x: f64) -> f64 {
    let mut v = x % TWOPI;
    if v < 0.0 {
        v += TWOPI;
    }
    v
}

/// Gravity constants, corresponds to Python getgravconst().
#[derive(Debug, Clone, Copy)]
pub struct GravConst {
    pub tumin: f64,
    pub mu: f64,
    pub radiusearthkm: f64,
    pub xke: f64,
    pub j2: f64,
    pub j3: f64,
    pub j4: f64,
    pub j3oj2: f64,
}

/// Port of python `getgravconst(whichconst)`.
pub fn get_grav_const(whichconst: &str) -> GravConst {
    match whichconst {
        "wgs72old" => {
            let mu = 398600.79964_f64;
            let radiusearthkm = 6378.135_f64;
            let xke = 0.074_366_916_1_f64;
            let tumin = 1.0 / xke;
            let j2 = 0.001_082_616_f64;
            let j3 = -0.000_002_538_81_f64;
            let j4 = -0.000_001_655_97_f64;
            let j3oj2 = j3 / j2;
            GravConst {
                tumin,
                mu,
                radiusearthkm,
                xke,
                j2,
                j3,
                j4,
                j3oj2,
            }
        }
        "wgs72" => {
            let mu = 398600.8_f64;
            let radiusearthkm = 6378.135_f64;
            let xke = 60.0 / (radiusearthkm * radiusearthkm * radiusearthkm / mu).sqrt();
            let tumin = 1.0 / xke;
            let j2 = 0.001_082_616_f64;
            let j3 = -0.000_002_538_81_f64;
            let j4 = -0.000_001_655_97_f64;
            let j3oj2 = j3 / j2;
            GravConst {
                tumin,
                mu,
                radiusearthkm,
                xke,
                j2,
                j3,
                j4,
                j3oj2,
            }
        }
        "wgs84" => {
            let mu = 398600.5_f64;
            let radiusearthkm = 6378.137_f64;
            let xke = 60.0 / (radiusearthkm * radiusearthkm * radiusearthkm / mu).sqrt();
            let tumin = 1.0 / xke;
            let j2 = 0.001_082_629_989_05_f64;
            let j3 = -0.000_002_532_153_06_f64;
            let j4 = -0.000_001_610_987_61_f64;
            let j3oj2 = j3 / j2;
            GravConst {
                tumin,
                mu,
                radiusearthkm,
                xke,
                j2,
                j3,
                j4,
                j3oj2,
            }
        }
        other => panic!("Unknown whichconst: {other} (use wgs72old, wgs72, wgs84)"),
    }
}

/// Vallado gstime(jdut1) as in the Python code.
pub fn gstime(jdut1: f64) -> f64 {
    let tut1 = (jdut1 - 2_451_545.0) / 36_525.0;
    let mut temp = -6.2e-6 * tut1 * tut1 * tut1
        + 0.093_104 * tut1 * tut1
        + (876_600.0 * 3600.0 + 8_640_184.812_866) * tut1
        + 67_310.548_41;
    // convert to radians:
    // 360/86400 = 1/240, so divide by 240 to go from sec -> deg, then deg->rad
    temp = (temp * DEG2RAD / 240.0) % TWOPI;
    if temp < 0.0 {
        temp += TWOPI;
    }
    temp
}

/// Old name for compatibility with the Python code style.
pub fn _gstime(jdut1: f64) -> f64 {
    gstime(jdut1)
}

/// Big port of python-sgp4's satrec object fields.
/// This is intentionally verbose and close to the original attributes.
#[derive(Debug, Clone)]
pub struct SatRec {
    // Earth constants
    pub tumin: f64,
    pub mu: f64,
    pub radiusearthkm: f64,
    pub xke: f64,
    pub j2: f64,
    pub j3: f64,
    pub j4: f64,
    pub j3oj2: f64,

    // General state / identification
    pub error: i32,
    pub error_message: Option<String>,
    pub operationmode: char,
    pub satnum_str: String,
    pub classification: char,

    // Input orbital elements
    pub bstar: f64,
    pub ndot: f64,
    pub nddot: f64,
    pub ecco: f64,
    pub argpo: f64,
    pub inclo: f64,
    pub mo: f64,
    pub no_kozai: f64,
    pub nodeo: f64,

    // Single averaged mean elements
    pub am: f64,
    pub em: f64,
    pub im: f64,
    pub Om: f64,
    pub om: f64,
    pub mm: f64,
    pub nm: f64,

    // Near-earth variables
    pub isimp: i32,
    pub method: char,
    pub aycof: f64,
    pub con41: f64,
    pub cc1: f64,
    pub cc4: f64,
    pub cc5: f64,
    pub d2: f64,
    pub d3: f64,
    pub d4: f64,
    pub delmo: f64,
    pub eta: f64,
    pub argpdot: f64,
    pub omgcof: f64,
    pub sinmao: f64,
    pub t: f64,
    pub t2cof: f64,
    pub t3cof: f64,
    pub t4cof: f64,
    pub t5cof: f64,
    pub x1mth2: f64,
    pub x7thm1: f64,
    pub mdot: f64,
    pub nodedot: f64,
    pub xlcof: f64,
    pub xmcof: f64,
    pub nodecf: f64,

    // Deep-space variables
    pub irez: i32,
    pub d2201: f64,
    pub d2211: f64,
    pub d3210: f64,
    pub d3222: f64,
    pub d4410: f64,
    pub d4422: f64,
    pub d5220: f64,
    pub d5232: f64,
    pub d5421: f64,
    pub d5433: f64,
    pub dedt: f64,
    pub del1: f64,
    pub del2: f64,
    pub del3: f64,
    pub didt: f64,
    pub dmdt: f64,
    pub dnodt: f64,
    pub domdt: f64,
    pub e3: f64,
    pub ee2: f64,
    pub peo: f64,
    pub pgho: f64,
    pub pho: f64,
    pub pinco: f64,
    pub plo: f64,
    pub se2: f64,
    pub se3: f64,
    pub sgh2: f64,
    pub sgh3: f64,
    pub sgh4: f64,
    pub sh2: f64,
    pub sh3: f64,
    pub si2: f64,
    pub si3: f64,
    pub sl2: f64,
    pub sl3: f64,
    pub sl4: f64,
    pub gsto: f64,
    pub xfact: f64,
    pub xgh2: f64,
    pub xgh3: f64,
    pub xgh4: f64,
    pub xh2: f64,
    pub xh3: f64,
    pub xi2: f64,
    pub xi3: f64,
    pub xl2: f64,
    pub xl3: f64,
    pub xl4: f64,
    pub xlamo: f64,
    pub zmol: f64,
    pub zmos: f64,
    pub atime: f64,
    pub xli: f64,
    pub xni: f64,

    // Derived quantities
    pub no_unkozai: f64,
    pub a: f64,
    pub alta: f64,
    pub altp: f64,
    pub init: char,

    pub epochyr: i32,
    pub epochdays: f64,
    pub jdsatepoch: f64,
    pub jdsatepoch_f: f64,
}

const MINUTES_PER_DAY: f64 = 1440.0;

impl SatRec {
    /// Convenience accessor mirroring the Python `no` property.
    /// In Vallado's code `no` is just `no_kozai`.
    pub fn no(&self) -> f64 {
        self.no_kozai
    }

    /// Construct a `SatRec` from two TLE lines, similar to Python
    /// `Satrec.twoline2rv()`.
    ///
    /// `whichconst` is typically "wgs72", "wgs84", or "wgs72old".
    /// This uses opsmode `'i'` (improved) by default.
    pub fn twoline2rv(line1: &str, line2: &str, whichconst: &str) -> SatRec {
        // This is the Rust port of Vallado's `twoline2rv` you already have.
        // It fills in all the fields and calls `sgp4init`.
        twoline2satrec(line1, line2, whichconst, 'i').expect("Cannot parse satrec")
    }

    /// Instance-style wrapper around the global `sgp4init()` function.
    ///
    /// This matches the Python `Satrec.sgp4init()` signature.
    /// `epoch` is in **days since 1950-01-0 00:00 UT** (Vallado SGP4 epoch).
    /// Instance-style wrapper around the global `sgp4init()` function.
    /// Matches Python `Satrec.sgp4init()`.
    /// `epoch` is in days since 1950-01-0 00:00 UT.
    pub fn sgp4init(
        &mut self,
        whichconst: &str,
        opsmode: char,
        satnum: &str,
        epoch: f64,
        bstar: f64,
        ndot: f64,
        nddot: f64,
        ecco: f64,
        argpo: f64,
        inclo: f64,
        mo: f64,
        no_kozai: f64,
        nodeo: f64,
    ) -> bool {
        // Split epoch into whole + fractional days since 1950-01-0
        let whole = epoch.floor();
        let mut fraction = epoch - whole;
        let whole_jd = whole + 2_433_281.5_f64; // 2433281.5

        // If epoch looks like a TLE-style decimal, respect 1e-8 precision
        if (epoch * 1.0e8).round() == epoch * 1.0e8 {
            fraction = (fraction * 1.0e8).round() / 1.0e8;
        }

        // *** Key difference vs your current code: do NOT add fraction into jdsatepoch ***
        self.jdsatepoch = whole_jd; // integer JD part
        self.jdsatepoch_f = fraction; // fractional part

        // Fill epoch year/day-of-year like Python
        let (y, m, d, hr, min, sec) = invjday(whole_jd);
        let jan0 = jday(y, 1, 0, 0, 0, 0.0);
        self.epochyr = (y % 100) as i32;
        self.epochdays = (whole_jd - jan0) + fraction;

        // Basic bookkeeping
        self.classification = 'U';
        self.operationmode = opsmode;
        self.satnum_str = satnum.to_string();

        // Copy element values
        self.bstar = bstar;
        self.ndot = ndot;
        self.nddot = nddot;
        self.ecco = ecco;
        self.argpo = argpo;
        self.inclo = inclo;
        self.mo = mo;
        self.no_kozai = no_kozai;
        self.nodeo = nodeo;

        // Call the core Vallado initializer (same as Python)
        sgp4init(
            whichconst, opsmode, satnum, epoch, bstar, ndot, nddot, ecco, argpo, inclo, mo,
            no_kozai, nodeo, self,
        )
    }

    /// Propagate to a given Julian Date / fractional day, mirroring
    /// Python `Satrec.sgp4(jd, fr)`.
    ///
    /// Requires that `self.jdsatepoch` and `self.jdsatepoch_f` have been
    /// set (by TLE parsing or `sgp4init()`).
    ///
    /// Returns `(error_code, r_km, v_km_s)`.
    pub fn sgp4(&mut self, jd: f64, fr: f64) -> (i32, [f64; 3], [f64; 3]) {
        let tsince = ((jd - self.jdsatepoch) * MINUTES_PER_DAY)
            + ((fr - self.jdsatepoch_f) * MINUTES_PER_DAY);

        self.sgp4_tsince(tsince)
    }

    /// Propagate by minutes since epoch, mirroring Python `sgp4_tsince()`.
    ///
    /// Returns `(error_code, r_km, v_km_s)`.
    pub fn sgp4_tsince(&mut self, tsince: f64) -> (i32, [f64; 3], [f64; 3]) {
        let (r, v) = sgp4(self, tsince);
        match (r, v) {
            (Some(r), Some(v)) => (self.error, r, v),
            _ => panic!("sgp4 error {}: {:?}", self.error, self.error_message),
        }
    }

    /// Vectorized propagation similar to Python `sgp4_array()`, but using
    /// plain `Vec` instead of NumPy.
    ///
    /// `jd` and `fr` must be the same length.
    ///
    /// Returns `(errors, positions_km, velocities_km_s)`:
    /// * `errors`: one error code per epoch
    /// * `positions_km`: `len × 3` position vectors
    /// * `velocities_km_s`: `len × 3` velocity vectors
    pub fn sgp4_array(
        &mut self,
        jd: &[f64],
        fr: &[f64],
    ) -> (Vec<i32>, Vec<[f64; 3]>, Vec<[f64; 3]>) {
        assert_eq!(jd.len(), fr.len(), "jd and fr must have the same length");

        let mut errors = Vec::with_capacity(jd.len());
        let mut rs = Vec::with_capacity(jd.len());
        let mut vs = Vec::with_capacity(jd.len());

        for (&jd_i, &fr_i) in jd.iter().zip(fr.iter()) {
            let (e, r, v) = self.sgp4(jd_i, fr_i);
            errors.push(e);
            rs.push(r);
            vs.push(v);
        }

        (errors, rs, vs)
    }
}

impl Default for SatRec {
    fn default() -> Self {
        SatRec {
            tumin: 0.0,
            mu: 0.0,
            radiusearthkm: 0.0,
            xke: 0.0,
            j2: 0.0,
            j3: 0.0,
            j4: 0.0,
            j3oj2: 0.0,
            error: 0,
            error_message: None,
            operationmode: 'i',
            satnum_str: String::new(),
            classification: 'U',
            bstar: 0.0,
            ndot: 0.0,
            nddot: 0.0,
            ecco: 0.0,
            argpo: 0.0,
            inclo: 0.0,
            mo: 0.0,
            no_kozai: 0.0,
            nodeo: 0.0,
            am: 0.0,
            em: 0.0,
            im: 0.0,
            Om: 0.0,
            om: 0.0,
            mm: 0.0,
            nm: 0.0,
            isimp: 0,
            method: 'n',
            aycof: 0.0,
            con41: 0.0,
            cc1: 0.0,
            cc4: 0.0,
            cc5: 0.0,
            d2: 0.0,
            d3: 0.0,
            d4: 0.0,
            delmo: 0.0,
            eta: 0.0,
            argpdot: 0.0,
            omgcof: 0.0,
            sinmao: 0.0,
            t: 0.0,
            t2cof: 0.0,
            t3cof: 0.0,
            t4cof: 0.0,
            t5cof: 0.0,
            x1mth2: 0.0,
            x7thm1: 0.0,
            mdot: 0.0,
            nodedot: 0.0,
            xlcof: 0.0,
            xmcof: 0.0,
            nodecf: 0.0,
            irez: 0,
            d2201: 0.0,
            d2211: 0.0,
            d3210: 0.0,
            d3222: 0.0,
            d4410: 0.0,
            d4422: 0.0,
            d5220: 0.0,
            d5232: 0.0,
            d5421: 0.0,
            d5433: 0.0,
            dedt: 0.0,
            del1: 0.0,
            del2: 0.0,
            del3: 0.0,
            didt: 0.0,
            dmdt: 0.0,
            dnodt: 0.0,
            domdt: 0.0,
            e3: 0.0,
            ee2: 0.0,
            peo: 0.0,
            pgho: 0.0,
            pho: 0.0,
            pinco: 0.0,
            plo: 0.0,
            se2: 0.0,
            se3: 0.0,
            sgh2: 0.0,
            sgh3: 0.0,
            sgh4: 0.0,
            sh2: 0.0,
            sh3: 0.0,
            si2: 0.0,
            si3: 0.0,
            sl2: 0.0,
            sl3: 0.0,
            sl4: 0.0,
            gsto: 0.0,
            xfact: 0.0,
            xgh2: 0.0,
            xgh3: 0.0,
            xgh4: 0.0,
            xh2: 0.0,
            xh3: 0.0,
            xi2: 0.0,
            xi3: 0.0,
            xl2: 0.0,
            xl3: 0.0,
            xl4: 0.0,
            xlamo: 0.0,
            zmol: 0.0,
            zmos: 0.0,
            atime: 0.0,
            xli: 0.0,
            xni: 0.0,
            no_unkozai: 0.0,
            a: 0.0,
            alta: 0.0,
            altp: 0.0,
            init: 'n',
            epochyr: 0,
            epochdays: 0.0,
            jdsatepoch: 0.0,
            jdsatepoch_f: 0.0,
        }
    }
}

// -----------------------------------------------------------------------------
// Deep-space periodic terms: dpper
// -----------------------------------------------------------------------------
fn dpper(
    satrec: &SatRec, // or &mut SatRec if you mutate it
    inclo: f64,
    init: char, // 'y' / 'n'
    mut ep: f64,
    mut inclp: f64,
    mut nodep: f64,
    mut argpp: f64,
    mut mp: f64,
    opsmode: char, // 'a' or 'i'
) -> (f64, f64, f64, f64, f64) {
    let e3 = satrec.e3;
    let ee2 = satrec.ee2;
    let peo = satrec.peo;
    let pgho = satrec.pgho;
    let pho = satrec.pho;
    let pinco = satrec.pinco;
    let plo = satrec.plo;
    let se2 = satrec.se2;
    let se3 = satrec.se3;
    let sgh2 = satrec.sgh2;
    let sgh3 = satrec.sgh3;
    let sgh4 = satrec.sgh4;
    let sh2 = satrec.sh2;
    let sh3 = satrec.sh3;
    let si2 = satrec.si2;
    let si3 = satrec.si3;
    let sl2 = satrec.sl2;
    let sl3 = satrec.sl3;
    let sl4 = satrec.sl4;
    let t = satrec.t;
    let xgh2 = satrec.xgh2;
    let xgh3 = satrec.xgh3;
    let xgh4 = satrec.xgh4;
    let xh2 = satrec.xh2;
    let xh3 = satrec.xh3;
    let xi2 = satrec.xi2;
    let xi3 = satrec.xi3;
    let xl2 = satrec.xl2;
    let xl3 = satrec.xl3;
    let xl4 = satrec.xl4;
    let zmol = satrec.zmol;
    let zmos = satrec.zmos;

    let zns = 1.19459e-5_f64;
    let zes = 0.01675_f64;
    let znl = 1.5835218e-4_f64;
    let zel = 0.05490_f64;

    // --- solar terms ---
    let mut zm = zmos + zns * t;
    if init == 'y' {
        zm = zmos;
    }
    let mut zf = zm + 2.0 * zes * zm.sin();
    let mut sinzf = zf.sin();
    let mut f2 = 0.5 * sinzf * sinzf - 0.25;
    let mut f3 = -0.5 * sinzf * zf.cos();
    let ses = se2 * f2 + se3 * f3;
    let sis = si2 * f2 + si3 * f3;
    let sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
    let sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
    let shs = sh2 * f2 + sh3 * f3;

    // --- lunar terms ---
    zm = zmol + znl * t;
    if init == 'y' {
        zm = zmol;
    }
    zf = zm + 2.0 * zel * zm.sin();
    sinzf = zf.sin();
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * zf.cos();
    let sel = ee2 * f2 + e3 * f3;
    let sil = xi2 * f2 + xi3 * f3;
    let sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
    let sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
    let shll = xh2 * f2 + xh3 * f3;

    let mut pe = ses + sel;
    let mut pinc = sis + sil;
    let mut pl = sls + sll;
    let mut pgh = sghs + sghl;
    let mut ph = shs + shll;

    if init == 'n' {
        pe -= peo;
        pinc -= pinco;
        pl -= plo;
        pgh -= pgho;
        ph -= pho;
        inclp += pinc;
        ep += pe;
        let sinip = inclp.sin();
        let cosip = inclp.cos();

        if inclp >= 0.2 {
            let mut ph_adj = ph / sinip;
            let mut pgh_adj = pgh - cosip * ph_adj;
            argpp += pgh_adj;
            nodep += ph_adj;
            mp += pl;
        } else {
            let sinop = nodep.sin();
            let cosop = nodep.cos();
            let mut alfdp = sinip * sinop;
            let mut betdp = sinip * cosop;
            let dalf = ph * cosop + pinc * cosip * sinop;
            let dbet = -ph * sinop + pinc * cosip * cosop;
            alfdp += dalf;
            betdp += dbet;

            // wrap nodep
            nodep = if nodep >= 0.0 {
                nodep % TWOPI
            } else {
                -((-nodep) % TWOPI)
            };

            if nodep < 0.0 && opsmode == 'a' {
                nodep += TWOPI;
            }

            let xls = mp + argpp + pl + pgh + (cosip - pinc * sinip) * nodep;
            let xnoh = nodep;
            nodep = alfdp.atan2(betdp);

            if nodep < 0.0 && opsmode == 'a' {
                nodep += TWOPI;
            }

            if (xnoh - nodep).abs() > PI {
                if nodep < xnoh {
                    nodep += TWOPI;
                } else {
                    nodep -= TWOPI;
                }
            }

            mp += pl;
            argpp = xls - mp - cosip * nodep;
        }
    }

    (ep, inclp, nodep, argpp, mp)
}

// -----------------------------------------------------------------------------
// dscom - deep-space common items
// -----------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn dscom(
    epoch: f64,
    ep: f64,
    argpp: f64,
    tc: f64,
    inclp: f64,
    nodep: f64,
    np: f64,
) -> (
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let zes = 0.01675_f64;
    let zel = 0.05490_f64;
    let c1ss = 2.986_479_7e-6_f64;
    let c1l = 4.796_806_5e-7_f64;
    let zsinis = 0.397_854_16_f64;
    let zcosis = 0.917_448_67_f64;
    let zcosgs = 0.194_590_5_f64;
    let zsings = -0.980_884_58_f64;

    let mut nm = np;
    let mut em = ep;
    let snodm = nodep.sin();
    let cnodm = nodep.cos();
    let sinomm = argpp.sin();
    let cosomm = argpp.cos();
    let sinim = inclp.sin();
    let cosim = inclp.cos();
    let emsq = em * em;
    let betasq = 1.0 - emsq;
    let rtemsq = betasq.sqrt();

    // initialize lunar solar terms
    let mut peo = 0.0;
    let mut pinco = 0.0;
    let mut plo = 0.0;
    let mut pgho = 0.0;
    let mut pho = 0.0;
    let day = epoch + 18_261.5 + tc / 1440.0;
    let mut xnodce = (4.5236020 - 9.242_202_9e-4 * day) % TWOPI;
    let stem = xnodce.sin();
    let ctem = xnodce.cos();
    let mut zcosil = 0.913_751_64 - 0.035_680_96 * ctem;
    let mut zsinil = (1.0 - zcosil * zcosil).sqrt();
    let mut zsinhl = 0.089_683_511 * stem / zsinil;
    let mut zcoshl = (1.0 - zsinhl * zsinhl).sqrt();
    let mut gam = 5.835_151_4 + 0.001_944_368_0 * day;
    let mut zx = 0.397_854_16 * stem / zsinil;
    let mut zy = zcoshl * ctem + 0.917_448_67 * zsinhl * stem;
    zx = zx.atan2(zy);
    zx = gam + zx - xnodce;
    let mut zcosgl = zx.cos();
    let mut zsingl = zx.sin();

    let mut zcosg = zcosgs;
    let mut zsing = zsings;
    let mut zcosi = zcosis;
    let mut zsini = zsinis;
    let mut zcosh = cnodm;
    let mut zsinh = snodm;
    let mut cc = c1ss;
    let xnoi = 1.0 / nm;

    let mut s1 = 0.0;
    let mut s2 = 0.0;
    let mut s3 = 0.0;
    let mut s4 = 0.0;
    let mut s5 = 0.0;
    let mut s6 = 0.0;
    let mut s7 = 0.0;

    let mut ss1 = 0.0;
    let mut ss2 = 0.0;
    let mut ss3 = 0.0;
    let mut ss4 = 0.0;
    let mut ss5 = 0.0;
    let mut ss6 = 0.0;
    let mut ss7 = 0.0;

    let mut sz1 = 0.0;
    let mut sz2 = 0.0;
    let mut sz3 = 0.0;
    let mut sz11 = 0.0;
    let mut sz12 = 0.0;
    let mut sz13 = 0.0;
    let mut sz21 = 0.0;
    let mut sz22 = 0.0;
    let mut sz23 = 0.0;
    let mut sz31 = 0.0;
    let mut sz32 = 0.0;
    let mut sz33 = 0.0;

    let mut z1 = 0.0;
    let mut z2 = 0.0;
    let mut z3 = 0.0;
    let mut z11 = 0.0;
    let mut z12 = 0.0;
    let mut z13 = 0.0;
    let mut z21 = 0.0;
    let mut z22 = 0.0;
    let mut z23 = 0.0;
    let mut z31 = 0.0;
    let mut z32 = 0.0;
    let mut z33 = 0.0;

    let mut e3 = 0.0;
    let mut ee2 = 0.0;
    let mut se2 = 0.0;
    let mut se3 = 0.0;
    let mut sgh2 = 0.0;
    let mut sgh3 = 0.0;
    let mut sgh4 = 0.0;
    let mut sh2 = 0.0;
    let mut sh3 = 0.0;
    let mut si2 = 0.0;
    let mut si3 = 0.0;
    let mut sl2 = 0.0;
    let mut sl3 = 0.0;
    let mut sl4 = 0.0;
    let mut xgh2 = 0.0;
    let mut xgh3 = 0.0;
    let mut xgh4 = 0.0;
    let mut xh2 = 0.0;
    let mut xh3 = 0.0;
    let mut xi2 = 0.0;
    let mut xi3 = 0.0;
    let mut xl2 = 0.0;
    let mut xl3 = 0.0;
    let mut xl4 = 0.0;
    let mut zmol = 0.0;
    let mut zmos = 0.0;

    for lsflg in 1..=2 {
        let a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        let a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        let a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        let a8 = zsing * zsini;
        let a9 = zsing * zsinh + zcosg * zcosi * zcosh;
        let a10 = zcosg * zsini;
        let a2 = cosim * a7 + sinim * a8;
        let a4 = cosim * a9 + sinim * a10;
        let a5 = -sinim * a7 + cosim * a8;
        let a6 = -sinim * a9 + cosim * a10;

        let x1 = a1 * cosomm + a2 * sinomm;
        let x2 = a3 * cosomm + a4 * sinomm;
        let x3 = -a1 * sinomm + a2 * cosomm;
        let x4 = -a3 * sinomm + a4 * cosomm;
        let x5 = a5 * sinomm;
        let x6 = a6 * sinomm;
        let x7 = a5 * cosomm;
        let x8 = a6 * cosomm;

        z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
        z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * emsq;
        z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * emsq;
        z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * emsq;
        z11 = -6.0 * a1 * a5 + emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
        z12 = -6.0 * (a1 * a6 + a3 * a5)
            + emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        z21 = 6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
        z22 = 6.0 * (a4 * a5 + a2 * a6)
            + emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
        z23 = 6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
        z1 = z1 + z1 + betasq * z31;
        z2 = z2 + z2 + betasq * z32;
        z3 = z3 + z3 + betasq * z33;
        s3 = cc * xnoi;
        s2 = -0.5 * s3 / rtemsq;
        s4 = s3 * rtemsq;
        s1 = -15.0 * em * s4;
        s5 = x1 * x3 + x2 * x4;
        s6 = x2 * x3 + x1 * x4;
        s7 = x2 * x4 - x1 * x3;

        if lsflg == 1 {
            ss1 = s1;
            ss2 = s2;
            ss3 = s3;
            ss4 = s4;
            ss5 = s5;
            ss6 = s6;
            ss7 = s7;
            sz1 = z1;
            sz2 = z2;
            sz3 = z3;
            sz11 = z11;
            sz12 = z12;
            sz13 = z13;
            sz21 = z21;
            sz22 = z22;
            sz23 = z23;
            sz31 = z31;
            sz32 = z32;
            sz33 = z33;
            zcosg = zcosgl;
            zsing = zsingl;
            zcosi = zcosil;
            zsini = zsinil;
            zcosh = zcoshl * cnodm + zsinhl * snodm;
            zsinh = snodm * zcoshl - cnodm * zsinhl;
            cc = c1l;
        }
    }

    zmol = (4.719_967_2 + 0.229_971_50 * day - gam) % TWOPI;
    zmos = (6.256_583_7 + 0.017_201_977 * day) % TWOPI;

    // solar terms
    se2 = 2.0 * ss1 * ss6;
    se3 = 2.0 * ss1 * ss7;
    si2 = 2.0 * ss2 * sz12;
    si3 = 2.0 * ss2 * (sz13 - sz11);
    sl2 = -2.0 * ss3 * sz2;
    sl3 = -2.0 * ss3 * (sz3 - sz1);
    sl4 = -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
    sgh2 = 2.0 * ss4 * sz32;
    sgh3 = 2.0 * ss4 * (sz33 - sz31);
    sgh4 = -18.0 * ss4 * zes;
    sh2 = -2.0 * ss2 * sz22;
    sh3 = -2.0 * ss2 * (sz23 - sz21);

    // lunar terms
    ee2 = 2.0 * s1 * s6;
    e3 = 2.0 * s1 * s7;
    xi2 = 2.0 * s2 * z12;
    xi3 = 2.0 * s2 * (z13 - z11);
    xl2 = -2.0 * s3 * z2;
    xl3 = -2.0 * s3 * (z3 - z1);
    xl4 = -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
    xgh2 = 2.0 * s4 * z32;
    xgh3 = 2.0 * s4 * (z33 - z31);
    xgh4 = -18.0 * s4 * zel;
    xh2 = -2.0 * s2 * z22;
    xh3 = -2.0 * s2 * (z23 - z21);

    return (
        snodm, cnodm, sinim, cosim, sinomm, cosomm, day, e3, ee2, em, emsq, gam, peo, pgho, pho,
        pinco, plo, rtemsq, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, s1, s2,
        s3, s4, s5, s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3, sz11, sz12, sz13,
        sz21, sz22, sz23, sz31, sz32, sz33, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4,
        nm, z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33, zmol, zmos,
    );
}

// -----------------------------------------------------------------------------
// dsinit
// -----------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn dsinit(
    xke: f64,
    cosim: f64,
    emsq: f64,
    argpo: f64,
    s1: f64,
    s2: f64,
    s3: f64,
    s4: f64,
    s5: f64,
    sinim: f64,
    ss1: f64,
    ss2: f64,
    ss3: f64,
    ss4: f64,
    ss5: f64,
    sz1: f64,
    sz3: f64,
    sz11: f64,
    sz13: f64,
    sz21: f64,
    sz23: f64,
    sz31: f64,
    sz33: f64,
    t: f64,
    tc: f64,
    gsto: f64,
    mo: f64,
    mdot: f64,
    no: f64,
    nodeo: f64,
    nodedot: f64,
    xpidot: f64,
    z1: f64,
    z3: f64,
    z11: f64,
    z13: f64,
    z21: f64,
    z23: f64,
    z31: f64,
    z33: f64,
    ecco: f64,
    eccsq: f64,
    mut em: f64,
    mut argpm: f64,
    mut inclm: f64,
    mut mm: f64,
    mut nm: f64,
    mut nodem: f64,
    mut irez: i32,
    mut atime: f64,
    mut d2201: f64,
    mut d2211: f64,
    mut d3210: f64,
    mut d3222: f64,
    mut d4410: f64,
    mut d4422: f64,
    mut d5220: f64,
    mut d5232: f64,
    mut d5421: f64,
    mut d5433: f64,
    mut dedt: f64,
    mut didt: f64,
    mut dmdt: f64,
    mut dnodt: f64,
    mut domdt: f64,
    mut del1: f64,
    mut del2: f64,
    mut del3: f64,
    mut xfact: f64,
    mut xlamo: f64,
    mut xli: f64,
    mut xni: f64,
) -> (
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    i32,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let q22 = 1.7891679e-6_f64;
    let q31 = 2.1460748e-6_f64;
    let q33 = 2.2123015e-7_f64;
    let root22 = 1.7891679e-6_f64;
    let root44 = 7.3636953e-9_f64;
    let root54 = 2.1765803e-9_f64;
    let rptim = 4.37526908801129966e-3_f64;
    let root32 = 3.7393792e-7_f64;
    let root52 = 1.1428639e-7_f64;
    let x2o3 = 2.0 / 3.0;
    let znl = 1.5835218e-4_f64;
    let zns = 1.19459e-5_f64;

    irez = 0;
    if nm > 0.0034906585 && nm < 0.0052359877 {
        irez = 1;
    }
    if nm >= 8.26e-3 && nm <= 9.24e-3 && em >= 0.5 {
        irez = 2;
    }

    let ses = ss1 * zns * ss5;
    let sis = ss2 * zns * (sz11 + sz13);
    let sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
    let sghs = ss4 * zns * (sz31 + sz33 - 6.0);
    let mut shs = -zns * ss2 * (sz21 + sz23);
    if !(inclm < 5.2359877e-2 || inclm > PI - 5.2359877e-2) && sinim != 0.0 {
        // leave shs as computed
    } else {
        shs = 0.0;
    }
    if sinim != 0.0 {
        shs /= sinim;
    }
    let sgs = sghs - cosim * shs;

    dedt = ses + s1 * znl * s5;
    didt = sis + s2 * znl * (z11 + z13);
    dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
    let sghl = s4 * znl * (z31 + z33 - 6.0);
    let mut shll = -znl * s2 * (z21 + z23);
    if inclm < 5.2359877e-2 || inclm > PI - 5.2359877e-2 {
        shll = 0.0;
    }
    domdt = sgs + sghl;
    dnodt = shs;
    if sinim != 0.0 {
        domdt -= cosim / sinim * shll;
        dnodt += shll / sinim;
    }

    let mut dndt = 0.0;
    let theta = (gsto + tc * rptim) % TWOPI;
    em += dedt * t;
    inclm += didt * t;
    argpm += domdt * t;
    nodem += dnodt * t;
    mm += dmdt * t;

    if irez != 0 {
        let aonv = (nm / xke).powf(x2o3);

        if irez == 2 {
            let cosisq = cosim * cosim;
            let emo = em;
            em = ecco;
            let emsqo = emsq;
            let mut emsq_local = eccsq;
            let eoc = em * emsq_local;
            let mut g201 = -0.306 - (em - 0.64) * 0.440;

            let (g211, g310, g322, g410, g422, g520);
            if em <= 0.65 {
                g211 = 3.616 - 13.2470 * em + 16.2900 * emsq_local;
                g310 = -19.302 + 117.3900 * em - 228.4190 * emsq_local + 156.5910 * eoc;
                g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq_local + 146.5816 * eoc;
                g410 = -41.122 + 242.6940 * em - 471.0940 * emsq_local + 313.9530 * eoc;
                g422 = -146.407 + 841.8800 * em - 1629.014 * emsq_local + 1083.4350 * eoc;
                g520 = -532.114 + 3017.977 * em - 5740.032 * emsq_local + 3708.2760 * eoc;
            } else {
                g211 = -72.099 + 331.819 * em - 508.738 * emsq_local + 266.724 * eoc;
                g310 = -346.844 + 1582.851 * em - 2415.925 * emsq_local + 1246.113 * eoc;
                g322 = -342.585 + 1554.908 * em - 2366.899 * emsq_local + 1215.972 * eoc;
                g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq_local + 3651.957 * eoc;
                g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq_local + 12422.520 * eoc;
                if em > 0.715 {
                    g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq_local + 31324.56 * eoc;
                } else {
                    g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq_local;
                }
            }

            let (g533, g521, g532);
            if em < 0.7 {
                g533 = -919.22770 + 4988.61 * em - 9064.77 * emsq_local + 5542.21 * eoc;
                g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq_local + 5337.524 * eoc;
                g532 = -853.66600 + 4690.25 * em - 8624.77 * emsq_local + 5341.4 * eoc;
            } else {
                g533 = -37995.78 + 161616.52 * em - 229838.2 * emsq_local + 109377.94 * eoc;
                g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq_local + 146349.42 * eoc;
                g532 = -40023.88 + 170470.89 * em - 242699.48 * emsq_local + 115605.82 * eoc;
            }

            let sini2 = sinim * sinim;
            let f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
            let f221 = 1.5 * sini2;
            let f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
            let f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
            let f441 = 35.0 * sini2 * f220;
            let f442 = 39.375 * sini2 * sini2;
            let f522 = 9.84375
                * sinim
                * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq)
                    + 1.0 / 3.0 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
            let f523 = sinim
                * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq)
                    + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
            let f542 = 29.53125
                * sinim
                * (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
            let f543 = 29.53125
                * sinim
                * (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));

            let xno2 = nm * nm;
            let ainv2 = aonv * aonv;
            let mut temp1 = 3.0 * xno2 * ainv2;
            let mut temp = temp1 * root22;
            d2201 = temp * f220 * g201;
            d2211 = temp * f221 * g211;
            temp1 *= aonv;
            temp = temp1 * root32;
            d3210 = temp * f321 * g310;
            d3222 = temp * f322 * g322;
            temp1 *= aonv;
            temp = 2.0 * temp1 * root44;
            d4410 = temp * f441 * g410;
            d4422 = temp * f442 * g422;
            temp1 *= aonv;
            temp = temp1 * root52;
            d5220 = temp * f522 * g520;
            d5232 = temp * f523 * g532;
            temp = 2.0 * temp1 * root54;
            d5421 = temp * f542 * g521;
            d5433 = temp * f543 * g533;
            xlamo = (mo + nodeo + nodeo - theta - theta) % TWOPI;
            xfact = mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
            em = emo;
            emsq_local = emsqo;
        }

        if irez == 1 {
            let g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
            let g310 = 1.0 + 2.0 * emsq;
            let g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
            let f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
            let f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
            let mut f330 = 1.0 + cosim;
            f330 = 1.875 * f330 * f330 * f330;
            del1 = 3.0 * nm * nm * aonv * aonv;
            del2 = 2.0 * del1 * f220 * g200 * q22;
            del3 = 3.0 * del1 * f330 * g300 * q33 * aonv;
            del1 = del1 * f311 * g310 * q31 * aonv;
            xlamo = (mo + nodeo + argpo - theta) % TWOPI;
            xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no;
        }

        xli = xlamo;
        xni = no;
        atime = 0.0;
        nm = no + dndt;
    }

    (
        em, argpm, inclm, mm, nm, nodem, irez, atime, d2201, d2211, d3210, d3222, d4410, d4422,
        d5220, d5232, d5421, d5433, dedt, didt, dmdt, dndt, dnodt, domdt, del1, del2, del3, xfact,
        xlamo, xli, xni,
    )
}

// -----------------------------------------------------------------------------
// dspace
// -----------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn dspace(
    irez: i32,
    d2201: f64,
    d2211: f64,
    d3210: f64,
    d3222: f64,
    d4410: f64,
    d4422: f64,
    d5220: f64,
    d5232: f64,
    d5421: f64,
    d5433: f64,
    dedt: f64,
    del1: f64,
    del2: f64,
    del3: f64,
    didt: f64,
    dmdt: f64,
    dnodt: f64,
    domdt: f64,
    argpo: f64,
    argpdot: f64,
    t: f64,
    tc: f64,
    gsto: f64,
    xfact: f64,
    xlamo: f64,
    no: f64,
    mut atime: f64,
    mut em: f64,
    mut argpm: f64,
    mut inclm: f64,
    mut xli: f64,
    mut mm: f64,
    mut xni: f64,
    mut nodem: f64,
    mut nm: f64,
) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let fasx2 = 0.13130908_f64;
    let fasx4 = 2.8843198_f64;
    let fasx6 = 0.37448087_f64;
    let g22 = 5.7686396_f64;
    let g32 = 0.95240898_f64;
    let g44 = 1.8014998_f64;
    let g52 = 1.0508330_f64;
    let g54 = 4.4108898_f64;
    let rptim = 4.37526908801129966e-3_f64;
    let stepp = 720.0_f64;
    let stepn = -720.0_f64;
    let step2 = 259200.0_f64;

    let mut dndt = 0.0;
    let theta = (gsto + tc * rptim) % TWOPI;
    em += dedt * t;
    inclm += didt * t;
    argpm += domdt * t;
    nodem += dnodt * t;
    mm += dmdt * t;

    let mut ft = 0.0;
    if irez != 0 {
        if atime == 0.0 || t * atime <= 0.0 || t.abs() < atime.abs() {
            atime = 0.0;
            xni = no;
            xli = xlamo;
        }

        let delt = if t > 0.0 { stepp } else { stepn };

        let mut iretn = 381_i32;

        let mut xndt = 0.0;
        let mut xldot = 0.0;
        let mut xnddt = 0.0;

        while iretn == 381 {
            if irez != 2 {
                xndt = del1 * (xli - fasx2).sin()
                    + del2 * (2.0 * (xli - fasx4)).sin()
                    + del3 * (3.0 * (xli - fasx6)).sin();
                xldot = xni + xfact;
                xnddt = del1 * (xli - fasx2).cos()
                    + 2.0 * del2 * (2.0 * (xli - fasx4)).cos()
                    + 3.0 * del3 * (3.0 * (xli - fasx6)).cos();
                xnddt *= xldot;
            } else {
                let xomi = argpo + argpdot * atime;
                let x2omi = xomi + xomi;
                let x2li = xli + xli;
                xndt = d2201 * (x2omi + xli - g22).sin()
                    + d2211 * (xli - g22).sin()
                    + d3210 * (xomi + xli - g32).sin()
                    + d3222 * (-xomi + xli - g32).sin()
                    + d4410 * (x2omi + x2li - g44).sin()
                    + d4422 * (x2li - g44).sin()
                    + d5220 * (xomi + xli - g52).sin()
                    + d5232 * (-xomi + xli - g52).sin()
                    + d5421 * (xomi + x2li - g54).sin()
                    + d5433 * (-xomi + x2li - g54).sin();
                xldot = xni + xfact;
                xnddt = d2201 * (x2omi + xli - g22).cos()
                    + d2211 * (xli - g22).cos()
                    + d3210 * (xomi + xli - g32).cos()
                    + d3222 * (-xomi + xli - g32).cos()
                    + d5220 * (xomi + xli - g52).cos()
                    + d5232 * (-xomi + xli - g52).cos()
                    + 2.0
                        * (d4410 * (x2omi + x2li - g44).cos()
                            + d4422 * (x2li - g44).cos()
                            + d5421 * (xomi + x2li - g54).cos()
                            + d5433 * (-xomi + x2li - g54).cos());
                xnddt *= xldot;
            }

            if (t - atime).abs() >= stepp {
                iretn = 381;
            } else {
                ft = t - atime;
                iretn = 0;
            }

            if iretn == 381 {
                xli = xli + xldot * delt + xndt * step2;
                xni = xni + xndt * delt + xnddt * step2;
                atime += delt;
            }
        }

        nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
        let xl = xli + xldot * ft + xndt * ft * ft * 0.5;
        if irez != 1 {
            mm = xl - 2.0 * nodem + 2.0 * theta;
            dndt = nm - no;
        } else {
            mm = xl - nodem - argpm + theta;
            dndt = nm - no;
        }
        nm = no + dndt;
    }

    (atime, em, argpm, inclm, xli, mm, xni, nodem, dndt, nm)
}

// -----------------------------------------------------------------------------
// initl
// -----------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn initl(
    xke: f64,
    j2: f64,
    ecco: f64,
    epoch: f64,
    inclo: f64,
    mut no: f64,
    method: &mut char,
    opsmode: char,
) -> (
    f64,
    char,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let x2o3 = 2.0 / 3.0;
    let eccsq = ecco * ecco;
    let omeosq = 1.0 - eccsq;
    let rteosq = omeosq.sqrt();
    let cosio = inclo.cos();
    let cosio2 = cosio * cosio;

    let ak = (xke / no).powf(x2o3);
    let d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    let mut del_ = d1 / (ak * ak);
    let adel = ak * (1.0 - del_ * del_ - del_ * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0));
    del_ = d1 / (adel * adel);
    no = no / (1.0 + del_);

    let ao = (xke / no).powf(x2o3);
    let sinio = inclo.sin();
    let po = ao * omeosq;
    let con42 = 1.0 - 5.0 * cosio2;
    let con41 = -con42 - cosio2 - cosio2;
    let ainv = 1.0 / ao;
    let posq = po * po;
    let rp = ao * (1.0 - ecco);
    *method = 'n';

    let gsto = if opsmode == 'a' {
        let ts70 = epoch - 7305.0;
        let ds70 = (ts70 + 1.0e-8).floor();
        let tfrac = ts70 - ds70;
        let c1 = 1.720_279_169_407_036_39e-2_f64;
        let thgr70 = 1.732_134_385_650_937_4_f64;
        let fk5r = 5.075_514_194_322_694_42e-15_f64;
        let c1p2p = c1 + TWOPI;
        let mut gsto_local = (thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r) % TWOPI;
        if gsto_local < 0.0 {
            gsto_local += TWOPI;
        }
        gsto_local
    } else {
        gstime(epoch + 2_433_281.5)
    };

    (
        no, *method, ainv, ao, con41, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio,
        gsto,
    )
}

// -----------------------------------------------------------------------------
// sgp4init
// -----------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn sgp4init(
    whichconst: &str,
    opsmode: char,
    satn: &str,
    epoch: f64,
    xbstar: f64,
    xndot: f64,
    xnddot: f64,
    xecco: f64,
    xargpo: f64,
    xinclo: f64,
    xmo: f64,
    xno_kozai: f64,
    xnodeo: f64,
    satrec: &mut SatRec,
) -> bool {
    let temp4 = 1.5e-12_f64;

    // Zero near-earth vars
    satrec.isimp = 0;
    satrec.method = 'n';
    satrec.aycof = 0.0;
    satrec.con41 = 0.0;
    satrec.cc1 = 0.0;
    satrec.cc4 = 0.0;
    satrec.cc5 = 0.0;
    satrec.d2 = 0.0;
    satrec.d3 = 0.0;
    satrec.d4 = 0.0;
    satrec.delmo = 0.0;
    satrec.eta = 0.0;
    satrec.argpdot = 0.0;
    satrec.omgcof = 0.0;
    satrec.sinmao = 0.0;
    satrec.t = 0.0;
    satrec.t2cof = 0.0;
    satrec.t3cof = 0.0;
    satrec.t4cof = 0.0;
    satrec.t5cof = 0.0;
    satrec.x1mth2 = 0.0;
    satrec.x7thm1 = 0.0;
    satrec.mdot = 0.0;
    satrec.nodedot = 0.0;
    satrec.xlcof = 0.0;
    satrec.xmcof = 0.0;
    satrec.nodecf = 0.0;

    // Zero deep-space vars
    satrec.irez = 0;
    satrec.d2201 = 0.0;
    satrec.d2211 = 0.0;
    satrec.d3210 = 0.0;
    satrec.d3222 = 0.0;
    satrec.d4410 = 0.0;
    satrec.d4422 = 0.0;
    satrec.d5220 = 0.0;
    satrec.d5232 = 0.0;
    satrec.d5421 = 0.0;
    satrec.d5433 = 0.0;
    satrec.dedt = 0.0;
    satrec.del1 = 0.0;
    satrec.del2 = 0.0;
    satrec.del3 = 0.0;
    satrec.didt = 0.0;
    satrec.dmdt = 0.0;
    satrec.dnodt = 0.0;
    satrec.domdt = 0.0;
    satrec.e3 = 0.0;
    satrec.ee2 = 0.0;
    satrec.peo = 0.0;
    satrec.pgho = 0.0;
    satrec.pho = 0.0;
    satrec.pinco = 0.0;
    satrec.plo = 0.0;
    satrec.se2 = 0.0;
    satrec.se3 = 0.0;
    satrec.sgh2 = 0.0;
    satrec.sgh3 = 0.0;
    satrec.sgh4 = 0.0;
    satrec.sh2 = 0.0;
    satrec.sh3 = 0.0;
    satrec.si2 = 0.0;
    satrec.si3 = 0.0;
    satrec.sl2 = 0.0;
    satrec.sl3 = 0.0;
    satrec.sl4 = 0.0;
    satrec.gsto = 0.0;
    satrec.xfact = 0.0;
    satrec.xgh2 = 0.0;
    satrec.xgh3 = 0.0;
    satrec.xgh4 = 0.0;
    satrec.xh2 = 0.0;
    satrec.xh3 = 0.0;
    satrec.xi2 = 0.0;
    satrec.xi3 = 0.0;
    satrec.xl2 = 0.0;
    satrec.xl3 = 0.0;
    satrec.xl4 = 0.0;
    satrec.xlamo = 0.0;
    satrec.zmol = 0.0;
    satrec.zmos = 0.0;
    satrec.atime = 0.0;
    satrec.xli = 0.0;
    satrec.xni = 0.0;

    // Earth constants
    let gc = get_grav_const(whichconst);
    satrec.tumin = gc.tumin;
    satrec.mu = gc.mu;
    satrec.radiusearthkm = gc.radiusearthkm;
    satrec.xke = gc.xke;
    satrec.j2 = gc.j2;
    satrec.j3 = gc.j3;
    satrec.j4 = gc.j4;
    satrec.j3oj2 = gc.j3oj2;

    satrec.error = 0;
    satrec.operationmode = opsmode;
    satrec.satnum_str = satn.to_string();
    satrec.classification = 'U';

    satrec.bstar = xbstar;
    satrec.ndot = xndot;
    satrec.nddot = xnddot;
    satrec.ecco = xecco;
    satrec.argpo = xargpo;
    satrec.inclo = xinclo;
    satrec.mo = xmo;
    satrec.no_kozai = xno_kozai;
    satrec.nodeo = xnodeo;

    satrec.am = 0.0;
    satrec.em = 0.0;
    satrec.im = 0.0;
    satrec.Om = 0.0;
    satrec.om = 0.0;
    satrec.mm = 0.0;
    satrec.nm = 0.0;

    let ss = 78.0 / satrec.radiusearthkm + 1.0;
    let qzms2ttemp = (120.0 - 78.0) / satrec.radiusearthkm;
    let qzms2t = qzms2ttemp.powi(4);
    let x2o3 = 2.0 / 3.0;

    satrec.init = 'y';
    satrec.t = 0.0;

    let mut method = satrec.method;
    let (
        no_unkozai,
        method_out,
        ainv,
        ao,
        con41,
        con42,
        cosio,
        cosio2,
        eccsq,
        omeosq,
        posq,
        rp,
        rteosq,
        sinio,
        gsto,
    ) = initl(
        satrec.xke,
        satrec.j2,
        satrec.ecco,
        epoch,
        satrec.inclo,
        satrec.no_kozai,
        &mut method,
        satrec.operationmode,
    );
    satrec.no_unkozai = no_unkozai;
    satrec.method = method_out;
    satrec.con41 = con41;
    satrec.gsto = gsto;

    satrec.a = (satrec.no_unkozai * satrec.tumin).powf(-2.0 / 3.0);
    satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;

    if omeosq >= 0.0 || satrec.no_unkozai >= 0.0 {
        satrec.isimp = 0;
        if rp < 220.0 / satrec.radiusearthkm + 1.0 {
            satrec.isimp = 1;
        }
        let mut sfour = ss;
        let mut qzms24 = qzms2t;
        let perige = (rp - 1.0) * satrec.radiusearthkm;

        if perige < 156.0 {
            sfour = perige - 78.0;
            if perige < 98.0 {
                sfour = 20.0;
            }
            let qzms24temp = (120.0 - sfour) / satrec.radiusearthkm;
            qzms24 = qzms24temp.powi(4);
            sfour = sfour / satrec.radiusearthkm + 1.0;
        }

        let pinvsq = 1.0 / posq;
        let tsi = 1.0 / (ao - sfour);
        satrec.eta = ao * satrec.ecco * tsi;
        let etasq = satrec.eta * satrec.eta;
        let eeta = satrec.ecco * satrec.eta;
        let psisq = (1.0 - etasq).abs();
        let coef = qzms24 * tsi.powi(4);
        let coef1 = coef / psisq.powf(3.5);
        let cc2 = coef1
            * satrec.no_unkozai
            * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                + 0.375 * satrec.j2 * tsi / psisq
                    * satrec.con41
                    * (8.0 + 3.0 * etasq * (8.0 + etasq)));
        satrec.cc1 = satrec.bstar * cc2;
        let mut cc3 = 0.0;
        if satrec.ecco > 1.0e-4 {
            cc3 = -2.0 * coef * tsi * satrec.j3oj2 * satrec.no_unkozai * sinio / satrec.ecco;
        }
        satrec.x1mth2 = 1.0 - cosio2;
        satrec.cc4 = 2.0
            * satrec.no_unkozai
            * coef1
            * ao
            * omeosq
            * (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco * (0.5 + 2.0 * etasq)
                - satrec.j2 * tsi / (ao * psisq)
                    * (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                        + 0.75
                            * satrec.x1mth2
                            * (2.0 * etasq - eeta * (1.0 + etasq))
                            * (2.0 * satrec.argpo).cos()));
        satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

        let cosio4 = cosio2 * cosio2;
        let temp1 = 1.5 * satrec.j2 * pinvsq * satrec.no_unkozai;
        let temp2 = 0.5 * temp1 * satrec.j2 * pinvsq;
        let temp3 = -0.46875 * satrec.j4 * pinvsq * pinvsq * satrec.no_unkozai;
        satrec.mdot = satrec.no_unkozai
            + 0.5 * temp1 * rteosq * satrec.con41
            + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
        satrec.argpdot = -0.5 * temp1 * con42
            + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4)
            + temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
        let xhdot1 = -temp1 * cosio;
        satrec.nodedot = xhdot1
            + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
        let xpidot = satrec.argpdot + satrec.nodedot;
        satrec.omgcof = satrec.bstar * cc3 * satrec.argpo.cos();
        satrec.xmcof = if satrec.ecco > 1.0e-4 {
            -x2o3 * coef * satrec.bstar / eeta
        } else {
            0.0
        };
        satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
        satrec.t2cof = 1.5 * satrec.cc1;
        satrec.xlcof = if (cosio + 1.0).abs() > 1.5e-12 {
            -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
        } else {
            -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4
        };
        satrec.aycof = -0.5 * satrec.j3oj2 * sinio;
        let delmotemp = 1.0 + satrec.eta * satrec.mo.cos();
        satrec.delmo = delmotemp * delmotemp * delmotemp;
        satrec.sinmao = satrec.mo.sin();
        satrec.x7thm1 = 7.0 * cosio2 - 1.0;

        if 2.0 * PI / satrec.no_unkozai >= 225.0 {
            satrec.method = 'd';
            satrec.isimp = 1;
            let tc = 0.0;
            let mut inclm = satrec.inclo;

            let (
                snodm,
                cnodm,
                sinim,
                cosim,
                sinomm,
                cosomm,
                day,
                e3,
                ee2,
                em,
                emsq,
                gam,
                peo,
                pgho,
                pho,
                pinco,
                plo,
                rtemsq,
                se2,
                se3,
                sgh2,
                sgh3,
                sgh4,
                sh2,
                sh3,
                si2,
                si3,
                sl2,
                sl3,
                sl4,
                s1,
                s2,
                s3,
                s4,
                s5,
                s6,
                s7,
                ss1,
                ss2,
                ss3,
                ss4,
                ss5,
                ss6,
                ss7,
                sz1,
                sz2,
                sz3,
                sz11,
                sz12,
                sz13,
                sz21,
                sz22,
                sz23,
                sz31,
                sz32,
                sz33,
                xgh2,
                xgh3,
                xgh4,
                xh2,
                xh3,
                xi2,
                xi3,
                xl2,
                xl3,
                xl4,
                mut nm,
                z1,
                z2,
                z3,
                z11,
                z12,
                z13,
                z21,
                z22,
                z23,
                z31,
                z32,
                z33,
                zmol,
                zmos,
            ) = dscom(
                epoch,
                satrec.ecco,
                satrec.argpo,
                tc,
                satrec.inclo,
                satrec.nodeo,
                satrec.no_unkozai,
            );
            satrec.e3 = e3;
            satrec.ee2 = ee2;
            satrec.peo = peo;
            satrec.pgho = pgho;
            satrec.pho = pho;
            satrec.pinco = pinco;
            satrec.plo = plo;
            satrec.se2 = se2;
            satrec.se3 = se3;
            satrec.sgh2 = sgh2;
            satrec.sgh3 = sgh3;
            satrec.sgh4 = sgh4;
            satrec.sh2 = sh2;
            satrec.sh3 = sh3;
            satrec.si2 = si2;
            satrec.si3 = si3;
            satrec.sl2 = sl2;
            satrec.sl3 = sl3;
            satrec.sl4 = sl4;
            satrec.zmol = zmol;
            satrec.zmos = zmos;
            satrec.xgh2 = xgh2;
            satrec.xgh3 = xgh3;
            satrec.xgh4 = xgh4;
            satrec.xh2 = xh2;
            satrec.xh3 = xh3;
            satrec.xi2 = xi2;
            satrec.xi3 = xi3;
            satrec.xl2 = xl2;
            satrec.xl3 = xl3;
            satrec.xl4 = xl4;

            let (ecco, inclo, nodeo, argpo, mo) = dpper(
                &satrec,
                inclm,                // inclo
                satrec.init,          // init: char
                satrec.ecco,          // ep
                satrec.inclo,         // inclp
                satrec.nodeo,         // nodep
                satrec.argpo,         // argpp
                satrec.mo,            // mp
                satrec.operationmode, // opsmode
            );
            satrec.ecco = ecco;
            satrec.inclo = inclo;
            satrec.nodeo = nodeo;
            satrec.argpo = argpo;
            satrec.mo = mo;

            let mut argpm = 0.0;
            let mut nodem = 0.0;
            let mut mm = 0.0;

            let (
                em_out,
                argpm_out,
                inclm_out,
                mm_out,
                nm_out,
                nodem_out,
                irez,
                atime,
                d2201,
                d2211,
                d3210,
                d3222,
                d4410,
                d4422,
                d5220,
                d5232,
                d5421,
                d5433,
                dedt,
                didt,
                dmdt,
                dndt,
                dnodt,
                domdt,
                del1,
                del2,
                del3,
                xfact,
                xlamo,
                xli,
                xni,
            ) = dsinit(
                satrec.xke,
                cosim,
                emsq,
                satrec.argpo,
                s1,
                s2,
                s3,
                s4,
                s5,
                sinim,
                ss1,
                ss2,
                ss3,
                ss4,
                ss5,
                sz1,
                sz3,
                sz11,
                sz13,
                sz21,
                sz23,
                sz31,
                sz33,
                satrec.t,
                tc,
                satrec.gsto,
                satrec.mo,
                satrec.mdot,
                satrec.no_unkozai,
                satrec.nodeo,
                satrec.nodedot,
                xpidot,
                z1,
                z3,
                z11,
                z13,
                z21,
                z23,
                z31,
                z33,
                satrec.ecco,
                eccsq,
                em,
                argpm,
                inclm,
                mm,
                nm,
                nodem,
                satrec.irez,
                satrec.atime,
                satrec.d2201,
                satrec.d2211,
                satrec.d3210,
                satrec.d3222,
                satrec.d4410,
                satrec.d4422,
                satrec.d5220,
                satrec.d5232,
                satrec.d5421,
                satrec.d5433,
                satrec.dedt,
                satrec.didt,
                satrec.dmdt,
                satrec.dnodt,
                satrec.domdt,
                satrec.del1,
                satrec.del2,
                satrec.del3,
                satrec.xfact,
                satrec.xlamo,
                satrec.xli,
                satrec.xni,
            );
            satrec.ecco = em_out;
            satrec.inclo = inclm_out;
            satrec.nodeo = nodem_out;
            satrec.argpo = argpm_out;
            satrec.mo = mm_out;
            satrec.irez = irez;
            satrec.atime = atime;
            satrec.d2201 = d2201;
            satrec.d2211 = d2211;
            satrec.d3210 = d3210;
            satrec.d3222 = d3222;
            satrec.d4410 = d4410;
            satrec.d4422 = d4422;
            satrec.d5220 = d5220;
            satrec.d5232 = d5232;
            satrec.d5421 = d5421;
            satrec.d5433 = d5433;
            satrec.dedt = dedt;
            satrec.didt = didt;
            satrec.dmdt = dmdt;
            satrec.dnodt = dnodt;
            satrec.domdt = domdt;
            satrec.del1 = del1;
            satrec.del2 = del2;
            satrec.del3 = del3;
            satrec.xfact = xfact;
            satrec.xlamo = xlamo;
            satrec.xli = xli;
            satrec.xni = xni;
        }

        if satrec.isimp != 1 {
            let cc1sq = satrec.cc1 * satrec.cc1;
            satrec.d2 = 4.0 * ao * tsi * cc1sq;
            let temp = satrec.d2 * tsi * satrec.cc1 / 3.0;
            satrec.d3 = (17.0 * ao + sfour) * temp;
            satrec.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * satrec.cc1;
            satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
            satrec.t4cof =
                0.25 * (3.0 * satrec.d3 + satrec.cc1 * (12.0 * satrec.d2 + 10.0 * cc1sq));
            satrec.t5cof = 0.2
                * (3.0 * satrec.d4
                    + 12.0 * satrec.cc1 * satrec.d3
                    + 6.0 * satrec.d2 * satrec.d2
                    + 15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
        }
    }

    // Finally, propagate to 0.0 to prime everything
    let _ = sgp4(satrec, 0.0);

    satrec.init = 'n';
    true
}

// -----------------------------------------------------------------------------
// sgp4
// -----------------------------------------------------------------------------

pub fn sgp4(satrec: &mut SatRec, tsince: f64) -> (Option<[f64; 3]>, Option<[f64; 3]>) {
    let temp4 = 1.5e-12_f64;
    let x2o3 = 2.0 / 3.0;
    let vkmpersec = satrec.radiusearthkm * satrec.xke / 60.0;

    satrec.t = tsince;
    satrec.error = 0;
    satrec.error_message = None;

    let xmdf = satrec.mo + satrec.mdot * satrec.t;
    let argpdf = satrec.argpo + satrec.argpdot * satrec.t;
    let nodedf = satrec.nodeo + satrec.nodedot * satrec.t;
    let mut argpm = argpdf;
    let mut mm = xmdf;
    let t2 = satrec.t * satrec.t;
    let mut nodem = nodedf + satrec.nodecf * t2;
    let mut tempa = 1.0 - satrec.cc1 * satrec.t;
    let mut tempe = satrec.bstar * satrec.cc4 * satrec.t;
    let mut templ = satrec.t2cof * t2;

    if satrec.isimp != 1 {
        let delomg = satrec.omgcof * satrec.t;
        let delmtemp = 1.0 + satrec.eta * xmdf.cos();
        let delm = satrec.xmcof * (delmtemp * delmtemp * delmtemp - satrec.delmo);
        let temp = delomg + delm;
        mm = xmdf + temp;
        argpm = argpdf - temp;
        let t3 = t2 * satrec.t;
        let t4 = t3 * satrec.t;
        tempa -= satrec.d2 * t2 + satrec.d3 * t3 + satrec.d4 * t4;
        tempe += satrec.bstar * satrec.cc5 * (mm.sin() - satrec.sinmao);
        templ += satrec.t3cof * t3 + t4 * (satrec.t4cof + satrec.t * satrec.t5cof);
    }

    let mut nm = satrec.no_unkozai;
    let mut em = satrec.ecco;
    let mut inclm = satrec.inclo;

    if satrec.method == 'd' {
        let tc = satrec.t;
        let (atime, em_out, argpm_out, inclm_out, xli, mm_out, xni, nodem_out, dndt, nm_out) =
            dspace(
                satrec.irez,
                satrec.d2201,
                satrec.d2211,
                satrec.d3210,
                satrec.d3222,
                satrec.d4410,
                satrec.d4422,
                satrec.d5220,
                satrec.d5232,
                satrec.d5421,
                satrec.d5433,
                satrec.dedt,
                satrec.del1,
                satrec.del2,
                satrec.del3,
                satrec.didt,
                satrec.dmdt,
                satrec.dnodt,
                satrec.domdt,
                satrec.argpo,
                satrec.argpdot,
                satrec.t,
                tc,
                satrec.gsto,
                satrec.xfact,
                satrec.xlamo,
                satrec.no_unkozai,
                satrec.atime,
                em,
                argpm,
                inclm,
                satrec.xli,
                mm,
                satrec.xni,
                nodem,
                nm,
            );
        satrec.atime = atime;
        em = em_out;
        argpm = argpm_out;
        inclm = inclm_out;
        satrec.xli = xli;
        mm = mm_out;
        satrec.xni = xni;
        nodem = nodem_out;
        nm = nm_out;
    }

    if nm <= 0.0 {
        satrec.error_message = Some(format!("mean motion {} is less than zero", nm));
        satrec.error = 2;
        return (None, None);
    }

    let mut am = (satrec.xke / nm).powf(x2o3) * tempa * tempa;
    nm = satrec.xke / am.powf(1.5);
    em -= tempe;

    if em >= 1.0 || em < -0.001 {
        satrec.error_message = Some(format!(
            "mean eccentricity {} not within range 0.0 <= e < 1.0",
            em
        ));
        satrec.error = 1;
        return (None, None);
    }

    if em < 1.0e-6 {
        em = 1.0e-6;
    }
    mm += satrec.no_unkozai * templ;
    let mut xlm = mm + argpm + nodem;
    let emsq = em * em;

    nodem = if nodem >= 0.0 {
        nodem % TWOPI
    } else {
        -((-nodem) % TWOPI)
    };
    argpm = argpm % TWOPI;
    xlm = xlm % TWOPI;
    mm = (xlm - argpm - nodem) % TWOPI;

    satrec.am = am;
    satrec.em = em;
    satrec.im = inclm;
    satrec.Om = nodem;
    satrec.om = argpm;
    satrec.mm = mm;
    satrec.nm = nm;

    let sinim = inclm.sin();
    let cosim = inclm.cos();

    let mut ep = em;
    let mut xincp = inclm;
    let mut argpp = argpm;
    let mut nodep = nodem;
    let mut mp = mm;
    let mut sinip = sinim;
    let mut cosip = cosim;

    if satrec.method == 'd' {
        let (ep_out, xincp_out, nodep_out, argpp_out, mp_out) = dpper(
            satrec,
            satrec.inclo,
            'n',
            ep,
            xincp,
            nodep,
            argpp,
            mp,
            satrec.operationmode,
        );
        ep = ep_out;
        xincp = xincp_out;
        nodep = nodep_out;
        argpp = argpp_out;
        mp = mp_out;

        if xincp < 0.0 {
            xincp = -xincp;
            nodep += PI;
            argpp -= PI;
        }

        if ep < 0.0 || ep > 1.0 {
            satrec.error_message = Some(format!(
                "perturbed eccentricity {} not within range 0.0 <= e <= 1.0",
                ep
            ));
            satrec.error = 3;
            return (None, None);
        }
    }

    if satrec.method == 'd' {
        sinip = xincp.sin();
        cosip = xincp.cos();
        satrec.aycof = -0.5 * satrec.j3oj2 * sinip;
        satrec.xlcof = if (cosip + 1.0).abs() > 1.5e-12 {
            -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip)
        } else {
            -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4
        };
    }

    let axnl = ep * argpp.cos();
    let temp = 1.0 / (am * (1.0 - ep * ep));
    let aynl = ep * argpp.sin() + temp * satrec.aycof;
    let xl = mp + argpp + nodep + temp * satrec.xlcof * axnl;

    let mut u = (xl - nodep) % TWOPI;
    let mut eo1 = u;
    let mut tem5 = 9999.9_f64;
    let mut ktr = 1;

    while tem5.abs() >= 1.0e-12 && ktr <= 10 {
        let sineo1 = eo1.sin();
        let coseo1 = eo1.cos();
        tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if tem5.abs() >= 0.95 {
            tem5 = if tem5 > 0.0 { 0.95 } else { -0.95 };
        }
        eo1 += tem5;
        ktr += 1;
    }

    let sineo1 = eo1.sin();
    let coseo1 = eo1.cos();
    let ecose = axnl * coseo1 + aynl * sineo1;
    let esine = axnl * sineo1 - aynl * coseo1;
    let el2 = axnl * axnl + aynl * aynl;
    let pl = am * (1.0 - el2);

    if pl < 0.0 {
        satrec.error_message = Some(format!("semilatus rectum {} is less than zero", pl));
        satrec.error = 4;
        return (None, None);
    }

    let rl = am * (1.0 - ecose);
    let rdotl = am.sqrt() * esine / rl;
    let rvdotl = pl.sqrt() / rl;
    let betal = (1.0 - el2).sqrt();
    let temp = esine / (1.0 + betal);
    let sinu = am / rl * (sineo1 - aynl - axnl * temp);
    let cosu = am / rl * (coseo1 - axnl + aynl * temp);
    let mut su = sinu.atan2(cosu);
    let sin2u = (cosu + cosu) * sinu;
    let cos2u = 1.0 - 2.0 * sinu * sinu;
    let temp = 1.0 / pl;
    let temp1 = 0.5 * satrec.j2 * temp;
    let temp2 = temp1 * temp;

    if satrec.method == 'd' {
        let cosisq = cosip * cosip;
        satrec.con41 = 3.0 * cosisq - 1.0;
        satrec.x1mth2 = 1.0 - cosisq;
        satrec.x7thm1 = 7.0 * cosisq - 1.0;
    }

    let mrt = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) + 0.5 * temp1 * satrec.x1mth2 * cos2u;
    su -= 0.25 * temp2 * satrec.x7thm1 * sin2u;
    let xnode = nodep + 1.5 * temp2 * cosip * sin2u;
    let xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
    let mvt = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / satrec.xke;
    let rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u + 1.5 * satrec.con41) / satrec.xke;

    let sinsu = su.sin();
    let cossu = su.cos();
    let snod = xnode.sin();
    let cnod = xnode.cos();
    let sini = xinc.sin();
    let cosi = xinc.cos();
    let xmx = -snod * cosi;
    let xmy = cnod * cosi;
    let ux = xmx * sinsu + cnod * cossu;
    let uy = xmy * sinsu + snod * cossu;
    let uz = sini * sinsu;
    let vx = xmx * cossu - cnod * sinsu;
    let vy = xmy * cossu - snod * sinsu;
    let vz = sini * cossu;

    let mr = mrt * satrec.radiusearthkm;
    let r = [mr * ux, mr * uy, mr * uz];
    let v = [
        (mvt * ux + rvdot * vx) * vkmpersec,
        (mvt * uy + rvdot * vy) * vkmpersec,
        (mvt * uz + rvdot * vz) * vkmpersec,
    ];

    if mrt < 1.0 {
        satrec.error_message = Some(format!(
            "mrt {} is less than 1.0 indicating the satellite has decayed",
            mrt
        ));
        satrec.error = 6;
    }

    (Some(r), Some(v))
}
