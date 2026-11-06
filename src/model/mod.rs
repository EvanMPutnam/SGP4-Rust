//! Vallado's SGP4/SDP4 equations.
//!
//! This follows the 2020 revision of Vallado's C++ code and Brandon Rhodes's
//! python-sgp4 port line for line where it can, so each step can be checked
//! against the reference. Variable names keep Vallado's spelling for the
//! same reason. Distances are in Earth radii and times in minutes.
//! Constants are copied digit for digit from the reference, so some carry
//! more precision than an `f64` holds.
#![allow(clippy::excessive_precision)]

mod deep_space;
mod propagate;

use std::f64::consts::TAU;

use crate::elements::Elements;
use crate::gravity::GravityConstants;
use crate::satellite::OpsMode;
use crate::time::JulianDate;
use deep_space::DeepSpace;

const X2O3: f64 = 2.0 / 3.0;
/// Guard against dividing by zero for orbits inclined 180 degrees.
const TEMP4: f64 = 1.5e-12;
/// Julian date of 1950-01-00 00:00, the zero point of Vallado's `epoch`.
const JD_1950: f64 = 2_433_281.5;

/// Everything SGP4 precomputes from one element set.
#[derive(Debug, Clone)]
pub(crate) struct Model {
    pub(crate) gravity: GravityConstants,
    mode: OpsMode,

    // Mean elements at epoch.
    ecco: f64,
    argpo: f64,
    inclo: f64,
    mo: f64,
    nodeo: f64,
    bstar: f64,
    /// Brouwer mean motion, radians per minute.
    pub(crate) no_unkozai: f64,

    // Near-earth secular, drag, and long-period coefficients.
    /// Perigee below 220 km: drop the higher-order drag terms.
    isimp: bool,
    aycof: f64,
    con41: f64,
    cc1: f64,
    cc4: f64,
    cc5: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    delmo: f64,
    eta: f64,
    argpdot: f64,
    omgcof: f64,
    sinmao: f64,
    t2cof: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,
    x1mth2: f64,
    x7thm1: f64,
    mdot: f64,
    nodedot: f64,
    xlcof: f64,
    xmcof: f64,
    nodecf: f64,

    /// Semi-major axis, Earth radii.
    pub(crate) a: f64,
    /// Apogee altitude, Earth radii.
    pub(crate) alta: f64,
    /// Perigee altitude, Earth radii.
    pub(crate) altp: f64,

    /// Present when the period is 225 minutes or longer.
    deep: Option<DeepSpace>,
}

impl Model {
    /// Vallado's `sgp4init` (including `initl`), minus the final call to
    /// `sgp4` at epoch, which the caller makes.
    pub(crate) fn new(elements: &Elements, gravity: GravityConstants, mode: OpsMode) -> Model {
        let GravityConstants {
            radius_km,
            xke,
            j2,
            j3oj2,
            j4,
            tumin,
            ..
        } = gravity;
        let ecco = elements.eccentricity;
        let inclo = elements.inclination;
        let argpo = elements.arg_of_perigee;
        let mo = elements.mean_anomaly;
        let bstar = elements.bstar;
        // Rounded the way Vallado's code rounds it. The deep-space terms are
        // sensitive enough to this that the more precise
        // `(day - JD_1950) + fraction` moves some test vectors by 4e-6 km.
        let epoch = elements.epoch.value() - JD_1950;

        let ss = 78.0 / radius_km + 1.0;
        let qzms2t = ((120.0 - 78.0) / radius_km).powi(4);

        // ---- initl: un-Kozai the mean motion ----
        let eccsq = ecco * ecco;
        let omeosq = 1.0 - eccsq;
        let rteosq = omeosq.sqrt();
        let cosio = inclo.cos();
        let cosio2 = cosio * cosio;

        let ak = (xke / elements.mean_motion).powf(X2O3);
        let d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
        let mut del = d1 / (ak * ak);
        let adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
        del = d1 / (adel * adel);
        let no_unkozai = elements.mean_motion / (1.0 + del);

        let ao = (xke / no_unkozai).powf(X2O3);
        let sinio = inclo.sin();
        let po = ao * omeosq;
        let con42 = 1.0 - 5.0 * cosio2;
        let con41 = -con42 - cosio2 - cosio2;
        let posq = po * po;
        let rp = ao * (1.0 - ecco);
        let gsto = match mode {
            OpsMode::Afspc => afspc_gsto(epoch),
            OpsMode::Improved => JulianDate::new(epoch + JD_1950, 0.0).gmst(),
        };

        let a = (no_unkozai * tumin).powf(-X2O3);

        // ---- sgp4init ----
        let isimp = rp < 220.0 / radius_km + 1.0;
        let mut sfour = ss;
        let mut qzms24 = qzms2t;
        let perige = (rp - 1.0) * radius_km;
        // Lower the drag reference altitude for perigees under 156 km.
        if perige < 156.0 {
            sfour = if perige < 98.0 { 20.0 } else { perige - 78.0 };
            qzms24 = ((120.0 - sfour) / radius_km).powi(4);
            sfour = sfour / radius_km + 1.0;
        }

        let pinvsq = 1.0 / posq;
        let tsi = 1.0 / (ao - sfour);
        let eta = ao * ecco * tsi;
        let etasq = eta * eta;
        let eeta = ecco * eta;
        let psisq = (1.0 - etasq).abs();
        let coef = qzms24 * tsi.powi(4);
        let coef1 = coef / psisq.powf(3.5);
        let cc2 = coef1
            * no_unkozai
            * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                + 0.375 * j2 * tsi / psisq * con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
        let cc1 = bstar * cc2;
        let cc3 = if ecco > 1.0e-4 {
            -2.0 * coef * tsi * j3oj2 * no_unkozai * sinio / ecco
        } else {
            0.0
        };
        let x1mth2 = 1.0 - cosio2;
        let cc4 = 2.0
            * no_unkozai
            * coef1
            * ao
            * omeosq
            * (eta * (2.0 + 0.5 * etasq) + ecco * (0.5 + 2.0 * etasq)
                - j2 * tsi / (ao * psisq)
                    * (-3.0 * con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                        + 0.75
                            * x1mth2
                            * (2.0 * etasq - eeta * (1.0 + etasq))
                            * (2.0 * argpo).cos()));
        let cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

        let cosio4 = cosio2 * cosio2;
        let temp1 = 1.5 * j2 * pinvsq * no_unkozai;
        let temp2 = 0.5 * temp1 * j2 * pinvsq;
        let temp3 = -0.46875 * j4 * pinvsq * pinvsq * no_unkozai;
        let mdot = no_unkozai
            + 0.5 * temp1 * rteosq * con41
            + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
        let argpdot = -0.5 * temp1 * con42
            + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4)
            + temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
        let xhdot1 = -temp1 * cosio;
        let nodedot = xhdot1
            + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
        let xpidot = argpdot + nodedot;
        let delmotemp = 1.0 + eta * mo.cos();

        let mut model = Model {
            gravity,
            mode,
            ecco,
            argpo,
            inclo,
            mo,
            nodeo: elements.raan,
            bstar,
            no_unkozai,
            isimp,
            aycof: -0.5 * j3oj2 * sinio,
            con41,
            cc1,
            cc4,
            cc5,
            d2: 0.0,
            d3: 0.0,
            d4: 0.0,
            delmo: delmotemp * delmotemp * delmotemp,
            eta,
            argpdot,
            omgcof: bstar * cc3 * argpo.cos(),
            sinmao: mo.sin(),
            t2cof: 1.5 * cc1,
            t3cof: 0.0,
            t4cof: 0.0,
            t5cof: 0.0,
            x1mth2,
            x7thm1: 7.0 * cosio2 - 1.0,
            mdot,
            nodedot,
            xlcof: xlcof(j3oj2, sinio, cosio),
            xmcof: if ecco > 1.0e-4 {
                -X2O3 * coef * bstar / eeta
            } else {
                0.0
            },
            nodecf: 3.5 * omeosq * xhdot1 * cc1,
            a,
            alta: a * (1.0 + ecco) - 1.0,
            altp: a * (1.0 - ecco) - 1.0,
            deep: None,
        };

        if TAU / no_unkozai >= 225.0 {
            model.isimp = true;
            model.deep = Some(DeepSpace::new(
                &model,
                deep_space::InitInputs {
                    epoch,
                    gsto,
                    xpidot,
                    eccsq,
                },
            ));
        }

        if !model.isimp {
            let cc1sq = cc1 * cc1;
            let d2 = 4.0 * ao * tsi * cc1sq;
            let temp = d2 * tsi * cc1 / 3.0;
            let d3 = (17.0 * ao + sfour) * temp;
            let d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * cc1;
            model.d2 = d2;
            model.d3 = d3;
            model.d4 = d4;
            model.t3cof = d2 + 2.0 * cc1sq;
            model.t4cof = 0.25 * (3.0 * d3 + cc1 * (12.0 * d2 + 10.0 * cc1sq));
            model.t5cof = 0.2
                * (3.0 * d4 + 12.0 * cc1 * d3 + 6.0 * d2 * d2 + 15.0 * cc1sq * (2.0 * d2 + cc1sq));
        }

        model
    }
}

/// Long-period periodic coefficient for mean longitude.
fn xlcof(j3oj2: f64, sinio: f64, cosio: f64) -> f64 {
    let denominator = if (cosio + 1.0).abs() > TEMP4 {
        1.0 + cosio
    } else {
        TEMP4
    };
    -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / denominator
}

/// Greenwich sidereal time at epoch the way AFSPC code computed it.
fn afspc_gsto(epoch: f64) -> f64 {
    let ts70 = epoch - 7305.0;
    let ds70 = (ts70 + 1.0e-8).floor();
    let tfrac = ts70 - ds70;
    let c1 = 1.720_279_169_407_036_39e-2;
    let thgr70 = 1.732_134_385_650_937_4;
    let fk5r = 5.075_514_194_322_694_42e-15;
    let c1p2p = c1 + TAU;
    let gsto = (thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r) % TAU;
    if gsto < 0.0 { gsto + TAU } else { gsto }
}

/// Reduces an angle with C's `fmod` sign convention, as Vallado does for
/// the node, so `-x` maps to `-(x mod 2π)`.
fn fmod_tau(x: f64) -> f64 {
    if x >= 0.0 { x % TAU } else { -((-x) % TAU) }
}

/// Mean elements partway through one propagation step.
#[derive(Debug, Clone, Copy)]
struct MeanElements {
    em: f64,
    argpm: f64,
    inclm: f64,
    mm: f64,
    nodem: f64,
    nm: f64,
}
