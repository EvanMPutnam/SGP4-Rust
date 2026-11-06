//! SDP4: lunar-solar perturbations and geopotential resonance for orbits
//! with periods of 225 minutes or more. Ports Vallado's `dscom`, `dsinit`,
//! `dpper`, and `dspace`.

use std::f64::consts::{PI, TAU};

use super::{MeanElements, Model, fmod_tau};
use crate::satellite::OpsMode;

const ZNS: f64 = 1.19459e-5;
const ZES: f64 = 0.01675;
const ZNL: f64 = 1.5835218e-4;
const ZEL: f64 = 0.05490;
/// Earth rotation rate, radians per minute.
const RPTIM: f64 = 4.375_269_088_011_299_66e-3;
/// Inclinations within this many radians of 0 or π count as equatorial.
const EQUATORIAL: f64 = 5.2359877e-2;

/// Values `sgp4init` has on hand that the deep-space setup also needs.
pub(super) struct InitInputs {
    /// Days since 1950-01-00 00:00 UTC.
    pub epoch: f64,
    pub gsto: f64,
    pub xpidot: f64,
    pub eccsq: f64,
}

#[derive(Debug, Clone)]
pub(super) struct DeepSpace {
    lunar_solar: LunarSolar,
    /// Greenwich sidereal time at epoch, radians.
    gsto: f64,
    // Secular rates from lunar-solar gravity.
    dedt: f64,
    didt: f64,
    dmdt: f64,
    dnodt: f64,
    domdt: f64,
    resonance: Option<Resonance>,
}

/// Lunar and solar periodic coefficients, computed once by `dscom`.
#[derive(Debug, Clone)]
struct LunarSolar {
    e3: f64,
    ee2: f64,
    se2: f64,
    se3: f64,
    sgh2: f64,
    sgh3: f64,
    sgh4: f64,
    sh2: f64,
    sh3: f64,
    si2: f64,
    si3: f64,
    sl2: f64,
    sl3: f64,
    sl4: f64,
    xgh2: f64,
    xgh3: f64,
    xgh4: f64,
    xh2: f64,
    xh3: f64,
    xi2: f64,
    xi3: f64,
    xl2: f64,
    xl3: f64,
    xl4: f64,
    zmol: f64,
    zmos: f64,
}

/// Geopotential resonance, integrated numerically in `dspace`.
#[derive(Debug, Clone)]
struct Resonance {
    terms: ResonanceTerms,
    xfact: f64,
    xlamo: f64,
}

#[derive(Debug, Clone)]
enum ResonanceTerms {
    /// One revolution per day, as for geostationary orbits (Vallado's `irez == 1`).
    Synchronous { del1: f64, del2: f64, del3: f64 },
    /// Two revolutions per day with eccentricity of 0.5 or more, as for
    /// Molniya orbits (`irez == 2`).
    HalfDay {
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
    },
}

/// Outputs of `dscom` that `dsinit` reads, besides the lunar-solar terms.
struct Common {
    sinim: f64,
    cosim: f64,
    emsq: f64,
    s1: f64,
    s2: f64,
    s3: f64,
    s4: f64,
    s5: f64,
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
    z1: f64,
    z3: f64,
    z11: f64,
    z13: f64,
    z21: f64,
    z23: f64,
    z31: f64,
    z33: f64,
}

impl DeepSpace {
    /// `dscom` followed by `dsinit`, both evaluated at epoch.
    ///
    /// Vallado also calls `dpper` here, but with `init == 'y'` it returns its
    /// inputs unchanged, so it is left out. Its `peo`/`pinco`/`plo`/`pgho`/`pho`
    /// offsets are always zero for the same reason.
    pub(super) fn new(model: &Model, inputs: InitInputs) -> DeepSpace {
        let (common, lunar_solar) = dscom(
            inputs.epoch,
            model.ecco,
            model.argpo,
            model.inclo,
            model.nodeo,
            model.no_unkozai,
        );
        let Common {
            sinim,
            cosim,
            emsq,
            s1,
            s2,
            s3,
            s4,
            s5,
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
            z1,
            z3,
            z11,
            z13,
            z21,
            z23,
            z31,
            z33,
        } = common;

        // ---- dsinit ----
        let nm = model.no_unkozai;
        let em = model.ecco;
        let inclm = model.inclo;
        let equatorial = !(EQUATORIAL..=PI - EQUATORIAL).contains(&inclm);

        // Solar secular terms.
        let ses = ss1 * ZNS * ss5;
        let sis = ss2 * ZNS * (sz11 + sz13);
        let sls = -ZNS * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
        let sghs = ss4 * ZNS * (sz31 + sz33 - 6.0);
        let mut shs = if equatorial {
            0.0
        } else {
            -ZNS * ss2 * (sz21 + sz23)
        };
        if sinim != 0.0 {
            shs /= sinim;
        }
        let sgs = sghs - cosim * shs;

        // Lunar secular terms.
        let dedt = ses + s1 * ZNL * s5;
        let didt = sis + s2 * ZNL * (z11 + z13);
        let dmdt = sls - ZNL * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
        let sghl = s4 * ZNL * (z31 + z33 - 6.0);
        let shll = if equatorial {
            0.0
        } else {
            -ZNL * s2 * (z21 + z23)
        };
        let mut domdt = sgs + sghl;
        let mut dnodt = shs;
        if sinim != 0.0 {
            domdt -= cosim / sinim * shll;
            dnodt += shll / sinim;
        }

        let theta = inputs.gsto % TAU;
        let synchronous = nm > 0.003_490_658_5 && nm < 0.005_235_987_7;
        let half_day = (8.26e-3..=9.24e-3).contains(&nm) && em >= 0.5;

        let resonance = if half_day {
            let aonv = (nm / model.gravity.xke).powf(super::X2O3);
            let terms = half_day_terms(model.ecco, inputs.eccsq, sinim, cosim, nm, aonv);
            Some(Resonance {
                terms,
                xlamo: (model.mo + model.nodeo + model.nodeo - theta - theta) % TAU,
                xfact: model.mdot + dmdt + 2.0 * (model.nodedot + dnodt - RPTIM) - model.no_unkozai,
            })
        } else if synchronous {
            let aonv = (nm / model.gravity.xke).powf(super::X2O3);
            const Q22: f64 = 1.7891679e-6;
            const Q31: f64 = 2.1460748e-6;
            const Q33: f64 = 2.2123015e-7;
            let g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
            let g310 = 1.0 + 2.0 * emsq;
            let g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
            let f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
            let f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
            let f330 = 1.875 * (1.0 + cosim) * (1.0 + cosim) * (1.0 + cosim);
            let del1 = 3.0 * nm * nm * aonv * aonv;
            let del2 = 2.0 * del1 * f220 * g200 * Q22;
            let del3 = 3.0 * del1 * f330 * g300 * Q33 * aonv;
            let del1 = del1 * f311 * g310 * Q31 * aonv;
            Some(Resonance {
                terms: ResonanceTerms::Synchronous { del1, del2, del3 },
                xlamo: (model.mo + model.nodeo + model.argpo - theta) % TAU,
                xfact: model.mdot + inputs.xpidot - RPTIM + dmdt + domdt + dnodt - model.no_unkozai,
            })
        } else {
            None
        };

        DeepSpace {
            lunar_solar,
            gsto: inputs.gsto,
            dedt,
            didt,
            dmdt,
            dnodt,
            domdt,
            resonance,
        }
    }

    /// `dspace`: applies lunar-solar secular rates, then integrates the
    /// resonance terms from epoch to `t` minutes in 720-minute steps.
    ///
    /// Vallado caches the integrator between calls; this restarts it from
    /// epoch each time, which takes the same steps and gives identical results.
    pub(super) fn secular(&self, t: f64, mean: &mut MeanElements, model: &Model) {
        const FASX2: f64 = 0.13130908;
        const FASX4: f64 = 2.8843198;
        const FASX6: f64 = 0.37448087;
        const G22: f64 = 5.7686396;
        const G32: f64 = 0.95240898;
        const G44: f64 = 1.8014998;
        const G52: f64 = 1.0508330;
        const G54: f64 = 4.4108898;
        const STEPP: f64 = 720.0;
        const STEP2: f64 = 259_200.0;

        mean.em += self.dedt * t;
        mean.inclm += self.didt * t;
        mean.argpm += self.domdt * t;
        mean.nodem += self.dnodt * t;
        mean.mm += self.dmdt * t;

        let Some(resonance) = &self.resonance else {
            return;
        };
        let theta = (self.gsto + t * RPTIM) % TAU;
        let no = model.no_unkozai;
        let delt = if t > 0.0 { STEPP } else { -STEPP };

        let mut atime = 0.0;
        let mut xni = no;
        let mut xli = resonance.xlamo;
        // The derivatives at the last step, as (xndt, xldot, xnddt).
        let derivatives = |xli: f64, xni: f64, atime: f64| -> (f64, f64, f64) {
            let xldot = xni + resonance.xfact;
            match resonance.terms {
                ResonanceTerms::Synchronous { del1, del2, del3 } => {
                    let xndt = del1 * (xli - FASX2).sin()
                        + del2 * (2.0 * (xli - FASX4)).sin()
                        + del3 * (3.0 * (xli - FASX6)).sin();
                    let xnddt = del1 * (xli - FASX2).cos()
                        + 2.0 * del2 * (2.0 * (xli - FASX4)).cos()
                        + 3.0 * del3 * (3.0 * (xli - FASX6)).cos();
                    (xndt, xldot, xnddt * xldot)
                }
                ResonanceTerms::HalfDay {
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
                } => {
                    let xomi = model.argpo + model.argpdot * atime;
                    let x2omi = xomi + xomi;
                    let x2li = xli + xli;
                    let xndt = d2201 * (x2omi + xli - G22).sin()
                        + d2211 * (xli - G22).sin()
                        + d3210 * (xomi + xli - G32).sin()
                        + d3222 * (-xomi + xli - G32).sin()
                        + d4410 * (x2omi + x2li - G44).sin()
                        + d4422 * (x2li - G44).sin()
                        + d5220 * (xomi + xli - G52).sin()
                        + d5232 * (-xomi + xli - G52).sin()
                        + d5421 * (xomi + x2li - G54).sin()
                        + d5433 * (-xomi + x2li - G54).sin();
                    let xnddt = d2201 * (x2omi + xli - G22).cos()
                        + d2211 * (xli - G22).cos()
                        + d3210 * (xomi + xli - G32).cos()
                        + d3222 * (-xomi + xli - G32).cos()
                        + d5220 * (xomi + xli - G52).cos()
                        + d5232 * (-xomi + xli - G52).cos()
                        + 2.0
                            * (d4410 * (x2omi + x2li - G44).cos()
                                + d4422 * (x2li - G44).cos()
                                + d5421 * (xomi + x2li - G54).cos()
                                + d5433 * (-xomi + x2li - G54).cos());
                    (xndt, xldot, xnddt * xldot)
                }
            }
        };

        let (xndt, xldot, xnddt) = loop {
            let (xndt, xldot, xnddt) = derivatives(xli, xni, atime);
            if (t - atime).abs() < STEPP {
                break (xndt, xldot, xnddt);
            }
            xli += xldot * delt + xndt * STEP2;
            xni += xndt * delt + xnddt * STEP2;
            atime += delt;
        };

        let ft = t - atime;
        let nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
        let xl = xli + xldot * ft + xndt * ft * ft * 0.5;
        let dndt = nm - no;
        mean.nm = no + dndt;
        mean.mm = match resonance.terms {
            ResonanceTerms::Synchronous { .. } => xl - mean.nodem - mean.argpm + theta,
            ResonanceTerms::HalfDay { .. } => xl - 2.0 * mean.nodem + 2.0 * theta,
        };
    }

    /// `dpper`: adds lunar-solar periodics to the elements at time `t`.
    pub(super) fn periodics(&self, t: f64, mode: OpsMode, p: &mut Perturbed) {
        let ls = &self.lunar_solar;

        let zm = ls.zmos + ZNS * t;
        let zf = zm + 2.0 * ZES * zm.sin();
        let sinzf = zf.sin();
        let f2 = 0.5 * sinzf * sinzf - 0.25;
        let f3 = -0.5 * sinzf * zf.cos();
        let ses = ls.se2 * f2 + ls.se3 * f3;
        let sis = ls.si2 * f2 + ls.si3 * f3;
        let sls = ls.sl2 * f2 + ls.sl3 * f3 + ls.sl4 * sinzf;
        let sghs = ls.sgh2 * f2 + ls.sgh3 * f3 + ls.sgh4 * sinzf;
        let shs = ls.sh2 * f2 + ls.sh3 * f3;

        let zm = ls.zmol + ZNL * t;
        let zf = zm + 2.0 * ZEL * zm.sin();
        let sinzf = zf.sin();
        let f2 = 0.5 * sinzf * sinzf - 0.25;
        let f3 = -0.5 * sinzf * zf.cos();
        let sel = ls.ee2 * f2 + ls.e3 * f3;
        let sil = ls.xi2 * f2 + ls.xi3 * f3;
        let sll = ls.xl2 * f2 + ls.xl3 * f3 + ls.xl4 * sinzf;
        let sghl = ls.xgh2 * f2 + ls.xgh3 * f3 + ls.xgh4 * sinzf;
        let shll = ls.xh2 * f2 + ls.xh3 * f3;

        let pe = ses + sel;
        let pinc = sis + sil;
        let pl = sls + sll;
        let pgh = sghs + sghl;
        let ph = shs + shll;

        p.inclp += pinc;
        p.ep += pe;
        let sinip = p.inclp.sin();
        let cosip = p.inclp.cos();

        if p.inclp >= 0.2 {
            let ph = ph / sinip;
            let pgh = pgh - cosip * ph;
            p.argpp += pgh;
            p.nodep += ph;
            p.mp += pl;
            return;
        }

        // Lyddane's modification for low inclinations.
        let sinop = p.nodep.sin();
        let cosop = p.nodep.cos();
        let alfdp = sinip * sinop + (ph * cosop + pinc * cosip * sinop);
        let betdp = sinip * cosop + (-ph * sinop + pinc * cosip * cosop);
        p.nodep = fmod_tau(p.nodep);
        if p.nodep < 0.0 && mode == OpsMode::Afspc {
            p.nodep += TAU;
        }
        let xls = p.mp + p.argpp + pl + pgh + (cosip - pinc * sinip) * p.nodep;
        let xnoh = p.nodep;
        p.nodep = alfdp.atan2(betdp);
        if p.nodep < 0.0 && mode == OpsMode::Afspc {
            p.nodep += TAU;
        }
        if (xnoh - p.nodep).abs() > PI {
            if p.nodep < xnoh {
                p.nodep += TAU;
            } else {
                p.nodep -= TAU;
            }
        }
        p.mp += pl;
        p.argpp = xls - p.mp - cosip * p.nodep;
    }
}

/// Elements that `dpper` perturbs.
#[derive(Debug, Clone, Copy)]
pub(super) struct Perturbed {
    pub ep: f64,
    pub inclp: f64,
    pub nodep: f64,
    pub argpp: f64,
    pub mp: f64,
}

/// `dscom` at `tc = 0`.
fn dscom(epoch: f64, ep: f64, argpp: f64, inclp: f64, nodep: f64, np: f64) -> (Common, LunarSolar) {
    const C1SS: f64 = 2.986_479_7e-6;
    const C1L: f64 = 4.796_806_5e-7;
    const ZSINIS: f64 = 0.397_854_16;
    const ZCOSIS: f64 = 0.917_448_67;
    const ZCOSGS: f64 = 0.194_590_5;
    const ZSINGS: f64 = -0.980_884_58;

    let nm = np;
    let em = ep;
    let snodm = nodep.sin();
    let cnodm = nodep.cos();
    let sinomm = argpp.sin();
    let cosomm = argpp.cos();
    let sinim = inclp.sin();
    let cosim = inclp.cos();
    let emsq = em * em;
    let betasq = 1.0 - emsq;
    let rtemsq = betasq.sqrt();
    let xnoi = 1.0 / nm;

    // Lunar orbit geometry at epoch.
    let day = epoch + 18_261.5;
    let xnodce = (4.523_602_0 - 9.242_202_9e-4 * day) % TAU;
    let stem = xnodce.sin();
    let ctem = xnodce.cos();
    let zcosil = 0.913_751_64 - 0.035_680_96 * ctem;
    let zsinil = (1.0 - zcosil * zcosil).sqrt();
    let zsinhl = 0.089_683_511 * stem / zsinil;
    let zcoshl = (1.0 - zsinhl * zsinhl).sqrt();
    let gam = 5.835_151_4 + 0.001_944_368_0 * day;
    let zx = 0.397_854_16 * stem / zsinil;
    let zy = zcoshl * ctem + 0.917_448_67 * zsinhl * stem;
    let zx = gam + zx.atan2(zy) - xnodce;
    let zcosgl = zx.cos();
    let zsingl = zx.sin();

    /// Terms for one perturbing body: the Sun on the first pass, the Moon on the second.
    struct Body {
        s1: f64,
        s2: f64,
        s3: f64,
        s4: f64,
        s5: f64,
        s6: f64,
        s7: f64,
        z1: f64,
        z2: f64,
        z3: f64,
        z11: f64,
        z12: f64,
        z13: f64,
        z21: f64,
        z22: f64,
        z23: f64,
        z31: f64,
        z32: f64,
        z33: f64,
    }

    let body =
        |zcosg: f64, zsing: f64, zcosi: f64, zsini: f64, zcosh: f64, zsinh: f64, cc: f64| -> Body {
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

            let z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
            let z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
            let z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
            let z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * emsq;
            let z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * emsq;
            let z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * emsq;
            let z11 = -6.0 * a1 * a5 + emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
            let z12 = -6.0 * (a1 * a6 + a3 * a5)
                + emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
            let z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
            let z21 = 6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
            let z22 = 6.0 * (a4 * a5 + a2 * a6)
                + emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
            let z23 = 6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
            let z1 = z1 + z1 + betasq * z31;
            let z2 = z2 + z2 + betasq * z32;
            let z3 = z3 + z3 + betasq * z33;
            let s3 = cc * xnoi;
            let s2 = -0.5 * s3 / rtemsq;
            let s4 = s3 * rtemsq;
            let s1 = -15.0 * em * s4;
            let s5 = x1 * x3 + x2 * x4;
            let s6 = x2 * x3 + x1 * x4;
            let s7 = x2 * x4 - x1 * x3;
            Body {
                s1,
                s2,
                s3,
                s4,
                s5,
                s6,
                s7,
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
            }
        };

    let sun = body(ZCOSGS, ZSINGS, ZCOSIS, ZSINIS, cnodm, snodm, C1SS);
    let moon = body(
        zcosgl,
        zsingl,
        zcosil,
        zsinil,
        zcoshl * cnodm + zsinhl * snodm,
        snodm * zcoshl - cnodm * zsinhl,
        C1L,
    );

    let lunar_solar = LunarSolar {
        // Solar.
        se2: 2.0 * sun.s1 * sun.s6,
        se3: 2.0 * sun.s1 * sun.s7,
        si2: 2.0 * sun.s2 * sun.z12,
        si3: 2.0 * sun.s2 * (sun.z13 - sun.z11),
        sl2: -2.0 * sun.s3 * sun.z2,
        sl3: -2.0 * sun.s3 * (sun.z3 - sun.z1),
        sl4: -2.0 * sun.s3 * (-21.0 - 9.0 * emsq) * ZES,
        sgh2: 2.0 * sun.s4 * sun.z32,
        sgh3: 2.0 * sun.s4 * (sun.z33 - sun.z31),
        sgh4: -18.0 * sun.s4 * ZES,
        sh2: -2.0 * sun.s2 * sun.z22,
        sh3: -2.0 * sun.s2 * (sun.z23 - sun.z21),
        // Lunar.
        ee2: 2.0 * moon.s1 * moon.s6,
        e3: 2.0 * moon.s1 * moon.s7,
        xi2: 2.0 * moon.s2 * moon.z12,
        xi3: 2.0 * moon.s2 * (moon.z13 - moon.z11),
        xl2: -2.0 * moon.s3 * moon.z2,
        xl3: -2.0 * moon.s3 * (moon.z3 - moon.z1),
        xl4: -2.0 * moon.s3 * (-21.0 - 9.0 * emsq) * ZEL,
        xgh2: 2.0 * moon.s4 * moon.z32,
        xgh3: 2.0 * moon.s4 * (moon.z33 - moon.z31),
        xgh4: -18.0 * moon.s4 * ZEL,
        xh2: -2.0 * moon.s2 * moon.z22,
        xh3: -2.0 * moon.s2 * (moon.z23 - moon.z21),
        zmol: (4.719_967_2 + 0.229_971_50 * day - gam) % TAU,
        zmos: (6.256_583_7 + 0.017_201_977 * day) % TAU,
    };

    let common = Common {
        sinim,
        cosim,
        emsq,
        s1: moon.s1,
        s2: moon.s2,
        s3: moon.s3,
        s4: moon.s4,
        s5: moon.s5,
        ss1: sun.s1,
        ss2: sun.s2,
        ss3: sun.s3,
        ss4: sun.s4,
        ss5: sun.s5,
        sz1: sun.z1,
        sz3: sun.z3,
        sz11: sun.z11,
        sz13: sun.z13,
        sz21: sun.z21,
        sz23: sun.z23,
        sz31: sun.z31,
        sz33: sun.z33,
        z1: moon.z1,
        z3: moon.z3,
        z11: moon.z11,
        z13: moon.z13,
        z21: moon.z21,
        z23: moon.z23,
        z31: moon.z31,
        z33: moon.z33,
    };

    (common, lunar_solar)
}

/// Resonance coefficients for 12-hour, eccentric orbits (the `irez == 2`
/// branch of `dsinit`).
fn half_day_terms(
    em: f64,
    eccsq: f64,
    sinim: f64,
    cosim: f64,
    nm: f64,
    aonv: f64,
) -> ResonanceTerms {
    const ROOT22: f64 = 1.7891679e-6;
    const ROOT44: f64 = 7.3636953e-9;
    const ROOT54: f64 = 2.1765803e-9;
    const ROOT32: f64 = 3.7393792e-7;
    const ROOT52: f64 = 1.1428639e-7;

    let emsq = eccsq;
    let cosisq = cosim * cosim;
    let eoc = em * emsq;
    let g201 = -0.306 - (em - 0.64) * 0.440;

    let (g211, g310, g322, g410, g422, g520);
    if em <= 0.65 {
        g211 = 3.616 - 13.2470 * em + 16.2900 * emsq;
        g310 = -19.302 + 117.3900 * em - 228.4190 * emsq + 156.5910 * eoc;
        g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc;
        g410 = -41.122 + 242.6940 * em - 471.0940 * emsq + 313.9530 * eoc;
        g422 = -146.407 + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
        g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
    } else {
        g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc;
        g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc;
        g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc;
        g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc;
        g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
        g520 = if em > 0.715 {
            -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc
        } else {
            1464.74 - 4664.75 * em + 3763.64 * emsq
        };
    }

    let (g533, g521, g532);
    if em < 0.7 {
        g533 = -919.22770 + 4988.61 * em - 9064.77 * emsq + 5542.21 * eoc;
        g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
        g532 = -853.66600 + 4690.25 * em - 8624.77 * emsq + 5341.4 * eoc;
    } else {
        g533 = -37995.78 + 161616.52 * em - 229838.2 * emsq + 109377.94 * eoc;
        g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
        g532 = -40023.88 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
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
    let f542 =
        29.53125 * sinim * (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
    let f543 =
        29.53125 * sinim * (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));

    let xno2 = nm * nm;
    let ainv2 = aonv * aonv;
    let mut temp1 = 3.0 * xno2 * ainv2;
    let mut temp = temp1 * ROOT22;
    let d2201 = temp * f220 * g201;
    let d2211 = temp * f221 * g211;
    temp1 *= aonv;
    temp = temp1 * ROOT32;
    let d3210 = temp * f321 * g310;
    let d3222 = temp * f322 * g322;
    temp1 *= aonv;
    temp = 2.0 * temp1 * ROOT44;
    let d4410 = temp * f441 * g410;
    let d4422 = temp * f442 * g422;
    temp1 *= aonv;
    temp = temp1 * ROOT52;
    let d5220 = temp * f522 * g520;
    let d5232 = temp * f523 * g532;
    temp = 2.0 * temp1 * ROOT54;
    let d5421 = temp * f542 * g521;
    let d5433 = temp * f543 * g533;

    ResonanceTerms::HalfDay {
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
    }
}
