//! Vallado's `sgp4` function.

use std::f64::consts::{PI, TAU};

use super::deep_space::Perturbed;
use super::{MeanElements, Model, X2O3, fmod_tau};
use crate::error::PropagationError;
use crate::satellite::State;

impl Model {
    /// The state `t` minutes after epoch, in the TEME frame (km and km/s).
    pub(crate) fn propagate(&self, t: f64) -> Result<State, PropagationError> {
        let g = &self.gravity;
        let vkmpersec = g.radius_km * g.xke / 60.0;

        // ---- secular gravity and atmospheric drag ----
        let xmdf = self.mo + self.mdot * t;
        let argpdf = self.argpo + self.argpdot * t;
        let nodedf = self.nodeo + self.nodedot * t;
        let mut argpm = argpdf;
        let mut mm = xmdf;
        let t2 = t * t;
        let nodem = nodedf + self.nodecf * t2;
        let mut tempa = 1.0 - self.cc1 * t;
        let mut tempe = self.bstar * self.cc4 * t;
        let mut templ = self.t2cof * t2;

        if !self.isimp {
            let delomg = self.omgcof * t;
            let delmtemp = 1.0 + self.eta * xmdf.cos();
            let delm = self.xmcof * (delmtemp * delmtemp * delmtemp - self.delmo);
            let temp = delomg + delm;
            mm = xmdf + temp;
            argpm = argpdf - temp;
            let t3 = t2 * t;
            let t4 = t3 * t;
            tempa = tempa - self.d2 * t2 - self.d3 * t3 - self.d4 * t4;
            tempe += self.bstar * self.cc5 * (mm.sin() - self.sinmao);
            templ = templ + self.t3cof * t3 + t4 * (self.t4cof + t * self.t5cof);
        }

        let mut mean = MeanElements {
            em: self.ecco,
            argpm,
            inclm: self.inclo,
            mm,
            nodem,
            nm: self.no_unkozai,
        };
        if let Some(deep) = &self.deep {
            deep.secular(t, &mut mean, self);
        }
        let MeanElements {
            mut em,
            mut argpm,
            inclm,
            mut mm,
            mut nodem,
            mut nm,
        } = mean;

        if nm <= 0.0 {
            return Err(PropagationError::MeanMotion(nm));
        }
        let am = (g.xke / nm).powf(X2O3) * tempa * tempa;
        nm = g.xke / am.powf(1.5);
        em -= tempe;
        if !(-0.001..1.0).contains(&em) {
            return Err(PropagationError::Eccentricity(em));
        }
        if em < 1.0e-6 {
            em = 1.0e-6;
        }
        mm += self.no_unkozai * templ;
        let xlm = mm + argpm + nodem;
        nodem = fmod_tau(nodem);
        argpm %= TAU;
        let xlm = xlm % TAU;
        mm = (xlm - argpm - nodem) % TAU;

        // ---- lunar-solar periodics ----
        let mut p = Perturbed {
            ep: em,
            inclp: inclm,
            nodep: nodem,
            argpp: argpm,
            mp: mm,
        };
        let (sinip, cosip, aycof, xlcof, con41, x1mth2, x7thm1);
        if let Some(deep) = &self.deep {
            deep.periodics(t, self.mode, &mut p);
            if p.inclp < 0.0 {
                p.inclp = -p.inclp;
                p.nodep += PI;
                p.argpp -= PI;
            }
            if p.ep < 0.0 || p.ep > 1.0 {
                return Err(PropagationError::PerturbedEccentricity(p.ep));
            }
            // Coefficients that depend on inclination, which is no longer constant.
            sinip = p.inclp.sin();
            cosip = p.inclp.cos();
            aycof = -0.5 * g.j3oj2 * sinip;
            xlcof = super::xlcof(g.j3oj2, sinip, cosip);
            let cosisq = cosip * cosip;
            con41 = 3.0 * cosisq - 1.0;
            x1mth2 = 1.0 - cosisq;
            x7thm1 = 7.0 * cosisq - 1.0;
        } else {
            sinip = inclm.sin();
            cosip = inclm.cos();
            aycof = self.aycof;
            xlcof = self.xlcof;
            con41 = self.con41;
            x1mth2 = self.x1mth2;
            x7thm1 = self.x7thm1;
        }
        let Perturbed {
            ep,
            inclp: xincp,
            nodep,
            argpp,
            mp,
        } = p;

        // ---- long-period periodics ----
        let axnl = ep * argpp.cos();
        let temp = 1.0 / (am * (1.0 - ep * ep));
        let aynl = ep * argpp.sin() + temp * aycof;
        let xl = mp + argpp + nodep + temp * xlcof * axnl;

        // ---- solve Kepler's equation ----
        let u = (xl - nodep) % TAU;
        let mut eo1 = u;
        let mut tem5: f64 = 9999.9;
        let (mut sineo1, mut coseo1) = (0.0, 0.0);
        let mut ktr = 1;
        // Like the reference code, the loop exits with sin/cos of the
        // second-to-last iterate; the result depends on that.
        while tem5.abs() >= 1.0e-12 && ktr <= 10 {
            (sineo1, coseo1) = eo1.sin_cos();
            tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
            tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
            if tem5.abs() >= 0.95 {
                tem5 = 0.95_f64.copysign(tem5);
            }
            eo1 += tem5;
            ktr += 1;
        }

        // ---- short-period periodics ----
        let ecose = axnl * coseo1 + aynl * sineo1;
        let esine = axnl * sineo1 - aynl * coseo1;
        let el2 = axnl * axnl + aynl * aynl;
        let pl = am * (1.0 - el2);
        if pl < 0.0 {
            return Err(PropagationError::SemiLatusRectum(pl));
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
        let temp1 = 0.5 * g.j2 * temp;
        let temp2 = temp1 * temp;

        let mrt = rl * (1.0 - 1.5 * temp2 * betal * con41) + 0.5 * temp1 * x1mth2 * cos2u;
        su -= 0.25 * temp2 * x7thm1 * sin2u;
        let xnode = nodep + 1.5 * temp2 * cosip * sin2u;
        let xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
        let mvt = rdotl - nm * temp1 * x1mth2 * sin2u / g.xke;
        let rvdot = rvdotl + nm * temp1 * (x1mth2 * cos2u + 1.5 * con41) / g.xke;

        if mrt < 1.0 {
            return Err(PropagationError::Decayed(mrt));
        }

        // ---- orientation vectors ----
        let (sinsu, cossu) = su.sin_cos();
        let (snod, cnod) = xnode.sin_cos();
        let (sini, cosi) = xinc.sin_cos();
        let xmx = -snod * cosi;
        let xmy = cnod * cosi;
        let ux = xmx * sinsu + cnod * cossu;
        let uy = xmy * sinsu + snod * cossu;
        let uz = sini * sinsu;
        let vx = xmx * cossu - cnod * sinsu;
        let vy = xmy * cossu - snod * sinsu;
        let vz = sini * cossu;

        let mr = mrt * g.radius_km;
        Ok(State {
            position: [mr * ux, mr * uy, mr * uz],
            velocity: [
                (mvt * ux + rvdot * vx) * vkmpersec,
                (mvt * uy + rvdot * vy) * vkmpersec,
                (mvt * uz + rvdot * vz) * vkmpersec,
            ],
        })
    }
}
