"""Deviations of the Mott relation (delta_M) and the Wiedemann-Franz law (delta_WF) for the table of the manuscript.

Self-contained (needs only numpy).  Valley xi=+, natural units, energies in units of v_x/a0, a0 = 0.142 nm.
PLL spectrum:  E/v_x = eta k + si sqrt(w L - rt^2 L^2/4),  L = 2|zeta|(n+1),  w = 1 + 3k/2,  rt = 3/2 + k/4,
zeta = (v_y/v_x) B~.  Bands are monotonic in |k| <= K, so sigma, alpha, ell are summed band by band:
    sigma = sum |v(kF)|/(2 pi) [C0(tb)-C0(ta)],  alpha = -sum g C1,  ell = sum g C2 / beta,   t = beta (E - mu),
where ta, tb are beta (E - mu) at the two ends of the window and C0, C1, C2 are cumulative integrals of
(-df/dE), (-df/dE) t, (-df/dE) t^2.
delta = sqrt(sum (x - y)^2 / sum y^2) over N points spread uniformly over mu in [MU_LO, MU_HI]; the Mott derivative
is a central difference with step h = FH*T, which must resolve the thermal width (hence FH = 0.1).

Usage:  python table.py [--N 60000] [--eta 0.3] [--K 0.3] [--T 1 2 5 10] [--B 1e-5 2e-5 4e-5] [--procs 2]
        T in millikelvin, B = a0^2 e B/hbar.  Output: one row per B with delta_M for each T, then delta_WF for each T.
Run time: about 20 minutes on two cores for the default table (N = 60000); use --N 6000 for a quick, noisier check.
"""
import argparse
import numpy as np
from multiprocessing import Pool

# ---- constants (hbar = c-based conversion as in the manuscript) ---------------------------------------
HBAR_C, A0_NM, KB_EV = 197.3269804, 0.142, 8.617333262e-5          # eV nm, nm, eV/K
B_UNIT_T = 1.054571817e-34/(1.602176634e-19*(A0_NM*1e-9)**2)      # hbar/(e a0^2) in tesla
V_GRAPHENE = 1e6/2.99792458e8                                      # v_x/c for v_x = 1e6 m/s
EPS0_EV = V_GRAPHENE*HBAR_C/A0_NM                                  # hbar v_x/a0 = 4.63 eV

MU_LO, MU_HI, RATIO, FH = 1e-3, 3e-3, 0.8, 0.1

# ---- Fermi-function cumulants -------------------------------------------------------------------------
_tg = np.linspace(-80, 80, 320001)
_kk = 0.25/np.cosh(_tg/2)**2
def _cum(f): return np.concatenate([[0], np.cumsum((f[1:] + f[:-1])/2*np.diff(_tg))])
_C0, _C1, _C2 = _cum(_kk), _cum(_kk*_tg), _cum(_kk*_tg**2)
def _tab(C, t): return np.interp(t, _tg, C)

# ---- PLL energy and group velocity --------------------------------------------------------------------
def _Ev(k, lam, eta, si):
    w = 1 + 1.5*k; r = 1.5 + k/4; arg = w*lam - r**2*lam**2/4
    E = eta*k + si*np.sqrt(arg)
    v = eta + si*(lam*1.5 - r*0.25*lam**2/2)/(2*np.sqrt(arg))
    return E, v

def compute(eta, Bt, K, beta, mus):
    """sigma, alpha, ell (tau-normalised) at the chemical potentials mus, and the max number of Fermi points."""
    z = abs(Bt*RATIO)
    nmax = int(((eta*K + np.max(np.abs(mus)) + 0.12)**2)/(2*z)) + 5
    lam = 2*z*(np.arange(nmax + 1) + 1)
    S = np.zeros(len(mus)); A = np.zeros(len(mus)); L = np.zeros(len(mus)); nfp = np.zeros(len(mus), int)
    for si in (1, -1):
        Ea, _ = _Ev(-K, lam, eta, si); Eb, _ = _Ev(K, lam, eta, si)
        for j0 in range(0, len(mus), 500):
            sl = slice(j0, j0 + 500); mu = mus[sl]
            ta = beta*(Ea[:, None] - mu[None, :]); tb = beta*(Eb[:, None] - mu[None, :])
            act = (ta < 80) & (tb > -80)                      # bands whose window ends lie within 80 T of mu
            idx = np.where(act.any(1))[0]
            if len(idx) == 0: continue
            lam_i = lam[idx][:, None]
            lo = np.full((len(idx), len(mu)), -K); hi = np.full((len(idx), len(mu)), K)
            for _ in range(42):                                # bisection for the Fermi point
                mid = (lo + hi)/2; E, _ = _Ev(mid, lam_i, eta, si)
                right = E < mu[None, :]; lo = np.where(right, mid, lo); hi = np.where(right, hi, mid)
            _, v = _Ev((lo + hi)/2, lam_i, eta, si); g = np.abs(v)/(2*np.pi)
            a = act[idx]; ta_ = ta[idx]; tb_ = tb[idx]
            c0, c1, c2 = (_tab(C, tb_) - _tab(C, ta_) for C in (_C0, _C1, _C2))
            S[sl] += np.sum(np.where(a, g*c0, 0), 0)
            A[sl] += -np.sum(np.where(a, g*c1, 0), 0)
            L[sl] += np.sum(np.where(a, g*c2, 0), 0)/beta
            nfp[sl] += np.sum(a & (ta_ < 0) & (tb_ > 0), 0)
    return S, A, L, nfp

def deviations(eta, Bt, K, T_mK, N):
    beta = EPS0_EV/(KB_EV*T_mK*1e-3); h = FH/beta
    mu = MU_LO + (np.arange(N) + 0.5)*((MU_HI - MU_LO)/N)
    S, A, L, nfp = compute(eta, Bt, K, beta, mu)
    Sp = compute(eta, Bt, K, beta, mu + h)[0]; Sm = compute(eta, Bt, K, beta, mu - h)[0]
    mott = -(np.pi**2/(3*beta))*(Sp - Sm)/(2*h); wf = (np.pi**2/(3*beta))*S
    L2 = lambda x, y: np.sqrt(np.sum((x - y)**2)/np.sum(y**2))
    return L2(A, mott), L2(L, wf), S.mean()

def _job(a): return deviations(*a)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--N', type=int, default=60000); p.add_argument('--eta', type=float, default=0.3)
    p.add_argument('--K', type=float, default=0.3)
    p.add_argument('--T', type=float, nargs='+', default=[1, 2, 5, 10])
    p.add_argument('--B', type=float, nargs='+', default=[1e-5, 2e-5, 4e-5])
    p.add_argument('--procs', type=int, default=2)
    a = p.parse_args()
    print('unit of energy %.3f eV; 1 mK -> %.3e; B unit %.0f T; eta=%g, |k|<=%g, N=%d'
          % (EPS0_EV, KB_EV*1e-3/EPS0_EV, B_UNIT_T, a.eta, a.K, a.N), flush=True)
    jobs = [(a.eta, B, a.K, T, a.N) for B in a.B for T in a.T]
    with Pool(a.procs) as pool: res = pool.map(_job, jobs)
    nT = len(a.T)
    print('B~  B[T]  | delta_M at T =', a.T, 'mK | delta_WF | mean sigma~ (window average)')
    for i, B in enumerate(a.B):
        r = res[i*nT:(i + 1)*nT]
        print('%.0e  %.2f  | ' % (B, B*B_UNIT_T) + ' '.join('%.3f' % x[0] for x in r)
              + ' | ' + ' '.join('%.2e' % x[1] for x in r)
              + ' | ' + ' '.join('%.4e' % x[2] for x in r), flush=True)
