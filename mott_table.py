"""
Mott-relation and Wiedemann-Franz deviations for the PLLs of a strained, tilted
Dirac cone (valley xi=+), reproducing Table tabdev of arxiv_red.tex.

Units: v_x = a_0 = hbar = k_B = e = tau = L_y = 1, so energies are in v_x/a_0.
Overall factors of 1/(2 pi) are omitted, since they cancel in the relative deviations.

Method
------
For xi=+ the PLL energies (Eq. eqevals, small-B form not used) are
    E = eta k + si sqrt(Q),
    Q = rt^2 A^2 + (k rt + 1 - k^2/4) B,
    rt = 3/2 + k/4,  A = zeta n + (b + zeta)/2,  B = 2 zeta n + b + zeta,
with b = a_0^2 e B and zeta = b v_y/v_x. The group velocity is the exact k-derivative.
Each band is monotonic in |k| <= KM, so it has one Fermi point whenever
lo_n < mu < hi_n, with lo_n, hi_n the band energies at k = -KM, +KM, and
    g(eps) = sum_n v_n(eps) [H(eps - lo_n) - H(eps - hi_n)].
sigma, alpha, ell are the Fermi-function moments of g. They are evaluated
semi-analytically: a smooth part (Sommerfeld) plus the exact thermal smearing of
every step of g inside the mu window. This resolves T of order mK, for which a
brute-force k grid is impractical. Validated against direct k-integration at T=1e-4.

Limitation: bands must be monotonic in |k|<=KM (true for eta >~ 0.2 here).
"""
import numpy as np
from scipy.special import spence, expit

PI2 = np.pi ** 2
KB_EV = 8.617333e-5      # eV/K
UNIT_EV = 4.635          # hbar v_x / a_0 in eV (v_x = 1e6 m/s, a_0 = 0.142 nm)


def Qfun(b, ratio, k, n):
    z = b * ratio
    rt = 1.5 + k / 4
    ss = 1 - k ** 2 / 4
    A = z * n + (b + z) / 2
    B = 2 * z * n + b + z
    Q = rt ** 2 * A ** 2 + (k * rt + ss) * B
    dQ = 2 * rt * 0.25 * A ** 2 + (rt + 0.25 * k - k / 2) * B
    return Q, dQ


def eps_v(b, ratio, eta, si, k, n):
    Q, dQ = Qfun(b, ratio, k, n)
    sq = np.sqrt(Q)
    return eta * k + si * sq, eta + si * dQ / (2 * sq)


def delta_E(b, ratio=0.8):
    k = np.array([0.0])
    e1, _ = eps_v(b, ratio, 0.0, 1, k, np.array([1.0]))
    e2, _ = eps_v(b, ratio, 0.0, 1, k, np.array([2.0]))
    return float((e2 - e1)[0])


# thermal moments of the kernel K = -f'(x)
def fermi(x, T):
    return expit(-x / T)


def Kker(x, T):
    u = np.minimum(np.abs(x) / (2 * T), 300)
    return 1 / (4 * T * np.cosh(u) ** 2)


def M1(x0, T):          # int_{x0}^inf x K dx (even in x0)
    x = np.abs(x0)
    return x * fermi(x, T) + T * np.logaddexp(0, -x / T)


def M2(x0, T):          # int_{x0}^inf x^2 K dx
    x = np.abs(x0)
    u = np.exp(-np.minimum(x / T, 700))
    li = spence(1 + u)  # Li_2(-u)
    m = x ** 2 * fermi(x, T) + 2 * (x * T * np.log1p(u) - T ** 2 * li)
    return np.where(x0 >= 0, m, PI2 * T ** 2 / 3 - m)


def run(b, T, eta=0.3, KM=0.3, ratio=0.8, mulo=1e-3, muhi=3e-3, N=None):
    """Return mu grid and sigma, d sigma/d mu, alpha, ell on it."""
    muc = 0.5 * (mulo + muhi)
    if N is None:
        N = int(((muhi + eta * KM * 1.2 + 0.01) / 0.7) ** 2 / (2 * b * ratio)) + 400
    n = np.arange(N, dtype=float)
    pts = np.array([mulo, muhi])
    edges = []
    for si in (+1, -1):
        for kk, s in ((-KM, +1), (KM, -1)):     # lo edge (s=+1), hi edge (s=-1)
            e, _ = eps_v(b, ratio, eta, si, np.full_like(n, kk), n)
            fs = []
            for mu in pts:                        # v at the Fermi point, by Newton
                k = np.full_like(n, kk)
                for _ in range(10):
                    ee, vv = eps_v(b, ratio, eta, si, k, n)
                    k = np.clip(k - (ee - mu) / vv, -0.6, 0.6)
                fs.append(eps_v(b, ratio, eta, si, k, n)[1])
            fs = np.array(fs)
            c = (fs[1] - fs[0]) / (pts[1] - pts[0])   # d v/d eps (linear in the window)
            a = fs[0] + c * (muc - pts[0])            # v at mu_c
            edges.append((e, s, a, c))
    eall = np.concatenate([E[0] for E in edges])
    ein = eall[(eall > mulo - 60 * T) & (eall < muhi + 60 * T)]
    coarse = np.linspace(mulo, muhi, 2001)
    fine = [e + np.arange(-60, 60.01, 0.2) * T for e in ein]
    mu = np.unique(np.concatenate([coarse] + fine))
    mu = mu[(mu >= mulo) & (mu <= muhi)]
    sig = np.zeros_like(mu); dsig = np.zeros_like(mu)
    alp = np.zeros_like(mu); ell = np.zeros_like(mu)
    for (e, s, a, c) in edges:
        order = np.argsort(e); es = e[order]
        cum_a = np.concatenate([[0], np.cumsum(a[order])])
        cum_c = np.concatenate([[0], np.cumsum(c[order])])
        idx = np.searchsorted(es, mu)
        A = cum_a[idx]; C = cum_c[idx]
        # smooth (step) part
        sig += s * (A + C * (mu - muc)); dsig += s * C
        alp += -s * (PI2 * T / 3) * C
        ell += s * (PI2 * T / 3) * (A + C * (mu - muc))
        # exact thermal smearing minus the step part, for edges inside the window
        for j in np.where((e > mulo - 60 * T) & (e < muhi + 60 * T))[0]:
            x0 = e[j] - mu; H = (x0 < 0).astype(float)
            ap = a[j] + c[j] * (mu - muc)
            m0 = fermi(x0, T); m1 = M1(x0, T); m2 = M2(x0, T)
            sig += s * ((ap * m0 + c[j] * m1) - ap * H)
            alp += s * (-(1 / T) * (ap * m1 + c[j] * m2) + (PI2 * T / 3) * c[j] * H)
            ell += s * ((1 / T) * ap * m2 - (PI2 * T / 3) * ap * H)
            dsig += s * (ap * Kker(x0, T) + c[j] * (m0 - H) + c[j] * x0 * Kker(x0, T))
    return mu, sig, dsig, alp, ell


def dev(b, T, **kw):
    """L2 relative deviations (delta_M, delta_WF) and mean sigma."""
    mu, sig, dsig, alp, ell = run(b, T, **kw)
    w = np.gradient(mu)
    nrm = lambda f: np.sqrt(np.sum(w * f ** 2))
    alp_M = -(PI2 * T / 3) * dsig          # Mott
    ell_WF = (PI2 * T / 3) * sig           # Wiedemann-Franz
    return nrm(alp - alp_M) / nrm(alp), nrm(ell - ell_WF) / nrm(ell), sig.mean()


def mK_to_units(mK):
    return mK * 1e-3 * KB_EV / UNIT_EV


if __name__ == "__main__":
    mKs = (1, 2, 5, 10)
    print("T (mK) -> T in units of v_x/a_0:", [f"{mK_to_units(m):.3g}" for m in mKs])
    print("\nTable: columns = delta_M at T = 1,2,5,10 mK | delta_WF at same T")
    for b in (1e-5, 2e-5, 4e-5):
        r = [dev(b, mK_to_units(m))[:2] for m in mKs]
        print(f"b={b:.0e}  DeltaE={delta_E(b):.2e}  dM:",
              " ".join(f"{x[0]:.3f}" for x in r), "| dWF:",
              " ".join(f"{x[1]:.2e}" for x in r))
    T5 = mK_to_units(5)
    print("\nTilt dependence at b=2e-5, T=5 mK (eta<~0.15 not covered)")
    for eta in (0.2, 0.3, 0.45):
        dM, dW, s = dev(2e-5, T5, eta=eta)
        print(f"eta={eta}: dM={dM:.3f} dWF={dW:.2e} <sigma>={s:.3g}")
    print("\nk_x cutoff dependence at b=2e-5, T=5 mK")
    for KM in (0.3, 0.1):
        dM, dW, s = dev(2e-5, T5, KM=KM)
        print(f"KM={KM}: dM={dM:.3f} dWF={dW:.2e} <sigma>={s:.3g}")
