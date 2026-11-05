# -*- coding: utf-8 -*-
"""
Codes_new.py -- all numerical figures of the manuscript (valley xi=+), corrected formulas.

Formulas (identical to Code-corrected.py, FULL_ENERGY=True, i.e. Eq. (12) of the manuscript):
    zeta = e B v_y/v_x ,  Lambda_n = 2|zeta|(n+1) ,  w = 1 + 3 a0 kx/2 ,  rt = a0(3/2 + a0 kx/4)
    E_n^si/v_x = eta kx + si sqrt(w Lambda_n - rt^2 Lambda_n^2/4) ,   n = 0, 1, 2, ...
    level n exists only if rt^2 |zeta| (n+1) <= w
    v = exact d/dkx of E.   sigma, alpha, ell = integrals of Eq. (alphaell), tau-normalised.
The kx integrals are evaluated with a fast histogram + Fermi-kernel convolution (checked against
scipy.quad of Code-corrected.py to ~1e-4).  The spectrum depends on B only through |zeta|,
so all coefficients are even in B and only B~>0 is plotted.

UNITS.  Natural units hbar = k_B = 1.  Dimensionless quantities
    mu~ = a0 mu/v_x ,  beta~ = v_x/(a0 T) ,  B~ = a0^2 e B/hbar ,  k~ = a0 kx
    sigma~ = L_y sigma/(tau v_x) ,  alpha~ = L_y alpha/(tau v_x) ,  ell~ = a0 L_y ell/(tau v_x^2)
Conversion to SI (v = v_x/c is the dimensionless velocity, a0 = 0.142 nm):
    eps0 = hbar v_x/a0 = v * 197.327 eV nm / 0.142 nm = v * 1389.6 eV     (energy unit)
    mu = mu~ eps0 ,   T = eps0/(beta~ k_B) ,   B = B~ * hbar/(e a0^2) = B~ * 3.264e4 T
    tau = 1e-13 s = 151.9 eV^-1 ,  hbar/tau = 6.58 meV  (tau only normalises the plotted quantities)

REALISTIC PARAMETERS USED (v = 0.0015 -> v_x = 4.5e5 m/s, eps0 = 2.08 eV):
    B~ = 2e-4 (6.5 T) ; beta~ = 1e3, 5e3, 1e4 (T = 24, 4.8, 2.4 K) ; beta~ = 200 (121 K) ;
    |mu~| <= 0.048 (|mu| <~ 100 meV) ; eta_x = 0.7 ; v_y/v_x = 0.8
    level spacing at kx=0 between n=0..6 and n+1: 15 -> 6.8 meV  (> hbar/tau = 6.6 meV)
KX WINDOWS (modelling choice, stated in the captions):
    [0, 1e-3]        for sigma/alpha/ell vs mu, sigma vs B, Mott/WF
    [-0.023, 0.023]  for sigma vs eta_x (sigma(eta_x=0) vanishes in the narrow window)
Usage:  python Codes_new.py [outdir] [fig1 fig2 fig3 ...]     (default: all, outdir = ./figs)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

# ------------------------------------------------------------------ constants / SI conversion
a0_nm, hbar_c, hbar_eVs, kB_eV = 0.142, 197.3269804, 6.582119569e-16, 8.617333262e-5
B_unit_T = 1.054571817e-34/(1.602176634e-19*(a0_nm*1e-9)**2)      # hbar/(e a0^2) [T]
def eps0_eV(v):  return v*hbar_c/a0_nm
def tau_nat_eV(tau_s): return tau_s/hbar_eVs

# ------------------------------------------------------------------ parameters
V_DIM  = 0.0015          # v_x/c
RATIO  = 0.8             # v_y/v_x
ETA    = 0.7
B0     = 2e-4            # B~ (6.5 T)
BETAS  = (1000, 5000, 10000)
MU_MAX = 0.048
WIN_A  = (0.0, 1e-3)
WIN_B  = (-0.023, 0.023)
COLS   = ['#1f77b4', '#E66100', '#74C476', '#9467bd']

# ------------------------------------------------------------------ PLL model (xi = +, a0 = 1)
def rt(k):   return 1.5 + k/4.0
def term(k): return k*rt(k) + (1 + k/2)*(1 - k/2)          # = 1 + 3k/2

def band(k, n, sign, eta, Bt, ratio=RATIO, full=True):
    z = abs(Bt*ratio); lam = 2*z*(n + 1.0)
    w = term(k)[None, :]; r = rt(k)[None, :]
    ok = (w > 0) & (r**2*z*(n + 1.0) <= w)
    arg = w*lam - (r**2*lam**2/4 if full else 0.0)
    ok &= arg > 0
    argc = np.where(ok, arg, 1.0)
    E = eta*k[None, :] + sign*np.sqrt(argc)
    darg = lam*1.5 - (r*0.25*lam**2/2 if full else 0.0)
    v = eta + sign*darg/(2*np.sqrt(argc))
    return np.where(ok, E, np.nan), np.where(ok, v, np.nan)

def spectral_weight(eta, Bt, kwin, Emax, dE, ratio=RATIO, nmax=None, full=True):
    """W(E) = sum_{n,si} int dk/(2 pi) v^2 delta(E - E_nk), on E_j = -Emax + (j+1/2) dE."""
    k1, k2 = kwin
    nb = int(round(2*Emax/dE)); W = np.zeros(nb)
    dk0 = dE/(2*(abs(eta) + 0.3))
    nk = int(np.ceil((k2 - k1)/dk0)); dk = (k2 - k1)/nk
    k = k1 + (np.arange(nk) + 0.5)*dk                           # midpoint rule
    z = abs(Bt*ratio)
    if nmax is None:
        Kmax = max(abs(k1), abs(k2))
        nmax = int(np.ceil((Emax + abs(eta)*Kmax + 0.02)**2/(2*z))) + 2
    chunk = max(1, int(4e6//nk))
    for sign in (1, -1):
        for n0 in range(0, nmax + 1, chunk):
            n = np.arange(n0, min(n0 + chunk, nmax + 1))[:, None]
            E, v = band(k, n, sign, eta, Bt, ratio, full)
            m = np.isfinite(E)
            idx = np.floor((E[m] + Emax)/dE).astype(np.int64)
            g = (idx >= 0) & (idx < nb)
            W += np.bincount(idx[g], weights=(v[m]**2*dk/(2*np.pi))[g], minlength=nb)
    Eg = -Emax + (np.arange(nb) + 0.5)*dE
    return Eg, W/dE, nmax

def responses(Eg, Wd, dE, beta, mus):
    """sigma~, alpha~, ell~ at chemical potentials mus (Fermi-kernel convolution)."""
    half = int(np.ceil(25.0/beta/dE)); x = np.arange(-half, half + 1)*dE
    K0 = beta/(4*np.cosh(np.clip(beta*x/2, -300, 300))**2)      # -df/dE
    Wb = Wd*dE
    s = fftconvolve(Wb, K0, mode='same')
    a = fftconvolve(Wb, beta*K0*x, mode='same')                 # alpha: +beta f' (E-mu) -> see sign check in text
    l = fftconvolve(Wb, beta*K0*x**2, mode='same')
    return tuple(np.interp(mus, Eg, q) for q in (s, a, l))

def curves(eta, Bt, kwin, beta, mus, Emax, fine=40):
    dE = min(1/(fine*beta), 1e-4)
    Eg, W, _ = spectral_weight(eta, Bt, kwin, Emax + 25/beta, dE)
    return responses(Eg, W, dE, beta, np.asarray(mus)), (Eg, W, dE)

# ------------------------------------------------------------------ plotting style (old style, with panel titles and legends)
def style(ax, xl, yl, fs=22):
    ax.set_xlabel(xl, fontsize=fs); ax.set_ylabel(yl, fontsize=fs)
    ax.tick_params(direction='in', top=True, right=True, labelsize=15); ax.minorticks_on()
    ax.grid(True, which='major', linestyle='--', linewidth=0.9, color='#444444', alpha=0.8)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.6, color='#666666', alpha=0.6)

def sci(x):
    """LaTeX string for x as 10^p or m x 10^p."""
    if abs(x) < 1e3: return r'%g' % x
    p = int(np.floor(np.log10(abs(x)) + 1e-9)); m = x/10**p
    return r'10^{%d}' % p if abs(m - 1) < 1e-9 else r'%g\times10^{%d}' % (m, p)

def legend_below(ax, ncol, fs=16, y=-0.24):
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, y), ncol=ncol, fontsize=fs,
              frameon=True, handlelength=2.2, columnspacing=1.2)

def legend_below_fig(fig, ax, ncol, fs=20):
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 0.11), ncol=ncol, fontsize=fs,
               frameon=True, handlelength=2.2, columnspacing=1.2)

# ------------------------------------------------------------------ Fig. transport: sigma, alpha, ell vs mu
def fig_transport(out):
    mus = np.linspace(-MU_MAX, MU_MAX, 4000); res = {}
    for b in BETAS:
        res[b], _ = curves(ETA, B0, WIN_A, b, mus, MU_MAX + 0.02)
    lab = [r'$\tilde{\sigma}_{xx}$', r'$\tilde{\alpha}_{xx}$', r'$\tilde{\ell}_{xx}$']
    for i, name in enumerate(('sigma', 'alpha', 'ell')):
        fig, ax = plt.subplots(figsize=(8, 4.6), facecolor='white')
        for b, c in zip(BETAS, COLS):
            ax.plot(mus, res[b][i], color=c, linewidth=1.8, label=r'$\tilde{\beta}=%s$' % sci(b))
        style(ax, r'$\tilde{\mu}$', lab[i], 20); ax.set_xlim(-MU_MAX, MU_MAX)
        ax.set_title(r'$\tilde{\mathcal{B}}=%s,\ \eta_x=%g$' % (sci(B0), ETA), fontsize=20)
        legend_below(ax, len(BETAS))
        if i > 0: ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.tight_layout(); plt.savefig(f'{out}/{name}.png', dpi=200, bbox_inches='tight'); plt.close()

# ------------------------------------------------------------------ sigma vs B
def fig_sigma_B(out):
    mus4 = [0.02, 0.03, 0.04, 0.05]; Bs = np.linspace(2e-5, 1e-3, 700)
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), facecolor='white')
    for ax, beta in zip(axs, (200, 1000)):
        S = np.array([curves(ETA, B, WIN_A, beta, mus4, 0.07)[0][0] for B in Bs])
        for m in range(4): ax.plot(Bs, S[:, m], color=COLS[m], linewidth=2, label=r'$\tilde{\mu}=%g$' % mus4[m])
        style(ax, r'$\tilde{\mathcal{B}}$', r'$\tilde{\sigma}_{xx}$'); ax.set_xlim(0, 1e-3)
        ax.set_title(r'$\tilde{\beta}=%s,\ \eta_x=%g$' % (sci(beta), ETA), fontsize=25)
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0)); ax.set_ylim(-0.02*S.max(), 1.04*S.max())
    plt.tight_layout(); plt.subplots_adjust(bottom=0.24); legend_below_fig(fig, axs[0], 4)
    plt.savefig(f'{out}/sigma_B.png', dpi=170, bbox_inches='tight'); plt.close()

# ------------------------------------------------------------------ sigma vs eta_x
def fig_sigma_eta(out):
    mus4 = [0.01, 0.02, 0.03, 0.04]; etas = np.linspace(-0.99, 0.99, 150)
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), facecolor='white')
    for ax, beta in zip(axs, (200, 1000)):
        S = np.array([curves(e, B0, WIN_B, beta, mus4, 0.08, fine=20)[0][0] for e in etas])
        for m in range(4): ax.plot(etas, S[:, m], color=COLS[m], linewidth=3, label=r'$\tilde{\mu}=%g$' % mus4[m])
        style(ax, r'$\eta_x$', r'$\tilde{\sigma}_{xx}$'); ax.set_xlim(-1, 1); ax.set_ylim(-0.02*S.max(), 1.04*S.max())
        ax.set_title(r'$\tilde{\beta}=%s$' % sci(beta), fontsize=25)
    plt.tight_layout(); plt.subplots_adjust(bottom=0.24); legend_below_fig(fig, axs[0], 4)
    plt.savefig(f'{out}/sigma_etax1.png', dpi=170, bbox_inches='tight'); plt.close()

# ------------------------------------------------------------------ Mott / Wiedemann-Franz
def fig_mott(out):
    beta = 1e4; lo, hi = 0.03, 0.048; h = 4e-6
    mu = np.linspace(lo, hi, 1800); dE = 1/(40*beta)
    Eg, W, _ = spectral_weight(ETA, B0, WIN_A, hi + 0.01 + 25/beta, dE)
    s, a, l = responses(Eg, W, dE, beta, mu)
    sp = responses(Eg, W, dE, beta, mu + h)[0]; sm = responses(Eg, W, dE, beta, mu - h)[0]
    mott = -(np.pi**2/(3*beta))*(sp - sm)/(2*h); wf = (np.pi**2/(3*beta))*s
    for arr, ref, name, yl in ((a, mott, 'mot1', r'$\tilde{\alpha}_{xx}$'), (l, wf, 'mot2', r'$\tilde{\ell}_{xx}$')):
        fig, ax = plt.subplots(figsize=(10, 7), facecolor='white')
        second = 'Mott' if name == 'mot1' else 'WF'
        ax.plot(mu, arr, color='blue', linewidth=2, label='Integral')
        ax.plot(mu, ref, color='red', linestyle='--', linewidth=2, label=second)
        style(ax, r'$\tilde{\mu}$', yl); ax.set_xlim(lo, hi)
        legend_below(ax, 2, 22, y=-0.16)
        if name == 'mot2': ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.tight_layout(); plt.savefig(f'{out}/{name}.png', dpi=150, bbox_inches='tight'); plt.close()
    L2 = lambda x, y: np.sqrt(np.sum((x - y)**2)/np.sum(y**2))
    print('relative L2 deviation: Mott %.3f, WF %.3f' % (L2(a, mott), L2(l, wf)))

def print_si():
    e0 = eps0_eV(V_DIM)
    print('eps0 = %.3f eV, v_x = %.2e m/s, B~=%g -> %.2f T, hbar/tau = %.2f meV' %
          (e0, V_DIM*2.998e8, B0, B0*B_unit_T, 1e3/tau_nat_eV(1e-13)))
    print('beta~ -> T [K]:', {b: round(e0/b/kB_eV, 2) for b in (200, 1000, 5000, 10000)}, '; mu_max = %.0f meV' % (MU_MAX*e0*1e3))

if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else 'figs'
    which = sys.argv[2:] or ['fig_transport', 'fig_sigma_B', 'fig_sigma_eta', 'fig_mott']
    os.makedirs(out, exist_ok=True); print_si()
    for w in which: print('running', w); globals()[w](out)
