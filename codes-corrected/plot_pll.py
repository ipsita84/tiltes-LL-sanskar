#!/usr/bin/env python3
"""Regenerates figPLL.pdf and figPLL2.pdf with the corrected PLL spectrum (arxiv_ss.tex).

Dimensionless variables: k = a0*kx, E = a0*E/vx, B = a0^2*e*B, rat = vy/vx.
Eq. (12):  E/vx = eta*k + si*sqrt( w*Lam - rt^2*Lam^2/4 ),
           Lam = 2*|zeta|*(n+1),  |zeta| = |B|*rat,  n = 0,1,2,...
           w = xi*k*rt + s_{-xi}*s_xi = 1 + 3*xi*k/2,  rt = 3/2 + xi*k/4.
Level n exists only where rt^2*|zeta|*(n+1) <= w  (otherwise the curve ends).
The red line is the zero mode E = eta*k (exists in both valleys).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------ parameters ------------------------------
REGIME = "valid"      # "valid": inside the validity window of the manuscript
                      # "old-range": old k-windows and B=1e-2, only for comparison
TITLES = True         # brown panel titles with the parameters (as in the manuscript figures)
FIGS = {              # name: (vy/vx, label of vy/vx in the title, draw zero-mode line)
    "figPLL":  (1/9, "1/9", False),
    "figPLL2": (1.0, "1", True),
}
ETA = 0.9
if REGIME == "valid":
    B, NLEV, KMAX = 2e-4, 20, 0.3          # |a0 kx| <= KMAX, levels n = 0..NLEV-1
    KWIN = {+1: (-KMAX, KMAX), -1: (-KMAX, KMAX)}
    NL = {+1: NLEV, -1: NLEV}
else:
    B = 1e-2
    KWIN = {+1: (-2.0, 4.0), -1: (-3.5, 4.0)}
    NL = {+1: 21, -1: 36}
OUTDIR = os.path.dirname(os.path.abspath(__file__))

# ------------------------------ model -----------------------------------
def rt(xi, k):  return 1.5 + xi * k / 4
def wxi(xi, k): return 1 + 1.5 * xi * k            # = xi*k*rt + s_{-xi}*s_xi

def en(xi, eta, k, n, sign, rat, B):
    r, w = rt(xi, k), wxi(xi, k)
    z = abs(B) * rat
    lam = 2 * z * (n + 1)
    ok = (w > 0) & (r**2 * z * (n + 1) <= w)
    arg = np.where(ok, w * lam - r**2 * lam**2 / 4, np.nan)
    return eta * k + sign * np.sqrt(arg)

# ------------------------------ style (as the Mathematica figures) -------
BLUE, ORANGE, RED = (0, 0, 1), (1, 0.5, 0), (1, 0, 0)
BROWN_DARK = (0.4, 0.27, 0.13)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "Times New Roman", "DejaVu Serif"],
    "font.weight": "bold",
    "mathtext.fontset": "stix",
    "axes.linewidth": 0.6,
})
GRID = dict(color="0.72", lw=0.35)

def blabel(B):
    m, e = f"{B:.0e}".split("e")
    m, e = int(float(m)), int(e)
    return r"10^{%d}" % e if m == 1 else r"%d\times10^{%d}" % (m, e)

def panel(ax, xi, eta, rat, ratlab, zero):
    lo, hi = KWIN[xi]
    k = np.linspace(lo, hi, 2001)
    for n in range(NL[xi]):
        ax.plot(k, en(xi, eta, k, n, +1, rat, B), color=BLUE, lw=0.9)
        ax.plot(k, en(xi, eta, k, n, -1, rat, B), color=ORANGE, lw=0.9)
    if zero:
        ax.plot(k, eta * k, color=RED, lw=1.0)
    ax.grid(True, **GRID)
    ax.set_axisbelow(True)
    ax.set_xlabel(r"$\mathbfit{a}_{\mathbf{0}}\,\mathbfit{k}_{\mathbfit{x}}$", fontsize=10)
    ax.set_ylabel(r"PLL energy $\mathbf{\times}\,\mathbfit{a}_{\mathbf{0}}/\mathbfit{v}_{\mathbfit{x}}$", fontsize=10)
    ax.tick_params(labelsize= 8, colors="0.3", length=2.5, width=0.5, pad=2)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontweight("bold")
    ax.tick_params(axis="both", which="both", top=False, right=False)
    if TITLES:
        sg = "+" if xi > 0 else "\u2212"
        e = f"{eta:g}".replace("-", "\u2212")
        ax.set_title(r"{$\xi$, $\eta_x$, $v_y/v_x$, $a_0^2\,$$\mathscr{B}$} = {%s, %s, %s, $%s$}" % (sg, e, ratlab, blabel(B)),
                     fontsize=10, color=BROWN_DARK, fontweight="bold")

def make(name, rat, ratlab, zero):
    fig, axs = plt.subplots(2, 2, figsize=(454/72, 323/72))   # same page size as the old PDFs
    for i, eta in enumerate((ETA, -ETA)):
        for j, xi in enumerate((+1, -1)):
            panel(axs[i, j], xi, eta, rat, ratlab, zero)
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.095, top=0.93,
                        wspace=0.34, hspace=0.55)
    fig.savefig(os.path.join(OUTDIR, name + ".pdf"))
    plt.close(fig)
    print(f"{name}: vy/vx={rat:.4g}, a0^2 eB={B:g}, |eta_x|={ETA}, "
          f"|a0 kx|<={KWIN[+1][1]:g}, levels n=0..{NL[+1]-1}, zero-mode line={zero}")

if __name__ == "__main__":
    for nm, (r, rl, z) in FIGS.items():
        make(nm, r, rl, z)
