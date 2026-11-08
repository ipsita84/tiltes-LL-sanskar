# -*- coding: utf-8 -*-
"""
sigma_B_inset.py -- sigma_xx vs B (both signs) with an inset zooming on small |B|.
Uses the model and plotting helpers of Codes.py (keep both files in the same folder).
Usage:  python sigma_B_inset.py [outdir]
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import Codes as C


def fig_sigma_B_inset(out):
    mus4 = [0.02, 0.03, 0.04, 0.05]
    # |B|: geometric grid down to 2e-6 (lower end of the validity window, Eq. (8)) + linear grid up to 1e-3
    Bp = np.unique(np.concatenate([np.geomspace(2e-6, 2e-5, 150), np.linspace(2e-5, 1e-3, 700)]))
    Bn = -Bp[::-1]                                   # negative-B branch (sigma is even in B)

    fig, axs = plt.subplots(1, 2, figsize=(14, 6), facecolor='white')
    for ax, beta in zip(axs, (200, 1000)):
        # sigma_xx(|B|) from Codes.py; the spectrum depends on B only through |zeta|,
        # so the value at -B equals the value at +B (checked: identical to machine precision)
        Sp = np.array([C.curves(C.ETA, B, C.WIN_A, beta, mus4, 0.07)[0][0] for B in Bp])
        Sn = Sp[::-1]

        ins = ax.inset_axes([0.60, 0.45, 0.37, 0.50])          # inset position: [x0, y0, width, height] in axes units
        for a, lw in ((ax, 2), (ins, 1.6)):
            for m in range(4):
                kw = dict(label=r'$\tilde{\mu}=%g$' % mus4[m]) if a is ax else {}
                a.plot(Bn, Sn[:, m], color=C.COLS[m], linewidth=lw, **kw)
                a.plot(Bp, Sp[:, m], color=C.COLS[m], linewidth=lw)

        # main panel
        C.style(ax, r'$\tilde{\mathcal{B}}$', r'$\tilde{\sigma}_{xx}$')
        ax.set_xlim(-1e-3, 1e-3); ax.set_ylim(-0.005, 0.29)
        ax.set_title(r'$\tilde{\beta}=%s,\ \eta_x=%g$' % (C.sci(beta), C.ETA), fontsize=20)
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        ax.xaxis.get_offset_text().set_fontsize(18)             # "1e-3"

        # inset
        ins.set_xlim(-2e-5, 2e-5); ins.set_ylim(0, 2.9)
        ins.tick_params(direction='in', labelsize=11, top=True, right=True)
        ins.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        ins.xaxis.get_offset_text().set_fontsize(14)            # "1e-5"
        ins.grid(True, linestyle='--', linewidth=0.7, color='#444444', alpha=0.7)
        ins.set_facecolor('white')
        for s in ins.spines.values(): s.set_linewidth(1.2)

    plt.tight_layout(); plt.subplots_adjust(bottom=0.24)
    C.legend_below_fig(fig, axs[0], 4)
    plt.savefig(f'{out}/sigma_B.png', dpi=170, bbox_inches='tight'); plt.close()


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else 'figs'
    os.makedirs(out, exist_ok=True)
    fig_sigma_B_inset(out)
