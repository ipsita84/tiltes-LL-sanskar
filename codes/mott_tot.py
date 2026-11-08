# -*- coding: utf-8 -*-
"""mott_side.py -- Mott relation (left) and Wiedemann-Franz law (right) side by side in one figure.
Uses Codes.py (same folder).  Usage: python mott_side.py [outdir]"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import Codes as C

LAB, TICK, OFF, LEG = 34, 22, 20, 28          # font sizes: axis labels, ticks, "1e-5", legend

def fig_mott_side(out):
    beta = 1e4; lo, hi = 0.03, 0.048; h = 4e-6
    mu = np.linspace(lo, hi, 1800); dE = 1/(40*beta)
    Eg, W, _ = C.spectral_weight(C.ETA, C.B0, C.WIN_A, hi + 0.01 + 25/beta, dE)
    s, a, l = C.responses(Eg, W, dE, beta, mu)
    sp = C.responses(Eg, W, dE, beta, mu + h)[0]; sm = C.responses(Eg, W, dE, beta, mu - h)[0]
    mott = -(np.pi**2/(3*beta))*(sp - sm)/(2*h); wf = (np.pi**2/(3*beta))*s

    fig, axs = plt.subplots(1, 2, figsize=(20, 8), facecolor='white')
    for ax, arr, ref, second, yl in ((axs[0], a, mott, 'Mott', r'$\tilde{\alpha}_{xx}$'),
                                     (axs[1], l, wf, 'WF', r'$\tilde{\ell}_{xx}$')):
        ax.plot(mu, arr, color='blue', linewidth=2.5, label='Integral')
        ax.plot(mu, ref, color='red', linestyle='--', linewidth=2.5, label=second)
        C.style(ax, r'$\tilde{\mu}$', yl, LAB); ax.set_xlim(lo, hi)
        ax.tick_params(labelsize=TICK, length=6)
        ax.set_xticks(np.arange(0.030, 0.0481, 0.004))                 # fewer ticks so the larger labels fit
        if second == 'WF':
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax.yaxis.get_offset_text().set_fontsize(OFF)                # "1e-5"
        C.legend_below(ax, 2, LEG, y=-0.22)
    plt.tight_layout(w_pad=3)
    plt.savefig(f'{out}/mott.png', dpi=150, bbox_inches='tight'); plt.close()

if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else 'figs'
    os.makedirs(out, exist_ok=True); fig_mott_side(out)
