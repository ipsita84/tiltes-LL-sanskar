#!/usr/bin/env python3
"""Combines figPLL (a) and figPLL2 (b) of plot_pll.py into one figure, fig_PLL_combined.{pdf,png}.
Uses the model, parameters and style of plot_pll.py unchanged (edit them there)."""
import os, sys
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_pll as P

OUT = os.path.dirname(os.path.abspath(__file__))
W, H = 454/72, 323/72                       # size of one block (as the old PDFs)
fig = plt.figure(figsize=(W, 2*H + 0.5))
subs = fig.subfigures(2, 1, height_ratios=[1, 1], hspace=0.0)
for sf, (name, (rat, ratlab, zero)), tag in zip(subs, P.FIGS.items(), ("(a)", "(b)")):
    axs = sf.subplots(2, 2)
    for i, eta in enumerate((P.ETA, -P.ETA)):
        for j, xi in enumerate((+1, -1)):
            P.panel(axs[i, j], xi, eta, rat, ratlab, zero)
    sf.subplots_adjust(left=0.095, right=0.985, bottom=0.115, top=0.93, wspace=0.34, hspace=0.55)
    sf.text(0.54, 0.012, tag, ha="center", va="bottom", fontsize=11, fontweight="normal")
#fig.savefig(os.path.join(OUT, "fig_PLL_combined.pdf"))                 # vector
fig.savefig(os.path.join(OUT, "fig_PLL_combined.png"), dpi=600)        # raster, 600 dpi
print("written fig_PLL_combined.png")
