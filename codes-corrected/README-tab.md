# Numerical code for "Linear response from tilted Dirac cones under strain-induced pseudomagnetic fields"

Valley xi = +, PLL spectrum of the manuscript (full expression with the finite upper bound on n), natural units hbar = k_B = 1, energies in units of v_x/a0 with a0 = 0.142 nm.

## Files

| File | Purpose |
|---|---|
| `Codes.py` | Generates the figures (sigma, alpha, ell versus mu; sigma versus B and eta_x; Mott and Wiedemann-Franz comparison). |
| `table.py` | Computes the table of delta_M and delta_WF. Self-contained (needs only numpy); does not import `Codes.py`. |
| `arxiv_ss.tex`, `arxiv_ss_marked.tex` | Manuscript with the updated table and text; the changes are in red in the marked version. |

`Codes.py` needs numpy, scipy and matplotlib.

## Usage

Figures (default output directory `figs`):

    python Codes.py [outdir] [fig_transport fig_sigma_B fig_sigma_eta fig_mott]

Table (T in millikelvin, B = a0^2 e B/hbar):

    python table.py                                   # full table, N = 60000, about 20 minutes on two cores
    python table.py --N 6000                          # quick, noisier check
    python table.py --T 5 --B 2e-5                    # a single entry
    python table.py --eta 0.45 --K 0.3 --T 5 --B 2e-5 # other tilt or cutoff

Options: `--N` (sample points), `--eta`, `--K` (cutoff on |a0 k_x|), `--T`, `--B`, `--procs`.

## Parameters of the table

xi = +, eta_x = 0.3, v_y/v_x = 0.8, |a0 k_x| <= 0.3, window 1e-3 <= mu <= 3e-3. For graphene, v_x = 1e6 m/s gives hbar v_x/a0 = 4.63 eV, so 1 mK corresponds to 1.86e-8 in units of v_x/a0, and hbar/(e a0^2) = 3.26e4 T converts a0^2 e B to tesla.

## Method

- Each PLL band is monotonic in the window, so sigma, alpha and ell are summed band by band. Each band contributes |v(k_F)|/(2 pi) times a Fermi-function integral between its energies at the two window ends.
- delta_M and delta_WF are relative L2 norms, sqrt(sum (x - y)^2 / sum y^2), between the direct integral and the right-hand side of the Mott relation and of the Wiedemann-Franz law, evaluated at N points spread uniformly over the window. The Mott derivative is a central difference with step h = T/10.
- The Mott deviation sits near the steps of g(epsilon), whose width is about T. The derivative step must therefore be well below T, and N must be large enough for points to fall inside the thermal width around the steps. Small N gives a noisy delta_M (at N = 6000 it scatters between about 0.107 and 0.137 at 1 to 10 mK).

## Results (N = 60000)

| a0^2 e B | delta_M at 1, 2, 5, 10 mK | delta_WF at 1, 2, 5, 10 mK (x 1e-5) |
|---|---|---|
| 1e-5 | 0.134, 0.133, 0.135, 0.136 | 1.2, 1.6, 2.6, 3.7 |
| 2e-5 | 0.133 for all | 1.6, 2.3, 3.7, 5.2 |
| 4e-5 | 0.133 for all | 2.4, 3.3, 5.3, 7.5 |

## Variation with tilt and cutoff (full window, N = 60000, T = 5 mK, a0^2 e B = 2e-5)

| Case | delta_M | delta_WF | mean sigma~ |
|---|---|---|---|
| eta_x = 0.2, K = 0.3 | 0.133 | 6.9e-5 | 10.9 |
| eta_x = 0.3, K = 0.3 | 0.133 | 3.7e-5 | 36.5 |
| eta_x = 0.45, K = 0.3 | 0.135 | 2.0e-5 | 124 |
| eta_x = 0.3, K = 0.1 | 0.133 | 2.0e-4 | 2.84 |

Reducing the cutoff from 0.3 to 0.1 lowers the window-averaged sigma~ by a factor of 12.9.

## Caveats

- Sampling noise in delta_WF is about 2 to 5 per cent at the lowest temperatures, so quote two significant figures at most.
- Temperatures of 1 K and above fall outside the regime in which g is probed on a single step and are not part of the table.
- eta_x = 0 has no Fermi points in the window, and for eta_x = 0.1 some levels are not monotonic over the retained range of k_x, so the deviations are not evaluated there.
