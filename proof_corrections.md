# Proof corrections: Transport Phenomena, tp-2026-0129

The proof was typeset from an earlier version of the manuscript. Part A lists what must change so that it matches the current `arxiv_ss.tex`. Part B lists typos and wording slips (marked "also in the .tex" where the same error is in the source, so the .tex should be fixed too). Part C covers the publisher queries. Page numbers are the printed ones (PDF page 5 is p. 1). "Tex" means `arxiv_ss.tex`.

## A. Changes needed to match the current .tex

### Figures
1. Figure 3 (p. 8): replace the 3 × 2 panels (B̃ = ±10⁻⁴) with the three new single-column plots σ̃, α̃, ℓ̃ against μ̃ (B̃ = 2×10⁻⁴, η_x = 0.7, β̃ = 10³, 5×10³, 10⁴), as in the tex. Caption: "Strained tilted Dirac cone with ξ=+: σ̃_xx, α̃_xx, and ℓ̃_xx as functions of μ̃, normalised by τ. The ratio of the velocities is v_y/v_x=0.8." Delete the bracketed "[Parameters to be completed from the plotting scripts…]".
2. Figure 4 (p. 10, δσ for B̃ = ±0.5) is replaced by two figures: σ̃_xx against η_x (tex `sigma_etax1`, B̃=2×10⁻⁴, v_y/v_x=0.8, |a₀k_x| ≤ 0.023) and δσ against η_x on a symmetric logarithmic scale (tex `sigma_etax_delta`). Captions are in the tex.
3. Figure 5 (p. 11, σ̃ against B̃ for −0.02 ≤ B̃ ≤ 0.02) becomes σ̃_xx against B̃ for 0 ≤ B̃ ≤ 10⁻³ (tex `sigma_B`). The caption lists v_y/v_x=0.8 and 0 ≤ a₀k_x ≤ 10⁻³.
4. Figure 6 (p. 12, Mott/WF comparison, μ̃ = 0.02 to 0.05, B̃ = 10⁻⁴) is replaced by the new `mot1`/`mot2` plots (μ̃ = 0.03 to 0.048, B̃ = 2×10⁻⁴, η_x = 0.7, β̃ = 10⁴). Use the tex caption and remove the placeholder.
5. Figure 2 caption (p. 6): "The valley is ξ = + in all panels" is wrong. The panel labels show ξ = + in the left column and ξ = − in the right. Tex: "The valley is ξ=+ in the left column and ξ=− in the right column of each panel."
6. The proof will therefore have seven figures, not six. Renumber the cross-references: Fig. 4 (σ̃ against η_x), Fig. 5 (δσ), Fig. 6 (σ̃ against B̃), Fig. 7 (Mott/WF). Query Q4 (alt texts "Figures 1–6") must then cover Figures 1–7.

### Section 2.1 (p. 4)
7. Eq. (8) text: "this window corresponds to 0.07 T ≪ |B| ≪ 23 T in conventional units" becomes "this window reads 2×10⁻⁶ ≪ a₀² e|B| ≪ 7×10⁻⁴ in terms of the dimensionless field used in the plots".

### Section 2.2 (p. 5 and p. 7)
8. Eq. (12): replace the energy by E^ξ_{n,si}/v_x = η_x k_x + si √( 2|ζ|(n+1)(ξk_x r̃_ξ + s_{−ξ}s_ξ) − r̃_ξ² ζ²(n+1)² ). The index range becomes n ∈ {0, 1, 2, …, ⌊8(2+3ξa₀k_x)/(a₀²|ζ|(6+ξa₀k_x)²)⌋ − 1}.
9. Eq. (13): the small-B form is E^ξ_{n,si}/v_x ≃ η_x k_x + si √( 2|ζ|(n+1)(ξk_x r̃_ξ + s_{−ξ}s_ξ) ).
10. The sentence after Eq. (13) needs "⌊…⌋ − 1" as the upper bound. The condition for E_{n,−} < 0 becomes η_x a₀ k_x > √( 2|ζ|(n+1)(ξk_x r̃_ξ + s_{−ξ}s_ξ) ).
11. Eq. (14): E⁺_{n,si}/v_x ≃ η_x k_x + si √( 2|ζ|(n+1)(1 + 3a₀k_x/2) ).

### Section 3 (pp. 7 to 12)
12. Introduce Eq. (18) with "follows from differentiating Eq. (13). We define Λ_n ≡ 2|ζ|(n+1), w_ξ ≡ ξk_x r̃_ξ + s_{−ξ}s_ξ = 1 + 3ξa₀k_x/2, where we used s_ξ s_{−ξ} = 1 − a₀²k_x²/4. The derivative of w_ξ is ∂_{k_x}w_ξ = 3ξa₀/2." Then Eq. (18) becomes v^{ξ,si}_x(n;k_x) ≡ ∂_{k_x}E^ξ_{n,si} = v_x[η_x + si Λ_n ∂_{k_x}w_ξ/(2√(Λ_n w_ξ))] = v_x[η_x + si (3ξa₀/4)√(Λ_n/w_ξ)].
13. Replace "We note that the tilt-induced contribution … modifications of the bandstructure" with the tex sentences beginning "This expression is the derivative of the small-B form, and it is accurate when r̃_ξ²Λ_n ≪ 4w_ξ. The tilt enters the group-velocity as an additive term v_xη_x. The squared velocity …".
14. Delete the sentence "The response coefficients are asymmetric with respect to B = 0 …". Insert: "The spectrum in Eq. (12) depends on B only through |ζ|. All three coefficients are therefore even functions of B, and Fig. 3 shows B̃ > 0 only."
15. Unit paragraph: replace "Energies are measured in units of ℏv_x/a₀, which is about 4.6 eV … k_BT = (ℏv_x/a₀)/β̃" with the tex text: "Energies and temperatures are measured in units of v_x/a₀. The temperature then follows from T = (v_x/a₀)/β̃." Change "are now stated in its caption" to "are stated in its panel labels and caption".
16. Insert the two new paragraphs from the tex that begin "For Fig. 3, we take v_x=1.5×10⁻³ …" and "The integration over k_x is restricted to 0 ≤ a₀k_x ≤ a₀k_c …".
17. Replace the paragraph "The sensitivity on the values of η_x has been depicted via Figure 4 …" up to "… (which is the integration variable), decreases with rising |B|" with the tex paragraphs beginning "The dependence on η_x is depicted in Fig. 4" and "Furthermore, we also illustrate the dependence of the strength of σ_xx on B through Fig. 6". Those include the added sentence on oscillations at the lower temperature.
18. "The features that survive at η_x = 0, which include the asymmetry … reversing the strain gradient" becomes "The features that survive at η_x=0 are due to the strain. Since the spectrum depends on B only through |ζ|, the response is even in B, and reversing the strain gradient does not change it." Also change "∝η_x and ∝η_x² in Figure 4" to refer to the new σ̃ against η_x figure (Fig. 4).
19. Mott/WF paragraph (p. 11): replace "which shows that the relations produce a pretty good fit" with the tex sentences beginning "which shows that the relations capture the positions and the signs of the features of the direct results", up to "…the Mott relation stays close to 0.13".
20. In the g(ε) paragraph, add "in the regime of Table 1" after "it is the more robust of the two".
21. Replace the whole paragraph "We have quantified the deviations … to the accuracy stated here" (pp. 11 and 12) with the tex paragraph that begins "We have quantified the deviations for ξ=+, v_y/v_x=0.8 and η_x=0.3 in the window 10⁻³ ≤ μ̃ ≤ 3×10⁻³". It uses the updated values (δ_M = 0.133 to 0.136, and δ_M = 0.133, 0.133, 0.135 and δ_WF = 6.9×10⁻⁵, 3.7×10⁻⁵, 2.0×10⁻⁵ for η_x = 0.2, 0.3, 0.45).
22. Replace Table 1 with the tex table. It has different columns and a window of 10⁻³ ≤ μ̃ ≤ 3×10⁻³, and the temperatures are given in mK (1, 2, 5, 10). Use the tex caption.

### Section 4 (pp. 12 and 13)
23. Observables: the second item becomes "the dependence of the conductivity on the tilt, σ̃_xx of Fig. 4, which has a term linear in η_x for μ ≠ 0", with no "no such term at μ = 0". The third item becomes "the invariance of the response under the reversal of the sign of the pseudomagnetic field". The Mott estimate becomes "S ∼ π²T/(3eδ). For T/δ ∼ 0.1, this is a few tens of μV/K on restoring k_B/e ≃ 86 μV/K".
24. Scale estimates: copy the tex paragraph that begins "We estimate the scales as follows" (e²/(2π) form, "E_n ≃ 36 meV √(n|B|/(1 tesla))", "a₀² e|B| ≃ 2×10⁻⁵ … about 0.7 tesla"). The proof uses e²/h in places.

### Appendix (pp. 14 and 15)
25. Eq. (22): the middle term should read "… 2k_x r̃ + 2 s_− s_+ − ζ r̃² …". The proof has "−eB r̃²", and its "2υ s_− s_+" looks like a mis-keyed 2 s_− s_+ (the stray υ).
26. Eq. (26): α = [2 r̃²|ζ| + 2(√Δ − s_−s_+ − k_x r̃)] / (2|ζ|r̃²).
27. Eq. (27): −n = [2 r̃²|ζ| + 2(√Δ − s_−s_+ − k_x r̃)] / (2|ζ|r̃²) ⇒ √Δ = −r̃²|ζ|(n+1) + k_x r̃ + s_−s_+. Delete the second line of the proof's version.
28. Eq. (28): ε² = 2|ζ|(n+1)(k_x r̃ + s_−s_+) − r̃²ζ²(n+1)².
29. Replace "Imposing the requisite physical condition that ε² ≥ 0 …" with "that √Δ ≥ 0, which reads r̃²|ζ|(n+1) ≤ k_x r̃ + s_−s_+, we obtain an upper bound for n as ⌊…⌋ − 1".
30. Eq. (30): same corrected energy as Eq. (12). Eq. (31): √Δ = −r̃_ξ²|ζ|(n+1) + ξk_x r̃_ξ + s_ξ s_{−ξ}.

## B. Typos and wording

1. p. 2, "One of the remarkable outcome is" → "outcomes are". Also in the .tex.
2. p. 3, Eq. (3): "a_r" should be "a_{r_A}". Also in the .tex.
3. p. 3, text above Eq. (4): "r_a + M_α" → "r_A + M_α". Also in the .tex.
4. p. 3, Eq. (4): "t_α" on the right-hand side should be "t", and "γ_a" should be "γ_α". Also in the .tex.
5. p. 3, "ensure that the strain to be applied in such a way" → "that the strain be applied in such a way". Also in the .tex.
6. p. 5, "which we will call as the PLLs" → "which we will call the PLLs". Also in the .tex.
7. p. 5, "pseudo magnetic field" → "pseudomagnetic". Also in the .tex.
8. p. 5, "The most intriguing nature that is revealed comprises …" is awkward. Suggest "The most intriguing feature is that the PLLs form closed curves …". "Cranked up" is informal.
9. p. 7, "the PLLs obtained here are apart from being dispersive are tilted as well" → "the PLLs obtained here, apart from being dispersive, are tilted as well". Also in the .tex.
10. p. 7, "set-ups shown in Figure 3": Figure 3 shows plots, and the set-up is in Fig. 1(b). Also in the .tex.
11. p. 8, Eq. (17), α_xx: "∂_{ε_n}" should be "∂_{ε_n^{ξ,si}}" as in Eq. (16). Also in the .tex.
12. p. 14, "The above form is in the form of a Sturm–Liouville (SSL) problem" → "is a Sturm–Liouville (SL) problem"; "singular SSL scenario" → "singular SL". Also in the .tex (where it is spelled "Liuouville").
13. p. 14, "the interval (−∞, 0] (or, [0, ∞)) infinite" → "is infinite". Also in the .tex.
14. p. 15, Eq. (29): the exponent uses r_+ and the power uses r_ξ², which is inconsistent for the ξ = + solution. Use one symbol throughout. Also in the .tex.
15. p. 10 (old paragraph, "A nonzero μ always gvies a parabola") → "gives". This paragraph is replaced by item A.17, so the typo goes with it.
16. p. 12, Fig. 6 caption: "Comparision" → "Comparison". The tex caption replaces it.
17. Mixed US/UK spelling (the journal uses British English elsewhere): "quantization" (Landau-level-quantization, p. 2), "renormalizes", "visualization" (p. 4), "generalized" (p. 1, p. 15), "analyze" (p. 3), "behavior" (p. 7), "nearest-neighbor" (p. 1) against "nearest-neighbour" (p. 3). All are also in the .tex.
18. "Wiedemann–Franz" (en dash) and "Wiedemann-Franz" (hyphen) are mixed (pp. 11 and 12). Choose one. Also in the .tex.
19. Fig. 3 caption: "β" should be "β̃" ("3 different values of β"). The tex caption replaces it.

## C. Author queries and declarations

- Q1: the running head "Linear response from tilted Dirac cones" is 39 characters, within the 75-character limit.
- Q2 and Q3: author names and affiliation as printed are correct (Sanskar Sharma; Ipsita Mandal, corresponding).
- Q4 and Q5: alt texts and a main title are needed for every figure. After A.1 to A.6 there are seven figures. A main title for Figure 1, for example "Tilted Dirac cone and transport set-up", is still to be supplied.
- Q6: issue numbers are needed for refs [1]–[21], [23], [25]–[27]. These have not been looked up.
- Q7: publisher locations. [24] Ziman, Electrons and Phonons, Oxford University Press: Oxford. [22] Arovas, CreateSpace Independent Publishing Platform: check the imprint page of the book (Scotts Valley, CA or North Charleston, SC are both seen in listings).
- Please confirm the declarations on p. 13. "Data availability: Not applicable" differs from the tex statement "All data generated or analysed during this study are included in this published article". "No Artificial Intelligence (AI) tools were used" should be checked against how the manuscript and code were prepared.
