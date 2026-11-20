/-!
# Round 532 — node N11 (GRH families and Siegel zeros), round 26, GRH harvest triage of 1 Jan – 30 Jun 2025 (arXiv API: 18 unseen ids, 5 read): two partial transfers from the Lean-checked half-plane (Farey fractions with k-free denominators; maximal sums of odd-order characters, Prop. Spsig) for a harvest agent; three blocked; three explorations, no transfer

* Node: N11 (26th round; R522 before). Round type: GRH API harvest (lesson 307) of the one half-year no harvest had covered (grhharv3: 1 Jul – 15 Sep 2026; grhharv4: 1 Jan – 30 Jun 2026; grhharv5: 1 Jul – 31 Dec 2025; R492/R493: 15 Sep – 20 Nov 2026). Query: abs GRH or generalized/generalised/extended Riemann hypothesis, submitted 1 Jan – 30 Jun 2025: 18 ids unseen in our files; 5 with a plausible Dirichlet-only GRH use read from TeX hit files (rg GRH|ERH|zero-free with line numbers).
* Offered input: L(s,χ) ≠ 0 for Re s > θ₁ = 0.8733458561067… for every Dirichlet χ (QRT PR 141, Lean), hence prime-counting errors x^{θ₁+ε} for q ≤ x^{small} and Lindelöf-type bounds right of θ₁; nothing on the line, nothing for non-abelian or modular L.
* (1) 2504.15658 (Brun's constant under GRH, B < 2.1594): B, partial: Thms 1–2 and Lemma 3 transfer only with z < x^{(1−θ₁)/2−ε} ≈ x^{0.0633} and x^{θ₁} for the √x terms: a much weaker finite bound, no new constant.
* (2) 2507.00228 (Farey fractions with k-free denominators): B, partial: Prop. mu's x^{1/2+ε} becomes x^{θ₁+ε}; thm2's GRH direction gives displacement ≪ Q^{θ₁+ε} (the converse stays blocked); beta2's second moment becomes ≪ Q^{2θ₁−2+ε} = Q^{−0.25330…+ε} in place of Q^{−1+ε}: a candidate unconditional power saving.
* (3) 2505.07651 (sharp bounds for maximal sums of odd-order characters): B, partial: Prop. Spsig (TeX l. 260–266) holds unconditionally for every η < 1 − θ₁ = 0.12665; Thm LMUpp and Thm oddOrdMag need GRH-quality control at Q = log q: blocked.
* (4) 2505.11226 (integral points in thin sets of type II): m ≥ 2 already unconditional; m = 1 needs GRH for non-abelian Dedekind zeta functions (TeX l. 2237–2249): blocked.
* (5) 2502.11987 (average least negative Hecke eigenvalue): Thm averagenf already unconditional (l. 1549); the GRH parts concern L(s,f) and L(s, Sym^m f): blocked.
* Explorations: Kapranov's tropical root counting (N6): |·| is archimedean, coefficient cancellation breaks the slope–root equality (z² + 10z + 1: predicted e^{±10}, actual 0.101, 9.899); Kalman–Bucy innovations (N1): the covariance cos((γ_i − γ_j)u) is a Gram matrix, PSD for any ordinates, blind to Re ρ; revenue equivalence (N3): needs an allocation environment absent for ζ: no transfer.
* Verdict: two B-items (2507.00228 beta2 and thm2's GRH direction; 2505.07651 Prop. Spsig) for a grhharv6 agent (referee each before COROLLARIES); nothing recorded yet.
-/

namespace Fayed.R532
end Fayed.R532
