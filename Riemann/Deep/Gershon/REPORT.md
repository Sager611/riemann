DISCOVERY: Gershon (Preprints 202604.1513 v1) Lemma 11.6 is false: C_1(50) - C_1(30) = 0.03318 (Arb-rigorous) where the lemma allows <= 9e-11; it carries the whole n >= 100 region of the RH claim, and that region's other routes are circular or contradicted (R(100) >= 0.3057 vs the claimed < 0.149).

# Audit: Gershon, "The De Bruijn-Newman Constant Is Zero" (Preprints.org, 22 Apr 2026)
Paper text: `.scratch/gershon/paper.txt` (via r.jina.ai; preprints.org blocks curl).

**Verdict: the global proof fails in Section 11, Region C1 (r >= 51, n >= 100) of Theorem 11.3.** That region rests on Prop. 11.12 (case n >= 100). The proposition moves a sum S, computed only at n = 30, to every n >= 100 using Lemma 11.6. Lemma 11.6 is false.

## Chain (only Section 11 carries the global claim)
1. gamma_k = int_0^oo Phi(u) u^{2k}/(2k)! du. D_r(n) = det(gamma_{n+i-j}) (Def. 11.1).
2. Thm 11.1 (Edrei-Schoenberg): D_r(n) > 0 for all r, n gives Xi in LP, so RH. **Correct.**
3. Thm 11.3 splits the (r, n) plane into regions. A: r <= 2, by Borell. **OK.** B: 3 <= r <= 50, by finite certificates plus Prop. 11.6. C2: n <= 99, by certificates plus a dominant-pole tail. **Plausible.** C1: n >= 100, by Prop. 11.12. **Fails.**
4. Not load-bearing: Parts I-III. These are the kernel log-concavity (Thm 1.1), Jensen d <= 22 / n <= 14 (Thm 8.1) and the Lorentz reformulation (Thm 9.1, conditional on [K]_m < 0). The brief's suspects (a) "log-concave implies LP" and (b) "finite Jensen checks extend to all degrees" are not used. After Thm 1.1 the paper itself says TP2 does not imply LP. Its example for that, e^{-t^4}, is wrong: Polya proved its transform has only real zeros. A correct example is K(t) = e^{-t^2}(1 + t^4/4) (`Riemann/Tree/R038.lean`).

## Failing step: Lemma 11.6, "spectral-gap factorisation"
- Statement: |C_s(n) - C_s^oo| <= K_s delta^n, with sum K_s < 2 and delta = 0.452, where C_s(n) = n^2 |log Theta_s(n)|. Prop. 11.12 (eq. 61) uses it to put R(n) = S/(a n) with S <= 19.41 from n = 30, then gets R(100) <= 0.149 < 1.
- **Refuted (Arb, `.scratch/gershon/arb_c1.py`, 400 bits, rigorous integral tails):** C_1(30) = 0.846747739437362382 and C_1(50) = 0.879930958891510116, so the difference is 0.0331832194541477343 +- 7e-23. The lemma allows at most 2(delta^30 + delta^50) = 9.0e-11. The paper says it checked "C_s(30) vs C_s(50): discrepancies < 1e-9". In fact C_1(n) approaches its limit polynomially, as the paper's own Lemma 11.2 (Theta_1 = 1 - 1/n^2 + O(1/n^3)) already implies.
- **Its proof is void:** it writes the Taylor coefficients of the entire g = sum gamma_k z^k as gamma_m = R_1 rho_1^m + R_2 rho_2^m + O(delta_3^m rho_1^m). That is the expansion of 1/g. For an entire series it forces R_1 = 0, so the bound "K_s = O(|R_2/R_1|^s)" means nothing.
- **Independent check of eq. (61) at n = 100** (`arb_r100.py`). Telescoping (53) gives sum_s |log Theta_s(n)| >= mu_1 - mu_oo. The paper's own Cor. 11.2 gives mu_oo = 2 log(t_101/t_100). With mu_1(100) = 0.0151295485, t_100 = 236.5242297 and t_101 = 237.7698205 (Arb zeta zeros), R(100) >= 0.30567 against the claimed < 0.149. The true S(100) is at least 45.9, not 19.41.
- **Root cause:** R(n) is at least 1 - mu_oo/mu_1, which depends on the size of the gap between consecutive zeros t_n and t_{n+1}. No fixed S can bound it. R(n) < 1 for all n is RH-strength in itself.
- The "analytical alternative" (Cor. 11.2, Remark 11.8) assumes the zeros -t_m^2 of g are real and ordered, which is RH. The paper admits this in Remark 11.5. It also proves only "D_r(n) > 0 for all r implies unitarity", which is the converse of what C1 needs. The circle is: C1 needs unitarity, and unitarity is derived from C1.

## Lean (`Riemann/Deep/Gershon/Paper.lean`; builds; no sorry; axioms are only propext, Classical.choice, Quot.sound)
- Real definitions: `Phi`, `gam` (gamma_k), `D`, `L1`, `eps`, `Theta1`, `C1`.
- `gershon_RH`: regions A, B, C2, C1 plus Edrei-Schoenberg give `RiemannHypothesis`.
- `Lemma116` is the failing lemma at s = 1. `Witness` is the Arb fact, stated as a cited Prop. `lemma116_false : Witness -> K <= 2 -> not Lemma116 gam Cinf K 0.452`.
- `expansion_forces_R1_zero`: an entire series cannot have the geometric coefficient expansion with R_1 != 0 (fully proved).
- `claim61_false`: telescoping plus mu_r -> mu_oo < (4/5) mu_1 refutes the drag bound (61) with R = 0.149.
- `gershon_gap` bundles these results.

## Novelty (web_search; no contact with anyone)
- BitConcepts/riemann-solver, docs/LANDSCAPE.md, flags only Part I: the eq. (16) perturbation bound goes the wrong direction, and the e^{-t^4} example is wrong. It says nothing on Section 11.
- Michalowski (arXiv 2602.20313) studies the PF order of the kernel, not this argument.
- Not found anywhere: the refutation of Lemma 11.6 or of eq. (61). **New, pending reviewer judgement.**

## Files
`.scratch/gershon/`: `paper.txt`, `abstract.txt`, `arb_c1.py` and `arb_c1.log` (about 2 s), `arb_r100.py` and `arb_r100.log` (about 1 s), `gammas.py` and `c1.py` (mpmath cross-check). Arb runs with `.scratch/weil/venv/bin/python`. Peak memory was under 60 MB.
