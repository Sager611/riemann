STATUS (2026-10-23, orchestrator): NOT NEW. The same theorem (N_R = 38,358,837,683 unconditionally) was posted by liqsweep on 14 Oct 2026 by another route (OpenAI's 7/8 zero-free half-plane); see NOVELTY.md. This repair of Padhi's route is a second proof that avoids the half-plane, but its certification fixes are incomplete (reviewers' list: .scratch history; outward rounding, the E2 small-argument region, smoothing constants for c = 75, 95, 105, upstream scripts) and the cover must be extended past log x = 3160.3 (Axler 2022 relies on the faulty case of Dudek–Platt Lemma 2.1). Work stopped.

DISCOVERY: Ramanujan's inequality pi(x)^2 < (e x/log x) pi(x/e) holds for every real x >= 38,358,837,683, unconditionally; so N_R = 38,358,837,683 (best proved before: N_R <= e^3158.442, Axler 2022).

What: a repaired and fully re-certified version of Padhi's proof (Zenodo 10.5281/zenodo.23440770 v2, 2026-10-17, whose published
certificate fails: factor-2 undercount of the far term, see Riemann/Deep/Ramanujan/). Credit: Padhi's framework plus our corrections.
Coverage: x <= e^103 Johnston 2022 (Thm 5.1); e^58 <= x <= e^3159 this work; x >= e^3158.442 Axler 2022 (Thm 1.1).
Certified result (ball arithmetic, python-flint): every one of the 317 scan intervals on [57, 3159) has q < 1, max q = 0.92768 on
[847, 857] (published, flawed: 0.99137); transfer to pi - li: max ratio 0.926851 at log y = 857 (< 1 needed), coverage complete.
How q encodes eps_R: |psi(t) - t| <= q eps_R(L2+1)/(1 + 2/L1) t + a(t) on [L1, L2]; the transfer then gives
|pi(y) - li(y)| < eps_R(log y + 1) y/log y and < eps_R(log y) y/log y, eps_R ~ 1/(3 L^4), which is Ramanujan by Padhi's Lemma 2.1.
Corrections: far term counted on both sides of each zero; residue factor (1 + (Y/X)^{2 nu})/2 instead of 1; rigorous
|G(1/2 + i tau)| <= 16 prod(1 + (sqrt p - 1)/p^2) W/log^2 W; exact Gamma in the pair sums; c1 tail redone.
Re-derived and certified inputs (NOTE.md section 3):
- zero-free region (replaces Padhi Thm 4.6): own BTY-type iteration from 1/5.558691, 5 steps, 137,828 boxes, rerun identical;
  lambda(29) = 0.21494 (Padhi 0.21384; his c10 carries an extra factor kappa); E2, E3, G, C_w, I0, m_p, Lemma 5 all recomputed.
- repulsion (Padhi Thm 5.2): not used (d* = 0). Bellotti's density tables (N3): not used.
- density cells (Padhi Thm 6.1): all 3320 recomputed, log N+ down 0.42 on average (range -0.66 ... +0.04).
- Barban-Vehov mean square (Padhi Thm 6.5): k1 <= -0.59807, k0 <= 0.37066 (rigorous contour integral), A1 = 148.0, Ax = 7.77,
  AV = 10.74, AU = 10.35, AK = 17.66, remainders R1-R8 and S_sharp with own constants; (B1): |1/zeta| <= 14.29, |zeta| <= 10.59.
- Buthe explicit-formula constants (tail 0.16(x+1)/sinh c ..., |Theta| <= 2), smoothing constant 0.803185 < 0.810492,
  Gamma integrals 5.92164 <= 5.922 and 4.71363 <= 4.72, eps_R expansion (Lemma 2.2): rerun or re-derived.
- RH height 3e12: Platt-Trudgian, Bull. LMS 53 (2021).
Used as published (not re-derived): Mossinghoff-Trudgian-Yang 2024; BTY (arXiv 2603.21490) Lemmas 6 and 13; Bellotti's
log-free framework (arXiv 2405.12545, detector and Halasz-Montgomery step; its flawed (4.13), (4.15) replaced); Kadiri-Lumley-Ng
Table 1 (rows checked); Hasanalizade-Shen-Wong; Hiary-Patel-Yang, Patel-Yang, Hiary; Ramare (Moebius sums, sum mu^2/phi);
Rosser-Schoenfeld; Dusart; Montgomery-Vaughan; Buthe 2016/2018; Johnston 2022; Axler 2022.
Margin left: 7.2% (q), 7.3% (transfer); unused slack: crude box constant in R2 (A1box 1.7e8 -> 148 is valid), kernels not re-optimised.
Files: NOTE.md (derivations, quotes of the false step with line numbers), code/ (scripts), code/out/ (certified outputs, logs).
