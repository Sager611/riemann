REFUTATION: Lemma 6.4 with Lemma 6.2(v) (Sec. 6.2, far off-diagonal term of the new zero-density bound Thm 6.1) counts the far representatives (|gamma_r - gamma_s| > 1) on one side of each representative only, losing a factor 2 (copied from Bellotti arXiv 2405.12545 (4.13)); corrected, the paper's own certified pipeline gives q > 1 near log x = 875 (best 1.033 on [874.5, 875] after re-optimising the kernel), so Thm 1.1 (e^58 <= x <= e^3159) and N_R = 38,358,837,683 are not proved as published.

Record: Zenodo 10.5281/zenodo.23440770 v2 (Padhi, 2026-10-17), code github.com/jaideepsaipadhi/ramanujans-inequality.
How the gap e^103 < x < e^3158.442 is covered: Buthe's smoothed explicit formula for psi with RH to 3e12, a per-height
zero-free region (BTY inequality), a zero-repulsion lemma and a new log-free zero-density bound (Bellotti's method with
a sharp Barban-Vehov mean square). Needed accuracy: |pi - li| <= eps_R y/log y, eps_R(103) = 3.38e-9, eps_R(1000) = 3.38e-13,
eps_R(3158) = 3.37e-15 (3 L^4 eps_R -> 1). Binding range log x ~ 875 with published margin 0.9% (q = 0.99137).

Error and size (their Arb code rerun, all 3320 density cells; old values reproduced to 2e-14):
- far term is 9-23% of the density constant C; doubling it raises every cell's N+ by 9-23%.
- q on [872,877]: 0.99137 -> 1.0726 with their kernel, 1.0554 best kernel; [874.5,875] 1.0327; [880,880.5] 1.0413.
Also unjustified (minor): Lemma 6.2(iv) uses Bellotti's |G(1/2+it)| <= d5 W/log W, whose proof drops signs; rigorous bound
23.6 W/log^2 W is up to 2.06x larger (effect <= 5%). N3 relies on Bellotti's published tables, which share these flaws.

Repair (ours): |X^w - Y^w| <= (1 + (Y/X)^(2 nu)) X^Re(w) instead of 2 X^Re(w) shrinks C by 0.61-0.70. With the far term
doubled, this slack and the rigorous G bound, and without Bellotti's tables: max q 0.8957, transfer ratio 0.8949, all of
[57, 3159). So the theorem is plausibly true, but it rests on parts not re-derived here.

Checked and fine: eps_R values; Lemma 3.1 constant 0.80314 <= 0.80319; Buthe tail constant (arXiv 1410.7015 Prop. 3);
1/|zeta| <= 13.13 <= 14.57 and |zeta| <= 10.6 <= 11.89 on the box (B1); Gamma integrals 5.92164, 4.71363; Barban-Vehov constants
k1 = -0.6075, k0 = 0.3520 by own quadrature (paper -0.6071, 0.3503); residue algebra of Thm 6.5; pair LP; scan reproduced.
Not checked: zero-free-region table (Thm 4.6), repulsion table (Thm 5.2), Thm 6.5 remainders, BTY constants.
Small x: Johnston 2022 Thm 5.1 covers 38,358,837,683 <= x <= e^103 (Axler 2018 below e^43); Dudek-Platt/Galway give the last
counterexample 38,358,837,682. Large x: Axler 2022 Thm 1.1, x >= e^3158.442 (via Fiori-Kadiri-Swidinsky). Both confirmed.
New theorem? Not established: unconditional N_R <= e^3158.442 stays the best proved bound. Side finding: Bellotti
2405.12545 (4.13) and (4.15) have the same two flaws. Details: NOTE.md.

