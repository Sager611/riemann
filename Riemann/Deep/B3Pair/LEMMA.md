# Lemma C′ — B3Nu's Lemma C with an inner exponent pair from the model-phase class (b3pair, 20 Nov 2026)
Setting: B3Nu REPORT (B3N) §3 (Lemma S), §4 (Lemma C); B3Pair REPORT §1.1–1.2 (phases f_u, F_u, g_h, constants A, Λ_u, T_h, G_h). θ ∈ (1/3, 1), 0 < ν < θ/2, ε > 0 fixed after (θ, ν).
**Definitions (TTY = arXiv 2501.16779 v2).** Model phase function (Def phase-def, l. 196–202): smooth F: [1, 2] → ℝ (variable) with F^{(p+1)}(w) − (d/dw)^p w^{−ρ} = o(1) for all fixed p ≥ 0, some fixed ρ > 0.
Exponent pair (Def exp-pair-def, l. 309–320): fixed (k, l), 0 ≤ k ≤ ½ ≤ l ≤ 1, k + l ≤ 1, with Σ_{n∈I}e(TF(n/N)) ≪ (T/N)^{k+o(1)}N^{l+o(1)} for all model F, T ≥ N ≥ 1, intervals I ⊂ [N, 2N].
Non-asymptotic form (TTY l. 322–330, proved as Lemma beta-asymp l. 271–289): for all fixed ε, ρ > 0 there are δ > 0, P, C with |Σ_{n∈I}e(TF(n/N))| ≤ C(T/N)^{k+ε}N^{l+ε} whenever T ≥ N ≥ 1,
I ⊂ [N, 2N] and |F^{(p+1)}(w) − (d/dw)^p w^{−ρ}| ≤ δ for 0 ≤ p ≤ P, w ∈ [1, 2]. (TTY's printed (ii) omits "δ, ρ"; the quantifiers are those of Lemma beta-asymp.)

**Lemma C′.** Let (k′, l′) be an exponent pair (TTY sense) and (k, l) := A(k′, l′) = (k′/(2k′+2), (k′+l′+1)/(2k′+2)). For every fixed ε ∈ (0, θ − 2ν) (review 1 repair: was ε > 0), uniformly in g, coprime p, q (gp, gq ≤ y),
σ = ±1, the twist sign, the shifts, N ≤ N0T^ε, u ≤ √(2N) and every interval J′ ⊂ J_u ⊂ (V, 2V] (V = N/u):  |Σ_{v∈J′}e(F_u(v))| ≪_ε T^{3ε} Λ_u^k V^l  (T ≥ T0(ε, θ, ν)).
Consequently B3N Lemma C (statement and proof (a)–(d)) holds for this (k, l): its step (c) used only "|Σ_{J′}e(F_u)| ≪ T^{Cε}(L_u^kV_u^l + L_u^{−1})" (Λ_u = L_u/√(2π) ≥ 1).
**Proof.** Step 0. Λ_u ≥ 1 for T ≥ T0 (L_u ≥ T^{θ−2ν−ε}, B3N row 8). If V ≤ V0(ε): |Σ| ≤ V + 1 ≪ 1 ≤ Λ_u^kV^l. Let V > V0.
Step 1 (Weyl–van der Corput, GK Lemma 2.5). For every integer 1 ≤ Q ≤ V: |Σ_{J′}e(F_u)|² ≤ (1 + |J′|/Q)Σ_{|h|<Q}(1 − |h|/Q)S_h, S_h := Σ_{v,v+h∈J′}e(F_u(v + h) − F_u(v)), |S_{−h}| = |S_h|, |S_0| ≤ V + 1.
  By REPORT (1.2), for h ≥ 1: S_h = e(λ_uh − uh/2pq)·Σ_{v∈I_h}e(g_h(v)), I_h := J′ ∩ (J′ − h) ⊂ (V, 2V] an interval; the twist λ_u has left the phase.
Step 2 (the inner pair on g_h). Fix ε; take δ, P, C from the non-asymptotic form for (k′, l′) and ρ = 3/2. By (∗) of REPORT §1.2 (G_h^{(p+1)}(w) − (d/dw)^p w^{−3/2} = O_p(h/V + t),
  t ≤ 3T^{ε−2θ} as v + h ≤ 3V), for h < Q ≤ V^{2/3} the bound |G_h^{(p+1)} − (d/dw)^p w^{−3/2}| ≤ δ holds for p ≤ P once V ≥ V0 and T ≥ T0. Write g_h(v) = −sgn(A)T_hG_h(v/V), T_h = hΛ_u/2.
  (i) T_h ≥ V: |S_h| ≤ C(T_h/V)^{k′+ε}V^{l′+ε} ≤ CT^{2ε}(T_h/V)^{k′}V^{l′} (conjugate if sgn A = +1; T_h/V ≤ Λ_u ≤ √(Tu) ≤ T, V ≤ T; review 1 repair: the (T_h/V)^ε loss kept).
  (ii) T_h ≤ V/2: |g_h′| = (T_h/V)|G_h′| ∈ [0.3T_h/V, 0.51] on I_h (G_h′ = w^{−3/2}(1 + O(δ))), g_h′ monotone (g_h″ ≠ 0): Kusmin–Landau (GK Thm 2.1): |S_h| ≪ V/T_h.
  (iii) V/2 < T_h < V: |g_h″| ≍ T_h/V² ≍ 1/V: second-derivative test (GK Thm 2.2): |S_h| ≪ V^{1/2} + V^{1/2} ≪ (T_h/V)^{k′}V^{l′} (T_h/V > ½, l′ ≥ ½).
  So |S_h| ≪_ε T^{2ε}[(T_h/V)^{k′}V^{l′} + V/T_h] for 1 ≤ h < Q ≤ V^{2/3}.
Step 3 (sum over h, choice of Q). Σ_{1≤h<Q}|S_h| ≪ T^{2ε}[(Λ_u/V)^{k′}Q^{1+k′}V^{l′} + (V/Λ_u)(1 + log Q)]; with 1 + |J′|/Q ≤ 2V/Q (J′ = (a, b], b − a ≤ V):
  |Σ|² ≪ T^{2ε}[V²/Q + V^{1+l′−k′}Λ_u^{k′}Q^{k′} + V²Q^{−1}Λ_u^{−1}log V].  Q0 := V(Λ_u^{k′}V^{l′})^{−1/(1+k′)}.
  Λ_u^kV^l = V/√Q0 (as V²/Q0 = Λ_u^{2k}V^{2l}, below). If Q0 < 2: |Σ| ≤ V + 1 ≤ 2V ≤ 2√2·Λ_u^kV^l, the trivial bound (review 1 pass 2 repair: was Q0 < 1 with real Q = Q0).
  Else Q := ⌊Q0⌋ ∈ [Q0/2, Q0] ⊂ [1, V^{1−l′/(1+k′)}] ⊂ [1, V^{2/3}] (l′ ≥ ½ ≥ k′; an integer, as GK Lemma 2.5 needs), so Step 2 covers every h < Q and η = h/V ≤ V^{−1/3};
  V²/Q ≤ 2V²/Q0, V^{1+l′−k′}Λ_u^{k′}Q^{k′} ≤ V^{1+l′−k′}Λ_u^{k′}Q0^{k′} = V²/Q0, and the third term ≤ 2(V²/Q0)(1 + log V) (Λ_u ≥ 1): |Σ|² ≪ T^{3ε}V²/Q0 = T^{3ε}Λ_u^{k′/(1+k′)}V^{(1+k′+l′)/(1+k′)} = T^{3ε}Λ_u^{2k}V^{2l}, so |Σ| ≪ T^{3ε}Λ_u^kV^l. ∎
**Remarks.** (1) Steps 1–3 are the A-process (TTY Lemma exp-process(A), cited there from Ivić Lemmas 2.8–2.9) run inside the model class; the only class-specific facts are (∗) and
  η ≤ V^{−1/3} (GK Ch. 3 states the A-process for F(N, P, s, y, ε) with the same h ≤ εN restriction). (2) F_u is itself model (ρ = ½; REPORT §1.3(c)), so TTY's exponent pairs also
  apply to F_u directly; the A-route places the imported estimates on ρ = 3/2 phases. (3) ε-pairs: TTY pairs carry o(1) losses by definition; Lemma C(d) needs ν < ν1(θ) strictly, so T^{Cε} is absorbed.

## Imported estimates and their ranges (lesson 344; review 2 repairs applied)
| # | import (source, lines) | used for | needs (hypotheses as stated in the source) | holds here |
|---|---|---|---|---|
| 1 | TTY Def phase-def (l. 196–202), Def exp-pair-def (l. 309–320); non-asymptotic form (l. 322–330, quantifiers of Lemma beta-asymp l. 271–289) | δ, P, C | fixed ε, ρ > 0 ⇒ ∃ δ, P, C: T ≥ N ≥ 1, interval I ⊂ [N, 2N], F smooth on [1, 2], |F^{(p+1)} − (d/dw)^p w^{−ρ}| ≤ δ (p ≤ P) ⇒ |Σ_I e(TF(n/N))| ≤ C(T/N)^{k+ε}N^{l+ε} | ρ = 3/2, N = V, T = T_h ≥ V, I = I_h ⊂ (V, 2V], (∗) with h/V ≤ V^{−1/3}, t ≤ 3T^{ε−2θ}: δ-close for V ≥ V0, T ≥ T0 ✓ |
| 2 | TTY Thm new-exp-pair (l. 487–499): (89/1282, 997/1282) | Step 2(i) | as row 1 | ✓ |
| 3 | its proof: Lemma beta-duality (l. 362–368), Remark (l. 384–387), Table beta-table (l. 437–482); A, B (Ivić Lemmas 2.8–2.9), D (Sargos) processes (TTY l. 336–360), β-symmetry (l. 300–306) | pair from the table | line ≥ table on [0, ½]; [½, 1] by symmetry when l − k ≥ ½ | l − k = 908/1282 ✓; re-run exactly (logs/b3p_ttycheck.log): VALID; load-bearing rows below |
| 4 | Bourgain 2017 (3.18), Thm 4 = (3.19) (TTY rows on [1673/4038, ½]) | β on [0.414, 0.5] | only (3.1): min(|F″|, |F‴|, |F⁗|) > c1 (Bou l. 744–746; TY 2306.05599 l. 918–927) | GK class ⇒ (3.1) (TY l. 920–927) ✓ |
| 5 | D(13/84, 55/84) row = TTY Lemma D-process (l. 349–356; proof: Sargos 1995 Thm 7.1) applied to Bourgain's Thm 6 pair; TY l. 851–855 (same application, for F(N, P, σ, y, c)) | β on [0.394, 0.405) (load-bearing: 2.0·10⁻³ without) | (i) the input (13/84 + ε, 55/84 + ε) is an exponent pair for the class: Bourgain Thm 6 via (3.19) on [0, ½] from Thm 4 (α ≥ 17/42; (3.1) only), (4.1) = Huxley's Table 17.1 entry (4 + 103α)/128, (4.2) = classical (1/9, 13/18), α > ½ and subintervals by Bou §5; (ii) TTY's statement: "If (k, l) is an exponent pair, then β(α) ≤ max(k1 + α(l1 − k1), 1/12 + 2α/3) for 0 ≤ α ≤ 1", D(k, l) = ((5k + l + 2)/(8(5k + 3l + 2)), (29k + 21l + 10)/(8(5k + 3l + 2))) | (i): (4.1)'s hypotheses = Huxley Thm 17.1.4/17.4.2, verified for every σ > 0 by TY l. 929–966; (4.2) GK class ✓; (ii): D(13/84, 55/84) = (18/199, 593/796) exactly, and on [1508/3825, 62831/155153] the max is the D-line (1/12 + 2α/3 lower by ≥ 0.0022; logs/b3p_dcheck.log) ✓. Inherited, not re-read (Proc. LMS 1995): Sargos's own hypotheses of Thm 7.1, taken from TTY's model-phase statement and TY's GK-class use (review 2 pass 2 repair) |
| 6 | Trudgian–Yang A-row: A(4742/38463, 35731/51284), TY Lemma new_exponent_pairs_lem (TY l. 99–107; proof l. 836–985) | β on [0.2988, 0.3066) (load-bearing: 2.1·10⁻⁴ without, logs/b3p_ttycheck.log) | TY prove it for f ∈ F(N, P, σ, y, c), every fixed σ > 0, from the same rows (Huxley 17.1/19.2 with verified hypotheses, Bourgain (3.13)/Thm 4, D(Bourgain), classical pairs) | GK class ✓ (inherited: TY's computer step, not re-run) |
| 7 | Huxley 1996 Tables 17.1, 19.2 (TTY rows *), incl. the tangent row (89 + 908α)/1282 at α = 1/3 | β on its ranges | Thm 17.1.4, 17.4.2 (R = 1, 2), Lemma 19.2.1 | verified for the GK class by TY l. 929–970 ✓ (book not re-read) |
| 8 | Heath-Brown 2017 Thm 1, k = 5 (TTY Lemma heath-brown-2017, l. 410–425) | β on [0, 890/3277) | f smooth, 0 < λ5 ≤ f^{(5)} ≤ Aλ5 | after conjugation λ5 ≍ T_h/V⁵, A fixed ✓ |
| 9 | Weyl–van der Corput (GK Lemma 2.5) | Step 1 | integer Q ≥ 1, (a, b] of length ≤ V | Q = ⌊Q0⌋ (Q0 ≥ 2) ✓; Q0 < 2: trivial bound |
| 10 | Kusmin–Landau (GK Thm 2.1) | Step 2(ii) | F′ monotone, ‖F′‖ ≥ λ > 0 | r := T_h/V ≤ ½: 0.3r ≤ |g_h′| ≤ 0.51 ⇒ ‖g_h′‖ ≥ min(0.3r, 0.49) = 0.3r; |g_h″| > 0 by (∗) ✓ |
| 11 | second-derivative test (GK Thm 2.2) | Step 2(iii) | λ ≤ |F″| ≤ αλ, α fixed | V/2 < T_h < V: λ ≍ 1/V by (∗) ✓ |
| 12 | φ_σ(t) = 2σ√(2πt) − πt + t^{3/2}ψ_σ(√t) (B3N §4) | §1.1, (∗) | t small, ψ_σ analytic | v + h ≤ 3V ⇒ t ≤ 3T^{ε−2θ} ✓ |
| 13 | L_u ≥ T^{θ−2ν−ε} (B3N §5.1 row 8) | Λ_u ≥ 1, Step 0 | ν < θ/2, ε < θ − 2ν | ✓ (ε ∈ (0, θ − 2ν) in Lemma C′) |
| 14 | B3N Lemma C (a), (b), (d) and partial summation; B3N §5.1 rows 5–9 | unchanged | as audited there (B3N §5.1, §6) | ν/θ hypotheses re-checked at the new ν0: section 'Range check at the new ν0', R1–R18 ✓ |
| 15 | B3N Thm E1″ chain, §5.1 rows 1–25 (all ν/θ rows re-checked: R1–R18 below) | theorem at ν0 = ν1(θ_c) | Lemma C(d) needs ν < ν1(θ) strictly; row 19: ν0 < ν1(θ′), θ′ = θ − η > θ_c; ν1 < θ/2 on (1/3, 1); ν1 increasing | (7478θ − 2457)/10042: < θ/2 ⟺ 2457θ < 2457 ✓; slope > 0 ✓; ν0 = ν1(θ_c) < ν1(θ′) strictly ✓; Lemma C(c), (d) sums need 1 + k − l > 0 and 0 < 1 + l − k < 2: 1 + l − k = 5021/2742 ✓ |
| 16 | uniformity (TTY Prop auto, l. 180–191, or row 1) | all bounds | constants depend on (ε, k′, l′, ρ) and ψ_σ only | ✓ |

## Range check at the new ν0 (independent-referee repair, 20 Nov 2026; logs/b3p_rangecheck.log, exact rationals, script source printed in the log)
Data: θ_c = 4939/10000, ν0 = ν1(θ_c) = 6181921/50210000 ≈ 0.1231213 (ν1(θ) = (7478θ − 2457)/10042). The edge θ = θ_c is excluded (theorem: θ > θ_c). Theorem at θ ∈ (θ_c, 1): η = (θ − θ_c)/2, θ′ = θ − η ∈ (θ_c, θ); the moment (E1″)
is used at (θ′, ν0), the counting and transfer at θ. Order of choices (review repair): fix θ > θ_c, then θ′, δ; ε < min{θ_c − 2ν0, −E/C} (C the fixed loss constant of
Lemma C′, E < 0 from R7); then A > A_min, σ > (1 + A)/(1 − 2ν0), D > D_min; all finite and T-independent for fixed θ (no uniformity as θ ↓ θ_c is claimed or needed).
R1–R6, R8, R9, R14–R16: linear in θ′ (or θ) at ν = ν0 with positive margin at θ_c, hence on (θ_c, 1); the upper bound θ′ < 1 (R2, R12) is open and holds as θ′ < θ < 1;
R7: zero at θ_c, strictly increasing, so strict for θ′ > θ_c; R10, R11, R13: parameter choices, finite for fixed θ.
| # | hypothesis | source lines | at (θ_c, ν0), exact | for θ′ ∈ (θ_c, θ) |
|---|---|---|---|---|
| R1 | ν < ½ (Cor B′; Young's diagonal; X0 ≥ 2√(2π)T^{1/2−ν} → ∞) | B3F l. 55, 85; B3F l. 5, 17 (B3N row 4); E1 l. 41 | ½ − ν0 = 18923079/50210000 | ν-only ✓ |
| R2 | 1/3 < θ′ < 1 (Lemma S, Lemma C setting) | B3F l. 55; B3N §3, §4 | θ_c − 1/3 = 4817/30000 | increasing ✓ |
| R3 | θ′ ≥ 1/3 (off-band tail e^{−H²/400} ≤ e^{−T^{2/3}/400}) | B3A l. 63 | 4817/30000 | ✓ |
| R4 | ε < θ′ (c_r ≤ 2√(2π)T^{ε/2} ≤ H/46) | B3F l. 57; B3A l. 31 | ε ∈ (0, θ_c − 2ν0) = (0, 12434877/50210000 ≈ 0.24766) | θ′ − ε > 2ν0 ✓ |
| R5 | 0 < ν < θ′/2, ε < θ′ − 2ν (L_u ≥ T^{θ′−2ν−ε} ≥ 1; pq ≤ y² < H) | B3N l. 40, 79 (Lemma S, Lemma C: 0 < ν < θ/2), §4(c), §5.1 row 8, standing facts | ν0 = 6181921/50210000 > 0 (second-referee repair); θ_c/2 − ν0 = 12434877/100420000 | ≥ that ✓; ε as R4 ✓ |
| R6 | ν < (1 + θ′)/4 (side terms T^{2ν−(1+θ′)/2}) | B3N §4(d), row 9 | (1 + θ_c)/4 − ν0 = 10056207/40168000 | ✓ |
| R7 | ν < ν1(θ′) strictly (Lemma C(d) exponent E = (1 + l − k)ν + (k + l)/2 − (½ + l)θ′ < 0) | B3N §4(d); Lemma C′ | ν1 − ν0 = 0 at θ_c (edge) | ν1(θ′) − ν0 = (7478/10042)(θ′ − θ_c) = 3739(θ − θ_c)/10042 > 0; E = −(3739/2742)(θ′ − θ_c): fix Cε below it ✓ |
| R8 | ν < θ′ (N0 < T) | B3F l. 154 | θ_c − ν0 = 9308399/25105000 | ✓ |
| R9 | 2ν < 1 (pq ≤ T, N0 ≤ T²: Lemma B(ii) tail) | B3F l. 19, 52–63 (B3N row 6) | 1 − 2ν0 = 18923079/25105000 | ν-only ✓ |
| R10 | A(1 − θ′) > 1 + ν − θ′ + δ (off-band, representation error, ε = +1 term) | B3F l. 91; B3A l. 63; E1A l. 61, 76 | A_min(θ_c) = 31593202/25411281 ≈ 1.2433 (δ = 0) | A_min(θ′) = (1 + ν0 − θ′ + δ)/(1 − θ′) finite for θ′ < 1: fix A ✓ |
| R11 | σ(1 − 2ν) > 1 + A (K-Bessel branch T^{1−σ(1−2ν)} ≪ T^{−A}) | B3F l. 98–99; E1A l. 77 | 1 − 2ν0 > 0 | fix σ > (1 + A)/(1 − 2ν0) ✓ |
| R12 | θ′ < 1 (residue change o(H/L); supp w ⊂ [T/2, 2T]; K = T^{1−θ′} → ∞) | B3F l. 84; B3A l. 54; B3N rows 7, 10, 11 | 1 − θ_c = 5061/10000 | θ′ < θ < 1 ✓ |
| R13 | windows: ν0 < ν1(θ′) (= R7) and D > (1 + 2ν − θ + δ)/(2η) (second-referee repair: B3N l. 147, §5.1 row 19, stronger than B3F l. 110 / E1A l. 79's 1 + ν − θ + δ) | B3N l. 147; B3F l. 108–110; E1A l. 79 | — | D_min = (1 + 2ν0 − θ + δ)/(θ − θ_c) (2η = θ − θ_c > 0; edge excluded); δ = 0: 37770102/5021 ≈ 7522.426210 at θ_c + 10⁻⁴, 3777507279/5021 ≈ 752341.620992 at θ_c + 10⁻⁶ (logs/b3p_rangecheck3.log); numerator > 0, finite: fix D ✓ |
| R14 | θ > 0 (N(T, H) ~ HL/2π, Backlund–Jensen O(L)) | B3F l. 112; B3A l. 84; E1 l. 64 | θ_c > 0 | ✓ |
| R15 | 0 < λ < θ < 1 (Wang local pair, λ ↑ θ) | B3F l. 113; B3A l. 87; E1A l. 66 | ✓ | ✓ |
| R16 | 0 < θ < 1 (c(θ) defined, c′(θ) = ½cot²(θ/√2) > 0) | B3F l. 113, 177; E1 l. 83 | ✓ | ✓ |
| R16b | κ < 1 and c(θ) < 1 (Thm 5.1; c(θ) = 2 − θ/2 − cot(θ/√2)/√2) | E1 l. 67 | κ ≈ 0.0861981; c(θ_c) ≈ −0.1886575 (Arb, logs/b3p_rangecheck2.log) | c increasing, c(1) ≈ 0.6725007 < 1 ✓ |
| R17 | κ(ν0, 44/5) > k*(θ_c); k* decreasing ⇒ κ > k*(θ) | B3F l. 114, 179–180; logs/b3p_arb_tty_4000.log | κ − k* ≈ 2.6350934·10⁻⁷ (Arb ball, radius 7.35·10⁻²⁸) | ✓ |
| R18 | not used: SL l. 199 (ν < θ − ½), SL l. 536–538 (θ > ½), SL l. 731–732 (θ − η > ½), Taylor J > (1 + ν − θ + δ)/(2θ − 1) | B3F l. 110; B3A l. 74–80, 136–137; E1 l. 53; E1A l. 60, 75 | exact profile: no Taylor step | — |
Inventory (E1 l. 67 added after the scan's first pass): every inequality in θ or ν found by a pattern scan of B3F (20 lines), B3A (22), E1 (27), E1A (18) and B3N §§3–5 is a row above or one of B3N/B3F's replaced conditions
ν < (3θ − 1)/4 (B3F l. 102, 108–109; B3A l. 34, 80; E1 l. 11–12, 51, 55, 62, 87–89; E1A l. 74, 81, 103), now R7, or a remark (B3A l. 47; E1 l. 3, 9, 98–113: statements and literature).
Result (review repair): every row holds for every fixed θ ∈ (θ_c, 1) with θ′ = (θ + θ_c)/2 and ν = ν0; at the excluded edge θ′ = θ_c, R7 is an equality and R13 needs η > 0.
No number changes.

