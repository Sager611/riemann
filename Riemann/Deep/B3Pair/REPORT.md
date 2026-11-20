RECORD: θ > 4939/10000 = 0.49390 (inner pair Tao–Trudgian–Yang (89/1282, 997/1282), outer A(·) = (89/2742, 1184/1371); reviewer A on Lemma C′: CORRECT (pass 3); reviewer B on the range table: CORRECT (pass 3); inherits TTY's pair, i.e. Bourgain's Thm 6 (via Sargos's D-process) and Trudgian–Yang's pair for GK's class; Arb 4000 = 3000 bits identical, κ − k* = 2.635·10⁻⁷ > 0)

# B3Pair — Lemma C with model-phase exponent pairs (Bourgain; Tao–Trudgian–Yang): the short-window simple-zero threshold (node N8)
b3pair (Orca), 20 Nov 2026, 12:27– UTC. Notation: B3Nu REPORT (= B3N) §§3–5, B3Cert FULLPROOF (= B3F) §1. H = T^θ (window), y = T^ν, (p, q) = 1, gp, gq ≤ y, N0 = pqT/H², ε > 0 small.
Sources (read today, TeX): TTY = arXiv 2501.16779 v2 (Tao_Trudgian_Yang_v2.tex), Bou = arXiv 1408.5794 (zeta1.tex). Code code/, logs logs/ (memcap peaks in the logs).
To avoid clashes: Lemma S's sign σ ∈ {±1}; TTY's model exponent is written ρ; the Weyl differencing length is Q (B3N's "H ≤ ..." in the brief), h its variable.

## 1. Falsification first (task 1): the phases of Lemma C and their class
**1.1 The undifferenced phase.** B3N Lemma S + Lemma C(c): for σ = ±1, u ≤ √(2N), dyadic N ≤ N0T^ε, v ∈ J′ ⊂ J_u ⊂ (V, 2V], V = V_u = N/u ≥ √(N/2):
  F_u(v) = f_u(v) + λ_u v, f_u(v) = Φ_σ(uv/pq)/(2π) = (T/2π)φ_σ(cv), c := u/(pqT), λ_u := ±(p̄u mod q)/q ∈ [0, 1) (twist sign ±; p = 1 ⇒ p̄ = 1; q = 1 ⇒ λ_u = 0),
  φ_σ(t) = log(1 + v_σ) + v_σ/(1 + v_σ), v_σ² = 2πt(1 + v_σ), sgn v_σ = σ; φ_σ(t) = 2σ√(2πt) − πt + t^{3/2}ψ_σ(√t), ψ_σ analytic at 0 (B3N §4, checked: v = σ√(2π)s + πs² + …, s = √t).
  Hence f_u(v) = A v^{1/2} − (u/2pq) v + E_u(v), A := σ√(2Tu/(πpq)), E_u^{(r)}(v) ≪_r t·|A|v^{1/2−r} (r ≥ 0), t := cv ≤ 2N/(pqT) ≤ 2T^{ε−2θ}.
  Parameters and ranges: g ≤ y; coprime p, q with gp, gq ≤ y; σ = ±1; twist sign ±; u ≤ √(2N); N ≤ N0T^ε; h < Q. Shifts (α, β, z) enter only the amplitude g_u(v) (B3N §4(c)),
  not the phase. Λ_u := |A|/(2√V) = |f_u′(V)|(1 + O(√t)) ≍ L_u (B3N's L_u = √(Tu/(pqV)); Λ_u = L_u/√(2π)), and L_u ≥ T^{θ−2ν−ε} (B3N row 8).
**1.2 The differenced phases.** For 1 ≤ h < Q and v, v + h ∈ J′: F_u(v + h) − F_u(v) = g_h(v) + const, const = λ_uh − uh/(2pq) (unimodular factor; the twist is gone),
  g_h(v) := A[(v + h)^{1/2} − v^{1/2}] + E_u(v + h) − E_u(v),  g_h^{(r)}(v) = A h (½)_{r+1} v^{−1/2−r}[1 + O_r(h/V) + O_r(t)]  (r ≥ 0; (½)_m := ½(½ − 1)⋯(½ − m + 1)).
  Normalized: T_h := |A|h/(4√V) = hΛ_u/2, G_h(w) := −sgn(A)·g_h(Vw)/T_h (w ∈ [1, 2]); then for every fixed p ≥ 0, uniformly in all parameters above,
  G_h^{(p+1)}(w) − (d/dw)^p w^{−3/2} = O_p(η + t),  η := h/V < Q/V.   (∗)
  So G_h is a TTY model phase function with ρ = 3/2 (TTY Definition phase-def, l. 196–202: F^{(p+1)}(u) − (d/du)^p u^{−ρ} = o(1) for all fixed p) iff η + t = o(1) uniformly;
  t ≤ 2T^{ε−2θ} always; η = o(1) needs Q = o(V). Sign: Σe(g_h) = Σe(−sgn(A)T_hG_h(v/V)) is Σe(T_hG_h(v/V)) or its conjugate.
  GK form: g_h′(x) = y_h x^{−3/2}(1 + O(η + t)) and likewise to every fixed order: g_h ∈ GK's monomial class F(V, P, 3/2, y_h, ε0) (GK (3.3.3)) for every P, ε0 once η + t ≤ c(P, ε0).
**1.3 Where the class fails (smallest cases) and why it does not matter.**
  (a) η not small: G_h is then the fixed profile [(w + η)^{1/2} − w^{1/2}]·(−4/η), not asymptotic to −2w^{−1/2}: at η = ½, w = 1: G_h′(1) = (√(2/3) − 1)/(−η/2)… = 0.7340 vs 1 (27 % off, all p).
      Smallest case: V bounded (e.g. N ≤ 8, u = 1, V ≤ 8, h = 1, η ≥ 1/8): no differencing length Q ≥ 1 gives η = o(1). Remedy (standard, §2 Step 3): bounded V ⇒ trivial bound
      |Σ| ≤ V ≪ 1 ≤ Λ^kV^l; unbounded V ⇒ the optimal Q0 = V·(Λ^{k′}V^{l′})^{−1/(1+k′)} ≤ V^{1−l′/(1+k′)} ≤ V^{2/3} (l′ ≥ ½, k′ ≤ ½), so η ≤ V^{−1/3} → 0 uniformly.
  (b) T_h < V (h < 2V/Λ_u): TTY's pair definition needs T ≥ N (TTY l. 313); there Kusmin–Landau (T_h ≤ V/2) and the second-derivative test (V/2 < T_h < V) give |Σ| ≪ V/T_h + V^{1/2}: GK's "+L^{−1}" term.
  (c) Twist: λ_u enters F_u but not g_h. F_u itself is also model (ρ = ½): its normalized first derivative carries √(2π)λ_u/L_u ≤ √(2π)T^{−(θ−2ν−ε)} and √(π/2)·√t, both o(1).
  (d) σ = ±1 and the twist sign: only sgn A and λ_u; conjugation. Shifts: amplitude only. g, p, q: only through c = u/(pqT) (t ≤ 2T^{ε−2θ}) and λ_u.
  Verdict (TTY sense; ε ∈ (0, θ − 2ν)): along every parameter sequence with V → ∞, every normalized g_h with h < Q0 is a model phase function (ρ = 3/2), uniformly in (g, p, q, u, N, σ, ±,
  shifts, h); bounded V is not claimed model and is handled by the trivial bound (review 1 repair). No failing case beyond (a)–(b), which the A-process step handles.

## 2. Lemma C′ (task 2): LEMMA.md
Statement: for every TTY exponent pair (k′, l′) and (k, l) = A(k′, l′), |Σ_{v∈J′}e(F_u(v))| ≪_ε T^εΛ_u^kV^l uniformly in all parameters of §1; so B3N Lemma C, hence Theorem E1″,
holds for (k, l) (ν1(θ) = [θ(1 + 2l) − (k + l)]/(2(1 + l − k))). Proof: the A-process inside the model class (Weyl–van der Corput; inner pair on G_h for T_h ≥ V; Kusmin–Landau and the second-
derivative test for T_h < V; Q0 = V(Λ_u^{k′}V^{l′})^{−1/(1+k′)} ≤ V^{2/3}, so η ≤ V^{−1/3}); imported estimates and ranges: §5. Reviews: §6.

## 3. The inner pair (task 3a): choice and provenance
Gain: ν1(θ) = (3θ − 1)/4 + g(1 − θ)/4, g = (1 − s)/(3 + s), s = 2k′ + l′ = 3·[k′ + (l′ − k′)/3]: s/3 is the pair's β-line at α = 1/3 (TTY Lemma beta-duality, l. 362–368: (k, l) is a pair
iff β(α) ≤ k + (l − k)α on [0, 1]). So min s = 3·(best β-bound at α = 1/3 carried by a pair line). TTY Table beta-table (l. 437–482) at α = 1/3: the row 89/1282 + (454/641)α (Huxley 1996
Table 17.1, range [423/1295, 227/601)), value 1175/3846 (code/b3p_ttycheck.py, logs/b3p_ttycheck.log). TTY's pair (89/1282, 997/1282) is exactly this line, so s = 1175/1282 is optimal among pairs
derivable from TTY's table; A-, B-images move away (s(A(89/1282, 997/1282)) = 0.9285; B-images have l + 3k > 1). Any (k, l) applies to F_u directly (§1.3(c)), but the best non-A line dominating
the table has gain (1 − l − 3k)/(1 + l − k) = 0.02036 < 107/5021 (float LP, logs/b3p_gain.log). Candidates (thresholds θ1: float bisection, code/b3p_thr.py = copy of b3nu_thr.py, 4000 bits, logs/b3p_thr.log):
| inner (k′, l′) | outer (k, l) = A(k′, l′) | s | g | θ1 | rests on |
|---|---|---|---|---|---|
| (13/238, 97/119) | (13/502, 445/502) | 110/119 = 0.924370 | 9/467 = 0.019272 | 0.493970073 (control = B3N) | GK classical |
| A(13/84, 55/84) = (13/194, 76/97) | (13/414, 359/414) | 89/97 = 0.917526 | 2/95 = 0.021053 | 0.493908737; +10⁻⁶: 0.493908745 | Bourgain Thm 6 |
| (652397/9713986, 7599781/9713986) | A(·) | 0.916676 | 0.021274 | 0.493901114 | TTY |
| **(89/1282, 997/1282)** | **(89/2742, 1184/1371)** | **1175/1282 = 0.916537** | **107/5021 = 0.021310** | **0.493899866; +10⁻⁶: 0.493899893** | TTY Thm new-exp-pair |
Chosen: (k, l) = (89/2742, 1184/1371): ν1(θ) = (7478θ − 2457)/10042 = (3θ − 1)/4 + 107(1 − θ)/20084; ν1(θ) < θ/2 ⟺ 2457θ < 2457 (θ < 1).
**Provenance audit (load-bearing input).** TTY prove the pair by Lemma beta-duality from their table ("further computer calculation"). Reproduced exactly (b3p_ttycheck.log): the line dominates the
table + classical pairs on [0, ½] (max(env − L) = 0, tangent on the Huxley row); TTY's Remark (l. 384–387) covers [½, 1] by the B-symmetry (l − k = 908/1282 ≥ ½). Without the rows from Bourgain
2017 (AD, D of (13/84, 55/84), (3.18), the pair itself) and Trudgian–Yang's A-row, the line is NOT dominated: short by 2.017·10⁻³ at α ≈ 0.3943 (on TTY's printed row ranges). Without only the Trudgian–Yang (TY, arXiv
2306.05599) A-row: short by 2.15·10⁻⁴ at α ≈ 0.3065. So the chosen pair needs Bourgain's Theorem 6 in GK's class (through Sargos's D-process; LEMMA table row 5) and TY's pair (row 6);
A²(Bourgain) needs Theorem 6 directly. Bourgain's caveat (l. 971–974) concerns only (3.19) below α = 17/42: there it uses (4.1) = Huxley [H1] Thm 3, the bound (4 + 103α)/128, and (4.2) =
the classical pair (1/9, 13/18), which "assume additional hypotheses concerning F, beyond condition (3.1)", checked there for F = log x; Theorem 4 (α ∈ [17/42, ½]) and (3.18) need only (3.1).
Those hypotheses are Huxley 1996 Thm 17.1.4/17.4.2 (the bound is a Table 17.1 entry) and GK's class; TY (proof of Lemma new_exponent_pairs_lem, TeX l. 918–970) verify both for F(N, P, σ, y, c),
every fixed σ > 0 (P ≥ 8, c < 1/(1000(σ + 1))). The conditions are σ-generic (ratios (σ + r)/(σ + r − 1), 3(σ + R + 1)/(σ + R + 2), a 3×3 determinant): no exponent is exceptional, and
g_h ∈ F(V, P, 3/2, y_h, c) for every P, c once T ≥ T0 by (∗). Inherited, not re-derived: Bourgain §5 (α > ½ by the B-process; subintervals by Sargos 1995 p. 310, used the same way by TY),
TY's computer step (TTY's was re-run here, exactly).

  Cross-check (ANTEDB = github.com/teorth/expdb, blueprint beta.tex, read today): its A-process for β (Lemma vdca-beta) uses exactly our normalization, "−(1/σ)(N/j)(F(u + j/N) − F(u)) is a
  model phase function", with j ≤ T^h, h < α (j = o(N)); at α = 0.4 its non-Bourgain bounds (Huxley 1993 Thm 1: 0.35667; Huxley 2005: 0.35647; Huxley 1993/Sargos (4 + 103α)/128: 0.35313;
  Huxley–Kolesnik 1991: 0.35455; Robert 2016: 0.36667) all exceed the pair line (89 + 908α)/1282 = 0.35273; only Bourgain's D(13/84, 55/84) row (0.35226) lies below. Same verdict: Bourgain-dependent.
  Fallback without Bourgain: A²(Huxley 2005 pair (32/205, 269/410)) = (8/253, 1755/2024), θ1 = 0.493919950 (+10⁻⁶: 0.493919959; logs/b3p_thr_hux.log), certified at θ_c = 0.49392 (§4).

## 4. Certificate (task 3b)
**Theorem (simple zeros, b3pair).** liminf_{T→∞} S(T, T^θ)/N(T, T^θ) ≥ h(θ_c; κ) > 0 for every fixed θ ∈ (θ_c, 1), θ_c := 4939/10000 = 0.49390 (before: 0.49398, B3Nu).
Proof: B3N §5 verbatim with (k, l) = A(89/1282, 997/1282) = (89/2742, 1184/1371) in Theorem E1″ (Lemma C′ in place of the classical pair in Lemma C(c)); ν0 := ν1(θ_c); for θ > θ_c,
η = (θ − θ_c)/2, θ′ = θ − η: ν0 < ν1(θ′) (ν1 increasing); TTY's o(1) losses are absorbed by this strict inequality (fix ε after θ′); detector = SL/EXP-010 frozen (P, Q), R = 44/5.
Arb (code/b3p_arb.py = copy of b3nu_arb.py, construction executed verbatim from code/b3c_arb_copy.py; check E1 κ1 = 0.0844779079697484032032447429603 ± 4.05·10⁻³²; 4000 and 3000 bits identical):
| outer pair | θ_c | ν0 = ν1(θ_c) | κ − k*(θ_c) | h(θ_c; κ) | logs |
|---|---|---|---|---|---|
| (89/2742, 1184/1371) | 4939/10000 | 6181921/50210000 = 0.1231213 | 2.63509·10⁻⁷ ± 7·10⁻²⁸ | 1.86875·10⁻⁷ | b3p_arb_tty_4000.log, _3000.log |
| A(89/1282 + 10⁻⁶, 997/1282 + 10⁻⁶) = (44500641/1371001282, 592000641/685500641) | 4939/10000 | 0.1231212 | 2.10405·10⁻⁷ | 1.49214·10⁻⁷ | b3p_arb_ttyd_4000.log, _3000.log |
| fallback A²(Huxley + 10⁻⁶) = (6400041/202400246, 175500205/202400246) | 3087/6250 = 0.49392 | 0.1230623 | 8.14791·10⁻⁸ | 5.77817·10⁻⁸ | b3p_arb_huxd_4000.log, _3000.log |
Exact checks printed: ν0 − (3θ_c − 1)/4 = (1 − θ_c)(1 − l − 3k)/(4(1 + l − k)), ν0 < θ_c/2. Control: Cor B′'s range at θ_c gives κ − k* = −1.48·10⁻³ < 0. Float bisection (b3p_thr.log)
reproduces B3N's 0.493970073 and R515's 0.493908737; R515's 0.493908812 is the outer pair + 10⁻⁶ (reproduced, logs/b3p_thr_r515.log), the inner + 10⁻⁶ gives 0.493908745. Footprint: max RSS 35 MB (/usr/bin/time -l, logs/b3p_arb_footprint.log, b3p_thr.log); memcap peaks ≤ 20 MB.


## 5. Ranges (task 4): LEMMA.md "Imported estimates and their ranges" (16 rows; reviewer B passes, §6)
Every ν/θ hypothesis of the B3N/B3F chain (B3F, B3A, E1, E1A, B3N §§3–5) re-checked in exact rationals at (θ_c, ν0) and for θ′ ∈ (θ_c, θ): LEMMA.md section 'Range check at
the new ν0', rows R1–R18 (logs/b3p_rangecheck.log; independent-referee repair: 'only three rows depend on the pair' did not establish the extension). All hold; no number changes.

## 6. Reviews (reviewer subagents; texts pasted inline; session logs under the pi-subagents run dirs)
- Reviewer A (Lemma C′, REPORT §1). Pass 1 GAP: (1) "every ε > 0" must be ε ∈ (0, θ − 2ν) (else Λ_u < 1, t ≍ 1: θ = ½, ν = 1/5, ε = 1, N = T); (2) the inner pair's (T_h/V)^ε loss was dropped;
  (3) §1.3 claimed model phases for bounded V. Pass 2 GAP: real Q = Q0 ∈ (1, 2) breaks Σ_{h<Q}1/h ≪ log Q (integer Q needed). All repaired (marked in the text). Pass 3: **VERDICT: CORRECT**
  ("the pass-2 repairs close the remaining gap; the proof is complete within the stipulated scope"). Confirmed throughout: expansion, (∗), constants and signs, the three S_h cases, Q0 ≤ V^{2/3},
  V²/Q0 = Λ^{2k}V^{2l}, uniformity, and that Lemma C′ gives Lemma C(c)'s bound.
- Reviewer B (range table). Pass 1 GAP: hypotheses of the non-asymptotic form; integer Q; |g_h′| and ‖g_h′‖ in Kusmin–Landau; Heath-Brown's λ5 > 0; Bourgain/Huxley hypotheses not supplied;
  strictness ν < ν1 and the pair-dependent B3N rows; the Trudgian–Yang row missing. All repaired; rows 4–7 now carry TY's hypothesis verification (found after the pass). Pass 2 GAP on row 5 only
  ("checks the Bourgain input, not Sargos's D-process application"; "an audit gap, not a demonstrated false estimate"; every other row accepted, inherited items as marked): row 5 now states
  TTY's Lemma D-process and TY's use of Sargos Thm 7.1 verbatim, checks D(13/84, 55/84) = (18/199, 593/796) and the max on the load-bearing range (logs/b3p_dcheck.log). Pass 3: **VERDICT: CORRECT** ("no remaining gap identified within this audit's stated inheritance scope").

## 7. Novelty (adversarial) and class
arXiv listings math.NT, math.CA 2026-10, 2026-11 (.scratch/lst.py, regex exponent pair|exponential sum|van der Corput|Bombieri–Iwaniec|decoupling|short interval|simple zero|Levinson|mollif;
logs/novelty.log, 24 titles) and Zenodo ("exponent pair" zeta): no new exponent pair, no short-window mollified moment; ANTEDB (teorth/expdb main, beta.tex, exponent_pairs.tex) lists the same
best row (89 + 908α)/1282 at α = 1/3 and the same four TTY pairs. B3N's same-day scan covers the short-window theorem. Known ingredients: A-process in the model class (TTY/ANTEDB, Ivić),
TTY's pair. New here: the model-phase verification of Lemma C's dual phases (twist removed, uniform in g, p, q, u, N, σ, shifts), the provenance audit of TTY's pair (Bourgain/TY rows load-bearing,
Bourgain's log-only caveat closed by TY's GK-class verification), the certified threshold. Class: RECORD (range/constant improvement of B3Nu's node; RECORDS.md), not a discovery.

## 8. What remains (stakes; logs/b3p_thr_stake.log)
Lemma C's gain is now capped by β(1/3): s = 3·β(1/3) ≥ ½. Under the exponent pair conjecture: A(0, ½) = (0, ¾) gives θ1 = 0.490046; (0, ½) applied to F_u directly (F_u is model, §1.3(c))
gives θ1 = 0.485302, the ceiling of every 1D-pair version of Lemma C with this detector (B3N's floor ≈ 0.479 is the detector's). Of the 0.00867 between B3N's 0.49397 and that
ceiling, the best known pair gains 0.00007.
Further progress needs B3N §8's (U) (2D/Huxley-type bounds for the divisor-split sum u·v uniform in the twist) or (G) (Kloosterman bilinear bounds on the (g, p, q)-sum), or a better β(1/3).
Suggested RECORDS line (N8): simple critical zeros in every (T, T + T^θ], θ > 0.49390 (b3pair: Lemma C′ = Lemma C with TTY's model-phase pair A(89/1282, 997/1282) = (89/2742, 1184/1371),
ν < (7478θ − 2457)/10042; reviewer CORRECT; Arb 4000 = 3000 bits); was 0.49398 (b3nu). Fallback without Bourgain's Theorem 6: A²(Huxley 2005), θ > 0.49392.
TREE N8: exponent-pair route of the short-window dual saturates at β(1/3) (EPC ceiling 0.4853 with this detector); next (U)/(G). NOGO: none.

## Files
REPORT.md (this), LEMMA.md (Lemma C′, range table). code/: b3p_ttycheck.py (exact re-run of TTY's table check, Bourgain/TY dependence), b3p_gain.py (gain LP), b3p_thr.py (copy of
b3nu_thr.py, pair list replaced), b3p_arb.py (copy of b3nu_arb.py), b3c_arb_copy.py (copy of B3Nu/code/b3c_arb_copy.py). logs/: b3p_ttycheck.log, b3p_gain.log, b3p_thr.log, b3p_thr_hux.log,
b3p_thr_stake.log, b3p_arb_tty_{4000,3000}.log, b3p_arb_ttyd_{4000,3000}.log, b3p_arb_huxd_{4000,3000}.log, b3p_arb_footprint.log (max RSS 35 MB), b3p_thr_r515.log, b3p_dcheck.log, novelty.log.

RH relevance: N8
