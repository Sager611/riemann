MODEL: OpenAI's Lean-checked Lemma 18.1, z = 0 clause (ZeroAt, Energy/Bands.lean:22, delivered by terminal_certificate, Energy/CertifiedExistence.lean:42), contains Remark 12.18 for F = Q(√−3) at every bounded length, but only under the Lean hypothesis 51/100 ≤ β_* (false if GRH holds for the family; a packaging artefact); for F ⊃ μ12 the clause stays unrefereed in its prose constructions (G1), with no false step: the order-dependent step is re-derived and exactly tight; numerics to U = 10^8 (Q(√−3)) and 3·10^7 (Q(ζ12)) show no signal.

# Liu 2610.26835v2: Remark 12.18 and the z = 0 clause of M2 (liusext, 11 Nov 2026)

Brief: .scratch/orca/liusext-brief.md. "l." = line of the e-print number_fields_v4.tex (13,294 lines); "o.l." = line of OpenAI's paper.tex (openai/math 1581d4cb, 16,677 lines); Lean paths are under Riemann/Deep/QRH/OAI/NumberTheory/DirichletL/. Derivation: DERIVATION.md.

## Short answer
- **(i) Lean model: yes for Q(√−3), with one condition.** The Lean statement ZeroAt bounds Σ_{z∈O, kept} |(X1X2)^{−1/2} S_z(X1) S_z(X2)|² Φ(N z/Z^{rowWidth}) by Z^{width+ε}, for every X_i ≤ Z^L and any L. Here S_z(X) = Σ_I η(I) (m⁶z/I)_6 N(I)^{it} W(N I/X). The balanced case (rowWidth = 1, X_i = Z^{1/2}) is included, and so is every other bounded length. The condition: terminal_certificate assumes 51/100 ≤ β_*, the supremum of the real parts of zeros of the Hecke family of Q(√−3).
- **(ii) Referee, general F: no false or unjustified step on the z = 0 path.**
  - The whole Lindelöf strength rests on one balance: the exceptional rows have the sixth-power count Z^{(m′−f)/6}, and centering with L = M/4 saves 2L/3 = M/6. The two match exactly; I re-derived both.
  - Without centering, the same recursion gives U^{1+1/n}. That is the optimal large sieve for arbitrary coefficients (de Faveri).
  - Still open is G1: the prose constructions of Lemma 12.13 and App. C.4.
- **(iii) Numerics: no falsification signal.**
  - Over Q(√−3), M4 = (0.978 ± 0.005) × rows × E4(X) for every U from 10^6 to 10^8. Here E4 is the value under a random multiplicative model. The ratio drifts like U^{−0.0001}.
  - Over Q(ζ12), the first field with infinite units, the ratio is 1.00 ± 0.01 up to U = 3·10^7.
  - The fitted power U^{1.012} is just the log X growth of E4.
  - Sextic and cubic behave alike, as expected: these numerics cannot tell a proved bound from an open one.

## (i) Lemma 18.1 versus the balanced sextic fourth moment
Read in the Lean statements (not in their doc comments):

| object | Lean | content |
|---|---|---|
| ZeroAt Q a b bΦ Bmask L M ε Z degree S C | Energy/Bands.lean:22–27 | for every state s with fixedModulus = Q and width ≤ M, every profile pair p (supp ⊆ [a,b]), every t and 0 < X1, X2 ≤ Z^L: plainEnergy ≤ C · diagonalControl(Φ) · control² · (1+‖t‖)^degree · Z^{width+ε} |
| NaturalState | Energy/State.lean:17–36 | character η (modulus norm ≤ Z^{characterWidth}); puncture (norm ≤ Z^{Bmask}); Radial (keep, Φ, scale = Z^{rowWidth}); kept rows z ≠ 0 and ¬FixedInducingRow; width = rowWidth + characterWidth (paper: M = m + q) |
| plainEnergy | Energy/State.lean:46–50 | energy with an empty slot list |
| energy | Moments/InductionEnergy.lean:18–22 | Σ'_{z∈O} [keep z] ‖positiveSlotRow‖² Re Φ(‖z‖²/K) |
| positiveSlotRow | Moments/RetainedEnergy.lean:16–21 | (X1X2)^{−1/2} · rowTwistedSum(X1) · rowTwistedSum(X2) (the empty slot product is 1) |
| rowTwistedSum | Moments/HeckeHeight.lean:127–129 | Σ_I η(I) · idealRowHom(m⁶A₀z)(I) · N(I)^{it} · W(N I/X) |
| idealRowHom | Eisenstein/CompactEnergyFamilies.lean:468–486 | sextic residue symbol (actualSextic) at primes P ∌ λ; 0 at λ |
| FixedInducingRow | CenteredExceptionalProfile.lean:340–343 | the row character is induced by a primitive ψ whose conductor divides Q (with 72 ∣ Q) |
| terminal_certificate | Energy/CertifiedExistence.lean:42–52 | hypotheses: hκ0: 3/4 ≤ κ; hbeta: 51/100 ≤ β_*; hκ: 2β_* − 1 ≤ κ; 0 < a ≤ 1/4 ≤ 1 ≤ b; W smooth |
| certified_terminal | Energy/CappedWidthInduction.lean:96 | gives ZeroAt at width Mcap and length L, for all large Z |
| β_* | Hecke/ZeroSupremum.lean:12–18 | sSup({1/2} ∪ {Re ρ : L(χ, ρ) = 0, Re ρ > 0, ρ not the pole}) |

- **Specialisation.** Take:
  - η trivial (characterWidth → 0) and puncture 1;
  - rowWidth = 1 and Φ ≥ 1_{[0,1]};
  - X1 = X2 = Z^{1/2} and t = 0;
  - keep = nonexceptional.

  This gives Σ_{0<N z≤Z} |Z^{−1/4} S_z(Z^{1/2})|⁴ ≪ Z^{1+ε}.
- **Liu's family is covered.** Liu takes one allowed generator per ideal prime to S. His family is a finite union of sub-collections of OpenAI's family, each twisted by a fixed S-unit character; that twist is a fixed η. All terms are nonnegative. Rows that induce a fixed nonprincipal character are added back by Liu's Lemma 12.14 (l. 7516–7538); for z = 0 this is trivial.
- **So Lemma 18.1 implies Remark 12.18 over Q(√−3), as Lean statements.** This holds at every bounded length, including Liu's range m ≤ 1/2 + r_F/100. No regime is missing.
- **The one condition.**
  - κ = 1 satisfies hκ0 and hκ, the latter by beta_le_one.
  - hbeta is not provable. If GRH holds for the family, then β_* = 1/2 and the certificate says nothing.
- **Why the hypothesis is there.** Lean's zero stage actual_zero_stage (Energy/NaturalZeroStage.lean:26–54) consumes PositiveAt, the prime-slot clause, at the child width. Through actual_growth_physical (Energy/ZeroGrowthPhysical.lean:33) it routes the plain case through the general slot machinery with an empty slot set. The paper does it the other way round: all zero-slot bands come first (o.l. 12592), and Liu's rank does the same (l. 7622–7640, 7661–7663).
- **Static scan.** 72 declarations carry the hypothesis (any spelling). Exactly one uses it otherwise than by forwarding it: Moments/PrimeGlobal.lean:14–24, zeroMaximum_le_beta.
  - That lemma bounds the zero maximum of a finite character family by β_*. The zero maximum carries the sentinel 51/100, so the sentinel itself must lie below β_*.
  - Its only consumer is ray_prime_global_squared (PrimeGlobal.lean:30), the bound for prime sums. There 51/100 is also the lower end of the contour lines a ∈ [51/100, 1] (e.g. Hecke/PrimeDyadicControl.lean:71–95).
- **Verdict: almost certainly a packaging artefact.** At κ = 1, running the chain with b = max(β_*, 51/100) in place of β_* would discharge both hypotheses. That is a Lean refactor (make the chain parametric in b ≥ β_*, or specialise the slot type α := Empty) and needs a build; I did not do it.
- **Routine side bound (COROLLARIES type; reviewer: sound).**
  - If β_* < 51/100, every row L-function is zero-free on Re s > 0.51. The row characters are pole-free and their conductors are polynomial in U.
  - The classical lemma "zero-free half-plane ⇒ Lindelöf there", used uniformly in the conductor, then gives |S| ≪ X^{0.51+δ}U^ε, hence Σ ≪ U^{1.02+ε}.
  - Splitting into the two cases, this gives Σ ≪ U^{51/50+ε} unconditionally, but ineffectively. The argument is analytic; the Lean declarations certify only the case β_* ≥ 51/100.
  - The bound beats the large-sieve route (U^{7/6+ε}) but is weaker than the paper's claim. It is not a discovery.
- **Reviewer** (verdict-only, with the Lean statements pasted): "SOUND, with scope qualifications". The qualifications: fixed η of bounded conductor; removability of hbeta plausible, not established; the side bound needs conductor uniformity. All three are stated here.

## (ii) Referee of the z = 0 clause (Thm 12.16, Remarks 12.15 and 12.18, App. C.4)
- **The path.**
  - The rank (l. 7622–7640) treats plain inputs first.
  - z = 0 uses no amplifier (w = w_o = ℓ = 0, l. 6774–6776), no prime sum and no zero-free hypothesis.
  - Reflection: Lemma 12.6 (l. 6149). Centering: P − P_comp with Y1Y2 = X1X2 and L = M/4 (l. 7650–7660).
  - Two transforms (l. 6576–6981), then the exceptional rows (l. 7146–7327), then the add-back (l. 7516).
  - Liu's §12 transcribes OpenAI §18: the threshold 5M/6 (o.l. 12594), L = M/4 (o.l. 12944), v_p(h′) ≡ 4 mod 6 (o.l. 14356) and the maximum A − M at v = L (o.l. 14765–14767) are all OpenAI's.
- **The double transform is an identity up to one positive enlargement.**
  - T1 sends the rows to Gauss-sum frequencies h (first bridge, l. 6628). The h-range is then enlarged from K to K + g, with g = J₊ + σ. T2 sends h to j (second bridge, l. 6850).
  - The Gauss sums cancel in the correlation F(D2a, E2b; j) = F(D2, E2; j) χ_a(j) χ̄_b(−j) (Lemma 12.10). The child is therefore the same plain family at width M′ = M + J − g − g2 + t2 ≤ M − σ (l. 6922–6927).
  - Checked: the identities (LiuAudit) and the width and allowance accounting (LiuAudit, plus a one-line model check, DERIVATION §2).
- **Exceptional rows, the only step that depends on the order of the characters.**
  - Kummer theory forces v_p(h′) ≥ 4 (l. 7157–7174). The count is Z^{(m′−f)/6} (l. 7189) and the raw volume is Z^{a0−b2} (l. 7199–7200).
  - exc_orders.py derives, symbolically and for any order n: E_x = count + volume − (M′ + Δ′) = A − (1 − 1/n)M − F1^{(n)} − F2^{(n)}. At n = 6 this is exactly the printed F1, F2 (l. 7203–7214).
  - F2 ≥ 2b2/3, with equality at i = 1 nonunit. F1 ≥ (2/3)(c + w) − 3σ, also tight.
- **Centering** (l. 7257–7327). sup_v{−2v/3 − (L − v)₊} = −2L/3 = −M/6 at L = M/4, which cancels A − 5M/6 exactly at A = M.
  - The feasible window is L ∈ [M/4, 5M/12).
  - The paper's choice L = M/4 has slack 0; the best choice has slack M/9.
- **Consistency with the literature.**
  - Uncentered, E_x ≤ M/n, so Σ ≪ U^{1+1/n}. This is the balanced point of de Faveri's optimal large sieve for n-th order characters (arXiv 2610.04045, Thm 1.1: (AB)^ε(A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n})). Xu–Sang (2610.17675) prove it sharp for sextic characters.
  - Hence Remark 12.18 cannot follow from any large sieve with arbitrary coefficients. Centering of equal products of smooth plains is the one extra input.
- **Cubic.**
  - Needed saving: M/3.
  - Constraints: L < M/3; F2^{(3)} ≥ b2; and, if B_c is unchanged, F1^{(3)} ≥ (11/15)c.
  - The deficit is then at least 4M/45. Even with coefficient 1 the bound closes only at the boundary L = M/3. This agrees with Liu's "not decided". I did not re-derive B_c, the T1 brace or the Gauss tables for n = 3.
- **Steps specific to the field F** (table in DERIVATION §5): the lift to G/Γ, orbit Poisson, the finite decomposition into classes of U_S/U_S⁶ (l. 3537) and Θ ⊇ Υ_S invariance (l. 2027). Each is standard or already checked by LiuAudit.
- **First failing step: none.**
- **Open (G1, prose only):**
  - the positive enlargement of the h-rows to a smooth ball on G/Γ when the unit rank is ≥ 1;
  - the common Fourier separation of K1 and K2;
  - the measure bound δ_edge ≤ C_*(η + ξ + ε0) (Lemma 12.13, l. 7407–7514; App. C.4, l. 12741–13100).

  For Q(√−3) these are exactly the Lean-checked constructions.

## (iii) Numerical falsification
- **Rows and columns.**
  - Q(√−3): rows k = a + bω ≡ 1 mod 3 with N k odd and N k ≤ U, one generator per ideal prime to 6. Columns: all ideals prime to 6, weighted by W(x) = exp(1 − 1/(1 − log₂²x)) on (1/2, 2), with X = U^{1/2}.
  - The cubic symbol is the square of the sextic symbol. Principal rows (sixth powers, resp. cubes) are excluded.
  - E4(X) is the random model: iid sixth roots of unity on the prime ideals, zero with probability 1/Nπ, 4·10⁵ draws, Monte Carlo error about 0.5%.
- **Validation** (m4 test):
  - multiplicativity in k: 0 failures in 20,000;
  - cubic reciprocity for primary primes: 0 failures in 200,001.
- **Cost.** One process at a time, single-threaded; peak 12.6 MB (U = 10⁸, 1000 s CPU).

Q(√−3), X = U^{1/2} (logs/m4_ladder.log, m4_random.log, m4_prin.log, fit_m4.log):

| U | rows | M4 sextic | M4/rows | M4/(rows·E4) | M4 cubic | M4/(rows·E4) | excluded sixth powers / M4 | excluded cubes / M4 | s |
|---|---|---|---|---|---|---|---|---|---|
| 10⁴ | 3,027 | 224.79 | 0.07426 | 0.928 | 223.86 | 0.927 | 0.22 | 1.12 | 0.0 |
| 10⁵ | 30,235 | 2,475.0 | 0.08186 | 0.993 | 2,451.3 | 0.922 | 0.19 | 2.17 | 0.0 |
| 10⁶ | 302,299 | 25,240 | 0.08349 | 0.978 | 23,728 | 0.910 | 0.39 | 4.77 | 0.8 |
| 10⁷ | 3,022,972 | 257,583 | 0.08521 | 0.979 | 250,731 | 0.923 | 0.65 | 9.64 | 31.6 |
| 3·10⁷ | 9,069,008 | 784,671 | 0.08652 | 0.976 | 767,509 | 0.928 | 0.64 | 13.7 | 165.6 |
| 10⁸ | 30,229,965 | 2,664,636 | 0.08815 | 0.979 | 2,593,006 | 0.924 | 0.92 | 20.4 | 1000.3 |

- **Fits.**

  | quantity | sextic | cubic |
  |---|---|---|
  | M4 ∝ U^α, U ≥ 10⁶ | α = 1.0120 (RSS of log 4.7·10⁻⁵) | α = 1.0185 (RSS 2.4·10⁻⁴) |
  | M4/(rows·E4) ∝ U^β, U ≥ 10⁶ | β = −0.0001 | β = +0.0041 |
  | local exponents, top 3 | 1.012–1.015 | 1.011–1.018 |

- **Excluded principal rows** (the shadow of the paper's exceptional power count). Their contribution grows like U^{1.20} for sixth powers (7 rows at 10⁸; the heuristic is U^{1+1/6}) and like U^{1.337} for cubes (heuristic U^{1+1/3}).
  - At U = 10⁸ the sixth powers alone already equal 0.92 × M4, and the cubes equal 20 × M4.
  - This is the count behind A − (1 − 1/n)M at A = M. Excluding these rows (Liu: "ψ_k nonprincipal") is essential.
- **Other lengths, X = U^m** (logs/m4_lengths.log). M4/rows stays between 0.075 and 0.082:
  - m = 0.35, U up to 10⁷: 0.0750, 0.0796, 0.0818;
  - m = 0.65, U up to 3·10⁶: 0.0776, 0.0793, 0.0799.

**Q(ζ12)** (degree 4, unit rank 1; m4z12.c; logs/m4z12.log). This is the first field beyond the Lean model.
- **Rows.** One generator per ideal prime to 6 with N k ≤ U. The generator is taken in the unit fundamental domain |σ1k|/|σ5k| ∈ [ε^{−1/2}, ε^{1/2}), ε = 2 + √3, and fixed among the 12 roots of unity by the smallest residue code mod 3. This family is a finite union of unit-twisted sub-collections of Liu's family (Υ_S ⊂ Θ).
- **Columns.** All ideals prime to 6, with the same W and X = U^{1/2}.
- **Validation.**
  - multiplicativity: 0 failures in 20,000;
  - (ζ6/𝔭)_6 = ζ6^{(N𝔭−1)/6} at all 800 prime ideals;
  - norms integral to 10⁻⁹;
  - candidates = exactly 12 × rows.

| U | rows | M4 sextic | M4/rows | M4/(rows·E4) | M4 cubic | M4/(rows·E4) | s |
|---|---|---|---|---|---|---|---|
| 10⁵ | 24,048 | 1,368.5 | 0.05691 | 0.981 | 1,375.0 | 0.986 | 0.5 |
| 10⁶ | 240,791 | 14,200 | 0.05897 | 1.012 | 13,831 | 0.981 | 6.0 |
| 3·10⁶ | 722,290 | 42,490 | 0.05883 | 0.997 | 44,995 | 0.978 | 23.4 |
| 10⁷ | 2,406,991 | 139,962 | 0.05815 | 1.003 | 142,035 | 0.983 | 111.0 |
| 3·10⁷ | 7,221,387 | 431,105 | 0.05970 | 1.003 | 427,344 | 0.986 | 587.9 |

**Fits for Q(ζ12), U ≥ 10⁵.**

| quantity | sextic | cubic |
|---|---|---|
| M4 ∝ U^α | α = 1.0120 (RSS 4.0·10⁻³) | α = 1.0145 (RSS 9.5·10⁻³) |
| M4/(rows·E4) ∝ U^β | β = +0.0022 | β = −0.0004 |

Peak memory was 3.5 MB.

**Verdict.** No exponent above 1 beyond the random model's log growth. The off-diagonal is negative over Q(√−3) (−2.2% sextic, −7.6% cubic) and about zero over Q(ζ12) (+0.3%, −1.4%). In both fields it is constant in U. No exceptional excess appears among nonprincipal rows.

## Novelty: is a Lindelöf-on-average sextic fourth moment known?
- **No in the classical literature; yes in two 2026 preprints.**
  - Classical (scout sc1, arXiv): Heath-Brown 2000 and Blomer–Goldmakher–Louvel 2014 Thm 1.3 (1112.1650) give U^{4/3+ε}.
  - de Faveri 2610.04045 (2 Oct 2026) gives U^{1+1/n+ε}, i.e. 7/6 sextic and 4/3 cubic. That is optimal for arbitrary coefficients (sextic lower bound: Xu–Sang 2610.17675, which use OpenAI's unmarked moments).
  - U^{1+ε} for the sextic family: OpenAI's 7/8 manuscript (Sept 2026), Lemma 18.1 with z = 0, over Q(√−3). It is Lean-checked under 51/100 ≤ β_* (part (i)).
  - Liu's Remark 12.18 extends it to cyclotomic F ⊃ μ12 (unrefereed).
  - Cubic U^{1+ε}: open.
- **Scan of this and last month** (scout sc2: lst.py over 7 categories and Zenodo; logs/novelty_scan.log): no other claim. 2611.01452 (first moments under GRH) and 2610.15316 (non-vanishing over F_q[T]) were checked: neither has a fourth moment.
- **Our files:** R440, NOGO l. 447 and l. 499 (no phase-sensitive sextic or cubic estimate). Nothing contradicts.
- **Our own observations are consistency remarks, not discoveries:** uncentered M2 equals the optimal large sieve, and the cubic deficit is 4M/45.

RH relevance: N7 (and N11). Remark 12.18 is the strongest single input of Liu's number-field B_2 and B♭ (it feeds Prop. refined-numerator), and through them of the Artin, Chebotarev (B) and binary-form entries. Over Q(√−3) it has a Lean model, conditional on 51/100 ≤ β_*. Over F ⊃ μ12 it remains a human proof with no flaw found; keep the label "human, unrefereed".

Files (Riemann/Deep/LiuSext/):
- DERIVATION.md;
- m4.c (Q(√−3): modes run / rm / prin / test) and m4z12.c (Q(ζ12): run / rm / test). Compile with cc -O3, then run, e.g., ./m4 run 1e6 0.5.
- exc_orders.py (uv run --no-project --with sympy, 61 MB) and fit_m4.py (--with numpy, 39 MB);
- logs/{m4_ladder, m4_random, m4_prin, m4_lengths, m4z12, exc_orders, fit_m4, novelty_scan}.log.

## Follow-up (orchestrator questions, 11 Nov 2026, about 22:05 UTC)

**(1) Does OpenAI state the balanced moment explicitly? No.**
- **Paper.** paper.tex never singles out the balanced case and never says "Lindelöf" or "on average". (rg over all 16,677 lines for Lindel, on average, fourth moment, large sieve; the only hits are the ones quoted below.) The balanced moment is only an instance of Lemma 18.1:
  - o.l. 12552: "\ll Z^{M+\epsilon}";
  - o.l. 12557–12558: "\item \(z=0\).  There is no restriction on the bounded nonnegative lengths \(n_1,n_2\), and no lower bound on the conductor.";
  - o.l. 12481–12482: "The two plain polynomials have no length restriction when no prime polynomial is present".
- **Its only large-sieve statement** is the "Sextic large sieve" lemma (o.l. 4707–4723). It has strength (KD)^ε{K + D + (KD)^{2/3}} Σ|c_n|², which is U^{4/3} at the balanced point. The paper credits it to Blomer–Goldmakher–Louvel (BGL) Thm 1.3 and §3, and to Heath-Brown (o.l. 133–138). It never compares Lemma 18.1 with it.
- **Lean doc comments.** None: rg for fourth moment and large sieve in the port gives nothing. Lindelöf occurs only as Mathlib's IsLindelöf and Phragmén–Lindelöf.

**(2) Is it known in the literature? No, before OpenAI and Liu.**

| source | balanced normalised M4 | coefficients |
|---|---|---|
| Heath-Brown 2000, cubic over Q(ω) | U^{4/3+ε} (from the (MN)^{2/3} term) | arbitrary |
| Baier–Young 0804.2233, cubic and sextic Dirichlet characters over Q | large-sieve type, not U^{1+ε} | arbitrary |
| BGL 2014, Thm 1.3, n-th order over K ⊃ μ_n | U^{4/3+ε} | arbitrary |
| de Faveri 2610.04045, Thm 1.1 | U^{1+1/n+ε}: sextic 7/6, cubic 4/3 | arbitrary; optimal (sextic lower bound: Xu–Sang 2610.17675) |
| Dunn–Radziwiłł 2109.07463 | Patterson's bias of cubic Gauss sums, under GRH; not a fourth moment | — |
| OpenAI Lemma 18.1, z = 0 (Sept 2026) | U^{1+ε} over Q(√−3) | smooth plains; Lean-checked only under 51/100 ≤ β_* |
| Liu Remark 12.18 | U^{1+ε} over F ⊃ μ12 | smooth plains; unrefereed |

- The cubic U^{1+ε} appears nowhere.
- de Faveri's TeX contains neither "fourth moment" nor "Lindelöf".

**(3) The hbeta refactor.**
- **The Empty slot type alone does not remove hbeta.** The zero stages still call the prime-sum lemmas, with empty products, and their signatures demand hbeta. Those stages would have to be re-proved.
- **The cheaper route** replaces β_* everywhere by betaEff := max β_* (51/100). Then hbeta holds by le_max_right, and hκ: 2·betaEff − 1 ≤ 1 holds at κ = 1.
  - **Real proof change: only Moments/PrimeGlobal.lean.** zeroMaximum_le_beta becomes unconditional (the 51/100 sentinel case by le_max_right, the zero case by le_max_left). ray_prime_global_squared uses beta_le_one, which becomes betaEff_le_one, and its exponent (β − 1/2)r becomes (betaEff − 1/2)r.
  - **Elsewhere: a mechanical word-boundary substitution β → betaEff.** 64 files carry the hypothesis, with 160 occurrences of beta; port-wide 72 declarations carry it. No other β lemma (zero_re_le_beta, beta_le_one, half_le_beta, LFunction_ne_zero_of_beta_lt) is used in Energy/ or Moments/ outside PrimeGlobal (rg).
  - **New files:** one defining betaEff with three lemmas, and one unconditional wrapper theorem (terminal_certificate at κ = 1, specialised to the balanced ZeroAt).
- **Files that must be re-elaborated.** Everything downstream of PrimeGlobal inside the import closure of CertifiedExistence (2,231 files): 238 files with 32,151 lines (183 in Energy/, 55 in Moments/).
  - None imports Mathlib directly, so no Mathlib-heavy file is touched.
  - A fork under a new module prefix is clean, because no file of that closure outside the set imports one inside it. An in-place edit instead invalidates 250 port files.
- **Blocker: the QRH oleans are not in .lake.** There is no .lake/build/lib/lean/Riemann/Deep/QRH; .lake/build holds 5.2 MB, and only Mathlib's 6.0 GB of oleans are present. Any check therefore first builds the unmodified upstream closure, 2,231 − 238 ≈ 1,993 files.
  - **Time.** The port's own build (2,924 files) took 3.5 h of wall time at 4 Lean threads, about 17 thread-seconds per file. So the upstream part needs about 9.4 thread-hours: about 2.4 h at 4 threads, or 9–10 h single-threaded. The 238-file fork adds about 67 thread-minutes.
  - **Memory: an estimate, not measured here.** About 2–3 GB RSS per lean process, since Mathlib and the QRH oleans are memory-mapped. Four threads would therefore need 8–12 GB, next to two other Lean agents.
  - **Disk.** +1.4 GB of oleans, plus up to 8 GB of IR C files if Lake emits them as in the port build (these can be truncated); 56 GB is free.
- **Verdict: NO-GO for a 3-hour Orca agent now.** The prerequisite build alone takes at least 2.4 h at 4 threads, or at least 9 h single-threaded, on a shared .lake.
- **GO once the QRH closure oleans exist**, for example from a separate overnight single-threaded lake build of Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Energy.CertifiedExistence (9–10 h), or from a restored olean archive.
  - The job would then take about 1.5–2 h: about 45 min of scripted edits, about 70 min to build the 238 files single-threaded, plus fix-ups.
  - Mathematical risk is low, because only PrimeGlobal's proof changes. Mechanical risk is moderate: renaming 238 modules, and hidden uses of beta in tactic blocks would only surface as elaboration errors.
