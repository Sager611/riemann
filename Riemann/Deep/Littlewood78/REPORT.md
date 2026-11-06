RECORD: |L(1, χ)| ≥ 0.11/max(3, log log q) for every primitive χ mod q ≥ 3 (0.115/log log q once log log q ≥ 3.13, 0.116 once ≥ 4, tending to (π²/6)e^{−γ}(1 − B) = 0.11694); h(−d) ≥ 0.0350·√d/max(3, log log d) for every fundamental discriminant −d < −4 (0.0366·√d/log log d once log log d ≥ 3.13).

# Explicit Littlewood bound for L(1, χ) from the zero-free half-plane (node N11)

Revision 1 (2026-11-05, 20:58 CET): the reviewer found a gap in the tail lemma (Lemma T): it did not assume κ_1 < 1, and it
claimed that τ(x) decreases, which is false (τ jumps up at every prime). Fixed in PROOF.md: Lemma T now assumes
κ_1(Λ_N) < 1 and proves κ_1 < 1 on the whole tail (κ_1 decreases in Λ); Lemma 9(d) defines the decreasing majorant τ̄ ≥ τ,
which Prop. 10 and Lemma T use (cert.py always evaluated τ̄); the bracket that bounds W is certified negative for L ≥ 8; the
upper-bound tail got the same hypothesis. cert.py now prints κ_1(Λ_N) for every tail (≤ 9.6·10⁻¹⁸ at the worst) and the
bracket checks (section (0)). Both certificates (B = 0.87338196661320752 and B = 7/8) were rerun: all claims PASS and
every number is unchanged. A re-review of the revised Lemma T returned CORRECT.

Agent littlewood78, 2026-11-05, 19:08–20:40 CET (revision 1: 20:49–21:00 CET). Files in this folder:

- `PROOF.md`: statements and complete proof (Lemmas 1–9, Propositions 10–11, tail Lemma T).
- `cert.py`: the certificate (python-flint arb, 256-bit balls); `cert_output.txt`: its output (all claims PASS).
  Run: `uv run --no-project --with python-flint python3 cert.py cert_output.txt` (single thread, about 3 min, < 100 MB).
- `cert_output_7_8.txt`: the same certificate with B = 7/8 (OpenAI's theorem, Lean port in Riemann/Deep/QRH) as fallback.
- `xcheck.py` / `xcheck_output.txt`: independent mpmath re-implementation; `sanity.py` / `sanity_output.txt`: numerical
  check of the explicit formula on small conductors; `dusart_check.py` / `dusart_check_output.txt`: the two cited Dusart
  bounds re-checked numerically on y ≤ 10⁸.

## What this is

L(1, χ) is the value at s = 1 of a Dirichlet L-function; for χ = χ_{−d} it equals π h(−d)/√d, so lower bounds for it are
lower bounds for class numbers. Littlewood (1928) showed under GRH that L(1, χ) ≫ 1/log log q. The 2026 zero-free
half-plane (no zeros of any Dirichlet L-function with Re s > B = 0.87338196661320752, Lean-checked) is enough for the same
shape, as noted in OpenAI's quasi-RH paper and in arXiv 2610.21511 (remark at TeX l.486–499: "We have not made it explicit
under the half-plane"). This report makes it explicit, with every constant certified.

## Results (unconditional given the half-plane)

| statement | range | source |
|---|---|---|
| \|L(1, χ)\| ≥ 0.11/max(3, log log q) | every primitive χ, q ≥ 3 (so 0.11/log log q for q ≥ e^{e³} ≈ 5.3·10⁸) | Thm 1(1) |
| \|L(1, χ)\| ≥ 0.115/log log q | log log q ≥ 3.13 (q ≥ 8.6·10⁹) | Thm 1(2) |
| \|L(1, χ)\| ≥ 0.116/log log q | log log q ≥ 4 (q ≥ 5.2·10²³) | Thm 1(2) |
| \|L(1, χ)\| ≥ 0.115659 / 0.116424 / 0.116873 / 0.116933 over log log q | log log q ≥ 40 / 100 / 10³ / 10⁴ | Thm 1(3) |
| liminf \|L(1, χ)\| log log q ≥ (π²/6)e^{−γ}(1 − B) = 0.1169398 | q → ∞ | Thm 1(3) |
| h(−d) ≥ 0.0350·√d/max(3, log log d) | fundamental −d < −4 | Cor 2 |
| h(−d) > 100 / 10³ / 2000 / 10⁴ / 10⁵ | d ≥ 6.973·10⁷ / 7.292·10⁹ / 2.957·10¹⁰ / 7.627·10¹¹ / 7.976·10¹³ | Cor 2 |
| side result: \|L(1, χ)\| ≤ 15 log log q (constant → e^γ/(1 − B) = 14.07) | log log q ≥ 3 | Thm 3 |

Fallback with B = 7/8 (OpenAI's theorem, already ported to Lean here; `cert_output_7_8.txt`): |L(1, χ)| ≥ 0.105/max(3, log log q),
≥ 0.11/log log q for log log q ≥ 3, ≥ 0.114/log log q for log log q ≥ 4 (limit 0.11545); h(−d) > 2000 for d ≥ 3.047·10¹⁰,
h(−d) > 10⁵ for d ≥ 8.137·10¹³; |L(1, χ)| ≤ 15 log log q for log log q ≥ 3 (limit 14.25).

Pointwise certified values (optimised parameters): |L(1, χ)| ≥ 0.0424 at log log q = 1, 0.0373 at 3, 0.0327 at 4,
0.0278 at 5, 0.0135 at 10, 0.0063 at 20, 0.0030 at 40; the product with log log q peaks at 0.1408 near log log q = 6.
The bounds hold for complex χ as well (the worst case χ(p) = −1 is available to real characters, so nothing is lost).

## Before

- Explicit unconditional lower bounds were of size 1/log q only: L(1, χ) ≥ 0.0005744/log q for primitive real χ and
  h(−d) ≥ 0.0001828446364·√d/log d (2610.21511 Prop. L1 with the explicit Siegel-zero gap; R305, R320); for non-real χ,
  |L(1, χ)| ≥ 1/(9.69030 log(q/π)) (Mossinghoff–Starichkova–Trudgian, arXiv 2107.09230, classical inputs only; ours is
  larger once q ≥ 53); our own Lean bound e^{−101}/log q (QRHCor).
- Explicit 1/log log q bounds existed only under GRH: Lamzouri–Li–Soundararajan (Math. Comp. 2015, arXiv 1309.3595,
  Thm 1.5, Cor 1.3), Languasco–Trudgian (2011.08348), Lumley (1804.09850), Chirre–Molero (2602.06199).
- With the half-plane, 2610.21511 (remark) obtained only thresholds: h(−d) > 2000 for log d ≥ 27.2 and h(−d) > 10⁵ for
  log d ≥ 35.4.

Crossover with 0.0001828446364·√d/log d: there is none. Our bound is larger for every d > 4, by a factor > 120 (124 at
d = 7, 1277 at log d = 20, 1885 at log d = 35, 4158 at log d = 100). Class-number thresholds improve from log d ≥ 27.2 to 24.11
(h > 2000: d ≥ 2.957·10¹⁰ instead of 6.5·10¹¹) and from 35.4 to 32.01 (h > 10⁵: d ≥ 7.976·10¹³ instead of 2.4·10¹⁵).
So with the half-plane every imaginary quadratic field with h ≤ 2000 has d < 2.96·10¹⁰, deep inside the Mosunov–Jacobson
tables (d < 2⁴⁰; 2610.21511 had already noted that its 6.5·10¹¹ lies inside them), and for h ≤ 10⁵ the search range
shrinks 30-fold, from 2.4·10¹⁵ to d < 7.98·10¹³ (still beyond the tables; the sieve was not attempted).

## Method (details in PROOF.md)

1. Smoothed explicit formula (an exact identity): −L'/L(σ, χ) = Σ Λ(n)χ(n) v(log n/L) n^{−σ} + Σ_ρ V(ρ − σ) with the
   weight v = 1 on [0, 1], linear down to 0 at 1 + h, scale x = e^L, and V(w) = (e^{(1+h)Lw} − e^{Lw})/(hLw²).
   Integrating over σ ∈ [1, ∞) gives log L(1, χ) as a prime sum plus Σ_ρ G(ρ).
2. The zero sum needs no zero counting: |G(ρ)| ≤ g(β)/|1 − ρ|², and the exact Hadamard identity
   Z = Σ_ρ (1 − β)/|1 − ρ|² = Re L'/L(1, χ) + ½ log(q/π) + ½ ψ((1 + a)/2). Pairing each zero with its functional-equation
   partner 1 − ρ̄ gives Σ|G| ≤ (g(B)B/b + g(1/2)) Z; the worst height is γ = 0 (low zeros).
3. Re L'/L(1, χ) is bounded by the same explicit formula at s = 1 with a second, shorter weight (L_1, h_1), which closes
   the loop: Z ≤ (−Re P_1 + ½ log(q/π) − γ/2 + T_1)/(1 − κ_1).
4. Prime sums: worst case χ(p) = −1 for every p up to x^{1+h}, evaluated with Dusart's explicit Mertens theorems
   (arXiv 1002.0442, Thms 6.11 and 6.12); this gives the factor (π²/6)e^{−γ}/log x.
5. Choosing x ≈ (log q)^{1/(1−B)} makes the zero term small, so log x ≈ log log q/(1 − B) and the constant tends to
   (π²/6)e^{−γ}(1 − B). With B = 1/2 the same argument returns Littlewood's GRH constant π²/(12e^γ), a check of the shape.

## Certificate

- `cert.py` evaluates the master inequality (Prop. 10) in arb balls at exact rational parameters (chosen by a float
  optimiser, which plays no role in correctness). F decreases in log q, so each cell [Λ_i, Λ_{i+1}] of Λ = log log q is
  checked at its right end against c/Λ_i. Grids: 7400 cells on [3, 40] for 0.11 (worst margin 0.0174 in log), 7374 cells
  on [3.13, 40] for 0.115 (worst margin 0.00047), 9600 cells on [4, 100] for 0.116 (worst margin 0.0228), 7400 cells on
  [3, 40] for the upper bound 15 (worst margin 0.0188); tails by Lemma T: one evaluation at Λ_N covers [Λ_N, ∞), because
  κ_1 and every error term decrease in Λ (with τ̄ in place of the non-monotone τ). Side conditions certified in section (0)
  (brackets ≤ −1.7237 and ≤ −1.4097 at L = 8, h = 2) and per tail (κ_1(Λ_N) ≤ 9.6·10⁻¹⁸, 8.4·10⁻⁴⁴, 6.1·10⁻⁴³⁵,
  1.4·10⁻⁴³⁴³ for Λ_N = 40, 100, 10³, 10⁴).
- An independent re-implementation of Prop. 10 in mpmath, written from the PROOF.md text (`xcheck.py`, output
  `xcheck_output.txt`), agrees with cert.py at the checked cells, including the tightest cell of (2) (margin 0.00047119972
  in both).
- Numerical sanity check of the explicit formula (`sanity.py`, output `sanity_output.txt`): for D = 5, −4, −3, 8, −7, 12, −23, 13
  the identities of Lemmas 2–3 hold with differences inside the GRH-version bounds; Z for χ_{−4} (0.0778) matches the sum
  over its first zeros plus the tail estimate.
- The cited Dusart bounds hold on all primes up to 10⁸ (`dusart_check.py`): worst ratios 0.9412 (Thm 6.11, at p = 10463)
  and 0.9235 (Thm 6.12, at p = 3947), both ≤ 1; we only use them for y ≥ e⁸ and y ≥ e¹⁰.

## Self-review (two passes)

1. Branch of log L near s = 1: only log|L(1, χ)| = Re log L(1, χ) is used, with log L(1, χ) = −∫_1^∞ L'/L(σ, χ) dσ; the
   half-plane gives L(σ, χ) ≠ 0 on [1, ∞). No branch choice enters.
2. Real vs complex χ: the per-prime bounds use only |χ(p)| ≤ 1; for complex χ the zeros are still symmetric under
   ρ ↦ 1 − ρ̄ (the zeros of L(s, χ̄) are the conjugates of those of L(s, χ)), which is all Lemma 6 needs.
3. Smoothing error: the explicit formula is exact; trivial zeros have the favourable sign or are bounded by e^{−L}
   (Lemma 7); the pairing lemma was re-derived (linear in r = (α² + u)/(c² + u), checked at both endpoints).
4. Uniformity in q: Prop. 10 holds for every q ≥ 3; F decreases in log q, which makes the cell checks and the small-q
   case (log log q ≤ 3, one evaluation at q = e^{e³}) valid.
5. Region |t| small: no separate treatment is needed; the identity for Z weighs every height, and in Lemma 6 the worst
   height for every β is γ = 0. Pass 1 found the factor-B gain from pairing; pass 2 found that Lemma 6 needs no extra
   condition at r = 1 (K(β) ≥ g(β) + g(1 − β) holds automatically) and re-checked every line of cert.py against Prop. 10.
6. A `reviewer` agent (focused on Lemma 6 and the assembly (Z-bound) → (★) → Prop. 10) returned: pairing inequality and
   assembly correct; it flagged three imprecise remarks (the gain is B + b·g(1/2)/g(B), not exactly B; the real pair
   attains g(B)B/b + g(b)b/B < K; C_0 = −0.0795153…, printed as −0.0796) — all fixed in PROOF.md; none affects a bound.
7. Second reviewer (orchestrator): GAP in Lemma T (missing κ_1 < 1 on the tail; false monotonicity of τ). Fixed as in
   Revision 1 above; the same reviewer found Lemma 6 for complex χ, the Z identity, the loop closing and Lemma 8 correct.
   While fixing, the limit statement was re-proved with explicit choices (h = 1/Λ_N gives κ_1(Λ_N) ≤ 2Be^{−Λ_N} +
   2be^{−Λ_N/(2b)}, κ_0e^{Λ_N} ≤ 2b/Λ_N, E_h ≤ h/2), and Corollary 2's monotonicity step now states its range (ℓ ≥ 20).
   A focused `reviewer` re-check of the revised Lemma 9(d), Prop. 10, Lemma T, the limit paragraph and the upper-form tail
   returned VERDICT: CORRECT.

## Caveats

- Depends on the half-plane Re s > 0.87338196661320752 for every Dirichlet L-function (Bebel's Lean check, trusted as in
  TREE.md N7); the 7/8 fallback above loses only about 1.3% (factor (1 − 7/8)/(1 − B) = 0.987 asymptotically).
- Uses Dusart's Theorems 6.11 and 6.12 and standard facts from Davenport (Ch. 12, 16, 17; the latter only qualitatively).
- The method's ceiling is (π²/6)e^{−γ}(1 − B) = 0.11694; going beyond needs more than a zero-free half-plane (for example
  a zero-density input).
- Two review passes (remarks and the Lemma T gap fixed; no bound changed); not formalized in Lean.

## Novelty pass (2026-11-05)

- arXiv math.NT monthly listings 2026-09 (812 entries), 2026-10 (854), 2026-11 (140): title scan for L(1, class number,
  Littlewood, quasi, half-plane, 7/8, log log, Siegel, zero-free; abstracts read for 2610.18491, 2610.21511, 2610.31402,
  2610.15483, 2610.30476, 2611.03356, 2610.15291, 2610.12234, 2610.21099, 2609.35950, 2609.30207, 2610.16953. The TeX of
  2610.21511 and 2610.18491 was searched: the former states that the explicit half-plane version is open, the latter has
  no L(1, χ) bound.
- All 102 works citing LLS 2015 (Semantic Scholar): every explicit 1/log log q bound is under GRH or RH to a height
  (2310.06723, ζ only); the unconditional ones are of size 1/log q.
- arXiv API full-archive searches (quasi-Riemann hypothesis; zero-free half-plane + class number; Littlewood + L(1;
  half-plane + log log q; quasi-GRH; 7/8 + class number): nothing relevant beyond the papers above.
- Zenodo API (7 queries): nothing relevant then; harvest2 (6 Nov 2026) later found B. W. Daugherty, Zenodo 23400116 (15 Oct) and 23459046 (18 Oct, Lean 4): L(1, χ) ≥ 0.0617 and 0.0621/log log q for real primitive χ with q ≥ 10^100, from Re s > 7/8 (weaker than ours: real χ only, q ≥ 10^100, smaller constant). Our repo: only the e^{−101}/log q Lean bound (QRHCor) and the 1/log q
  constants of R305/R320.
- Verdict: new as far as we can see. The shape is the expected consequence of the half-plane; the new part is the
  explicit certified constants, the class-number thresholds and the side upper bound.

RH relevance: N11 (GRH families, Siegel zeros, class numbers): the first explicit unconditional Littlewood-type bound. It is
the strongest explicit lower bound for L(1, χ) valid uniformly in q: for real χ at every conductor against the earlier
uniform bounds (below d = 7342 the trivial h(−d) ≥ 1 is stronger), for complex χ once q ≥ 53.
