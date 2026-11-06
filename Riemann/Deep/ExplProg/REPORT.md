RECORD: from the Lean-checked zero-free half-plane (B = 0.8733819666132076, 1/(1−B) = 7.8978), for every modulus q and (a, q) = 1: (i) |ψ(x;q,a) − x/φ(q)| ≤ 1.5 x^B log²(qx) for all x ≥ 2 (0.174 for x ≥ 64, 0.0324 for x ≥ e^10); (ii) the least prime p ≡ a (mod q) is ≤ (q log q)^{7.8978} for every q ≥ 2 (also ≤ (φ(q) log q)^{7.8978} for q ≥ 3), and ≤ (0.72 q log q)^{7.8978} for q ≥ 10^4; (iii) every nonprincipal character mod q ≥ 3 has χ(n) ∉ {0, 1} for some prime n ≤ (log q)^{7.8978}; (iv) |L'/L(1, χ)| ≤ 7.8978 log log q + 10.53 for every primitive χ mod q ≥ 3.

## Claim
Assume H: ζ(s) and every Dirichlet L(s, χ) are zero-free on Re s > B. H is a theorem: Bebel's Lean-kernel-checked submission to
QuasiRiemannTracker gives Re s > B₂ + 2.2·10⁻²⁵⁴ with B₂ = (777 − √1081)/852 = 0.8733819666132075199… < B (certified).
Then, with explicit constants and no other unproved input:

1. **Uniform PNT in progressions.** For all q ≥ 1, (a, q) = 1, real x ≥ 2:
   |ψ(x; q, a) − x/φ(q)| ≤ c₁ x^B log²(qx), with c₁ = 1.5 (forced up to 1.4854 by q = 1, x = 2); c₁ = 0.174 for x ≥ 64;
   0.0324 for x ≥ e¹⁰; 0.0147 for x ≥ e²⁰; 0.0112 for x ≥ e¹⁰⁰; 0.0079 for x ≥ e¹⁰⁰⁰; limit ½(b/4π + b²/π) = 0.0076, b = 1 − B.
2. **Linnik with numbers.** p(a, q) ≤ (c₂ q log q)^{1/(1−B)} with c₂ = 1 for every q ≥ 2, and for q ≥ Q:
   c₂ = 0.79 (Q = 10²), 0.75 (10³), 0.72 (10⁴), 0.67 (10⁶), 0.62 (10¹⁰), 0.56 (10²⁰), 0.53 (10³⁰), 0.51 (10⁵⁰), 0.49 (10¹⁰⁰),
   0.46 (10¹⁰⁰⁰); limit c_∞ = 0.45331 (every c₂ > c_∞ works for q ≥ q₀(c₂), effectively). The same holds with φ(q) log q in
   place of q log q for every q ≥ 3 (Theorem 2', same table).
3. **Least character non-residue.** For every nonprincipal χ mod q ≥ 3, the least n with χ(n) ∉ {0, 1} (a prime) is
   ≤ (c₃ log q)^{1/(1−B)} with c₃ = 1 for every q ≥ 3, and c₃ = 0.62 (log q ≥ 25), 0.60 (30), 0.56 (50), 0.52 (100), 0.47 (1000),
   0.46 (10⁴); limit 0.45331. For real χ this is the least quadratic-type non-residue; it beats the trivial bound q − 1 for q > e^{25.6}.

4. **Logarithmic derivative at 1 (item iv).** For every primitive χ mod q ≥ 3, |L'/L(1, χ)| ≤ (1/(1−B)) log log q + C with
   C = 10.53, and C = 0.352 (log q ≥ 25), −1.10 (≥ 100), −1.99 (≥ 1000), −2.14 (≥ 10⁴). Imprimitive χ: add Σ_{p|q} log p/(p − 1).
   This one is routine (harvest candidate 3, "COROLLARIES level").

## Before
- Explicit bounds valid for every q were exponential: Bennett–Martin–O'Bryant–Rechnitzer (Illinois J. Math. 62 (2018), arXiv 1802.00085,
  Cor 1.8 and Thm 1.1) prove |θ(x; q, a) − x/φ(q)| < x/(160 log x) only for x ≥ exp(8√q log³ q) (exp(0.03√q log³ q) for q > 10⁵),
  which gives p(a, q) ≤ max(exp(8√q log³ q), exp(φ(q)/160)) (our deduction; an exponential bound of this kind is quoted from BMOR in
  "Computing generalized convolutions faster than brute force", found by the GXL evidence search). No power saving, and no polynomial
  Linnik bound with numbers, existed unconditionally.
- Polynomial Linnik bounds have unspecified constants: Linnik; Xylouris L = 5; Pollack (Hexagon 2610.00545, "L absolute, effectively
  computable"); Naslund (Hexagon 2610.00018, q^{3.99} for q ≥ q₀, "we do not compute q₀ or C"); ours L ≤ 7/3 with ≪_ε
  (COROLLARIES.md, ArtinHalfPlane). GRH: Bach–Sorenson p ≤ 2(q log q)² type bounds; 2404.08380.
- Power-saving PNT in progressions from the 2026 half-planes: non-explicit only (Zenodo 23608576: x^{7/8} log²(qx), constant unspecified;
  OpenAI, 5 Oct: x^{11/12}). Platt–Trudgian (arXiv 2610.29018) said the explicit version "could be improved almost instantly".
- |L'/L(1, χ)|: unconditional C₅ log log q with C₅ ≈ 10¹¹ by formula only (Hexagon 2610.00256v3, Harvest1 evaluation); GRH gives
  2 log log q + O(1) (Ihara–Murty–Shimura).
- Least character non-residue: explicit unconditional bounds only of Burgess strength q^{1/4+o(1)} (e.g. explicit Burgess for composite
  q ≥ 1.627·10¹⁹, arXiv 2610.22939); polylogarithmic only under GRH (Bach: 3 log² q) or with an unspecified constant (OpenAI QRH,
  30 Sep, Cor 1.2: C(log p)³²).

## Method
DERIVATION.md has every step. In short:
- Explicit formula for weights supported in [t₀, t₁] ⊂ (0, ∞): Σ Λ(n)χ(n)w(n/x) = 1_{χ=1} xW(1) − Σ_ρ x^ρ W(ρ) − (trivial zeros),
  no pole at s = 0 and no bound for L'/L on a left line is needed.
- The involution ρ ↦ 1 − ρ̄ pairs every zero with β > 1/2 (weight x^B) with one with β < 1/2 (weight x^{1/2}): the zero sum is at most
  ½ Σ_ρ (x^B G(|γ|) + x^{1/2} G_L(|γ|)), where G majorises sup_{β∈[1/2,B]} x^{β−B}|W(β + iγ)| (Lemmas P, K1, K2).
- Σ_ρ G(|γ|) by partial summation against the refereed zero count of Bennett–Martin–O'Bryant–Rechnitzer (Math. Comp. 2021), with all
  zeros below height 5/7 put at height 0; closed-form antiderivatives give Z(L) = L z_a + 2n₀ log(1 + L) + z_b + z_e/(1 + L).
  For ζ: Σ 1/γ² ≤ 0.046371 from Σ Re 1/ρ = 1 + γ/2 − ½ log 4π (Theorems 2, 3) and Rosser's 1941 count (Theorem 1).
- Theorem 4: explicit formula at s = 1 for the ramp weight v(n/x) (pole of V at 0 gives −L'/L(1, χ)), the same pairing and BMOR sums,
  Dusart's Σ_{p≤x} log p/p for the prime part, x = (0.3 log q)^{1/(1−B)}.
- Theorems 2 and 3: triangle kernel in log n (length 7, peak at 2.1 from the top), positivity of the weighted prime sum
  (resp. of Σ Λ(n)(1 − Re χ(n))w(n/x)), and one monotone-majorant "tail" inequality per table row that covers all q ≥ Q at once.
- Theorem 1: sandwich ψ(√x; q, a) + S₋ ≤ ψ(x; q, a) ≤ ψ(√x; q, a) + S₊ with trapezoid weights of ramp width δ ≈ θφ(q)x^{−(1−B)};
  7307 certified cells in (log q, log x) for 64 ≤ x ≤ e⁶⁰, q ≤ x; an explicit tail for x ≥ e⁶⁰; x < q directly; q ≤ x < 66 exactly.

## Certificate
One command (14 s, about 50 MB, one thread; python-flint Arb, 200 bits; every inequality a certified ball sign, else abort):

    cd Riemann/Deep/ExplProg && uv run --no-project --with python-flint --with mpmath python3 cert.py all

Output: logs/cert_output.txt. The only non-Arb proof step is the exact least-prime sieve for 2 ≤ q < 100 in sanity.py (integer
sieve; there log p(a,q) ≤ 0.426·log((q log q)^{7.8978}) and, for q ≥ 3, ≤ 0.313·log((φ(q) log q)^{7.8978}); run:
uv run --no-project --with numpy python3 sanity.py, 7 min, ~0.5 GB; logs/sanity_output.txt).

## Checks
- Exact least primes for all 2 ≤ q ≤ 10⁴ and all classes (sieve to 10⁸): max log p(a,q)/log q = 1.83 (q = 5, p = 19); worst ratio
  log p/log((q log q)^{7.8978}) = 0.426 (q = 2), and log p/log((φ(q) log q)^{7.8978}) = 0.313 (q = 3, p = 7); no failure.
- Exact ψ(x; q, a), q ≤ 100, all a, 597 points x ∈ [2, 10⁸]: max |E|/(x^B log² qx) = 1.4848 at q = 1, x = 2 (certified sup 1.4854);
  for x ≥ 64 the maximum is 0.0093 (certified bound 0.174); max |E|/(√x log² qx) = 0.036 for x ≥ 100.
- Exact max_χ n_χ for 3 ≤ q ≤ 10⁴ (least N with ⟨n ≤ N coprime to q⟩ = (ℤ/q)^×): 43 (q = 2706).
- |L'/L(1, χ)| for all 4135 nonprincipal characters of prime moduli 3 ≤ q ≤ 200 (mpmath, digamma and Stieltjes γ₁):
  max 2.168 (q = 163, real χ), at most 9.3% of the Theorem 4 bound (lsanity.py, logs/lsanity_output.txt).
- Independent float re-computation of all zero sums (xcheck.py): closed forms ≥ direct Stieltjes sums, within 0.02–1.6% for
  log q ≥ 9.2. Random tests of the three kernel lemmas (ksanity.py, 3 × 20000 samples, max ratio 0.9991 ≤ 1) and of inequality (1.3)
  at x = 10⁷ for all q ≤ 30 (deviation ≤ 0.44% of the bound).
- Self-review against the brief's adversarial list: x < q (§6.1), q ≤ 4 (q = 1, 2 separate; q = 3, 4 conductors < 9 use log 9),
  prime powers (only ω(q) and Σ_{p|q} log p enter), real χ (covered by H, no exceptional zero), even χ (zero at s = 0 kept);
  certified-monotonicity tails; all roundings upward (DERIVATION §8, §9). One rigor slip found and fixed during review: a max() used as a
  lower bound for φ(q) now takes the certified larger ball only.

## Novelty
Searched 2026-11-06: arXiv math.NT listings Oct–Nov 2026 (lst.py, regexes 'progression|Linnik|least prime|non-?residue' and
'explicit|quasi-Riemann|zero-free|half-plane|Siegel|Dirichlet L|primes in'); Hexagon full-text search (Linnik, least prime,
arithmetic progressions, nonresidue, non-residue, quasi-Riemann, Dirichlet L-functions, Siegel zero, explicit constant) and the 10 newest
items; the TeX of Hexagon 2610.00178 (GPT survey of GRH consequences of a half-plane: no numeric Linnik, AP-PNT or small-non-residue
bound); Zenodo (7 queries); GXL evidence API (2 claims); arXiv abstracts of 1604.01750, 2108.10878, 1506.01635, 1802.00085, 2404.08380;
Tao's optimisation-constants page (Linnik's constant: 5, and 2.4*) and expdb; our repo (COROLLARIES l. 3 and l. 12, ArtinHalfPlane,
Harvest1). Nothing gives a numeric polynomial Linnik bound for all q, an explicit power-saving PNT in progressions uniform in q, or an
explicit polylogarithmic non-residue bound. Best previous numeric all-q results: exponential (BMOR 2018, above). Theorem 4 is new only
as a constant (routine). Details: NOVELTY.md.

## Caveat
- Everything is conditional only on H, which is Lean-kernel-checked (tracker admission) but not refereed; with OpenAI's 7/8 the same code
  gives exponent 8 (not run).
- Not refereed, no external reviewer yet. The constants are far from optimal: 72% of the leading coefficient comes from the ≈ 0.455 log q
  zeros below height 5/7 that the BMOR count cannot localise, all placed at height 0.
- Theorem 1's c₁ = 1.5 is set by x = 2; it says nothing new for small x, where trivial bounds are better.
- Imprimitive moduli are covered; the bounds use φ(q) ≤ q (Theorem 2) and the conductor bound q* ≤ q.

RH relevance: N10, N11

## Files
- REPORT.md (this file), DERIVATION.md (all steps and constants), NOVELTY.md (searches), cert.py (Arb certificate, sections core,
  linnik, nonres, pnt, lprime), sanity.py (exact checks and the q < 100 sieve), lsanity.py (L'/L(1, χ) check), xcheck.py (zero sums
  recomputed), ksanity.py (kernel lemmas, inequality (1.3)), logs/ (cert, sanity, lsanity, xcheck, ksanity outputs).
