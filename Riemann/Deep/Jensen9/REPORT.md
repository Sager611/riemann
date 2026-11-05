RECORD: J^{d,n} hyperbolic for all n and all d ≤ 16,001,870,271,185,050,232,604,444 ≈ 1.6002·10^25 (was 9,001,052,027,541,590,755,840,000 ≈ 9.0·10^24, Kim–Lee 2021 / GORTTW Thm 1.2 with Platt–Trudgian; our TREE said 8)

# jensen9 report (node N6, Jensen polynomials of ξ), 2026-11-05

## Short answer

- **The brief's premise was wrong.** GORZ's "d ≤ 8 for all n" (2019, a proof sketch) was already beaten by
  Chasse (thesis 2011): d ≤ 2·10^17 for all n. The same method with Platt–Trudgian's verified height
  T = 3,000,175,332,800 gives d ≤ T² ≈ 9.0·10^24 (Kim–Lee 2021; GORTTW 2022 Thm 1.2; Dumas 2026). Details: SOTA.md.
- **New, a routine corollary:** feeding OpenAI's 2026 zero-free half-plane Re s > 7/8 into that argument narrows the
  zero sector by a factor 3/4 and multiplies the range by 16/9: **D = (16T² − 4)/9 ≈ 1.6002·10^25.**
  d = 9 is covered, as is every degree up to D.
- **Small computational extra:** FLINT puts the first zero above T on the critical line (16 s), so RH holds up to
  T' = 3,000,175,332,800.2865 and the range rises to D' = 16,001,870,271,188,106,411,210,124.
- **Explicit N(d):** none is published for any degree above the all-n range. For d = D + 1 every shift is open.

## Theorem

Let 8ξ(1/2+z) = (−1+4z²)Λ(1/2+z) = Σ γ(n) z^{2n}/n! and J^{d,n}(X) = Σ_{j≤d} C(d,j) γ(n+j) X^j. For every n ≥ 0
and every 1 ≤ d ≤ ⌊1/sin²δ⌋ with δ = 2 arctan(3/(8T)), J^{d,n} has only real zeros. Here
⌊1/sin²δ⌋ = ⌊16T²/9 + 1/2 + 9/(256T²)⌋ = (16T² − 4)/9 = 16,001,870,271,185,050,232,604,444.

## Proof

1. Θ(w) := 8ξ(1/2 + √w) = Σ γ(k) w^k/k! is real entire of order 1/2 with positive coefficients. Its zeros are
   w = z² for the zeros z = ρ − 1/2 of ξ(1/2+z).
2. If |Im z| ≤ T then Re z = 0 (Platt–Trudgian; conjugation covers negative heights; ξ has no real zeros), so w = −γ².
   If |Im z| = y > T then |Re z| = |x| ≤ 3/8 (no zero with Re s > 7/8, plus the functional equation), so
   |arg w − π| = 2 arctan(|x|/y) < δ = 2 arctan(3/(8T)) < π/2. All zeros of Θ lie in the closed sector S(π,δ).
3. The sector is convex, so every derivative Θ^(n) (genus 0) has its zeros in S(π,δ) too (Gauss–Lucas; Chasse thesis
   Thm 135). Θ^(n) has order 1/2 and positive coefficients γ(n+k)/k!.
4. J^{d,n} is the degree-d Jensen polynomial of Θ^(n). Chasse's sector lemma (thesis Thm 171 = paper Thm 3.6,
   from Obreschkoff 1961): order ≤ 1, positive coefficients, zeros in S(π,δ) ⇒ Jensen polynomials of degree
   d ≤ 1/sin²δ have only real zeros.
5. With u = tan(δ/2) = 3/(8T): 1/sin²δ = (1+u²)²/(4u²) = 16T²/9 + 1/2 + 9/(256T²). T ≡ 5 (mod 9), so
   16T² = 9q + 4 and the bound is q + 17/18 + 9/(256T²); its floor is q = (16T² − 4)/9. ∎

The old bound is the case |x| < 1/2: 1/sin²δ = T² + 1/2 + 1/(16T²), floor T².

## Inputs

| input | statement used | status |
|---|---|---|
| Platt–Trudgian, Bull. LMS 53 (2021), arXiv 2004.09765, Thm 1 | RH true up to height 3,000,175,332,800 (lowest 12,363,153,437,138 zeros) | refereed; N(T) recomputed by FLINT (part C) |
| OpenAI, "The Quasi-Riemann Hypothesis: a zero-free half-plane Re s > 7/8" (preprint, 30 Sep 2026) | ζ(s) ≠ 0 for Re s > 7/8 | unrefereed, but its Lean proof builds sorry-free here (Riemann/Deep/QRH, axioms propext, Classical.choice, Quot.sound) |
| Chasse, PhD thesis U. Hawaii 2011, Thm 171 and 135; Complex Var. Elliptic Equ. 58 (2013) 875–885 | sector lemma; closure under differentiation | refereed; used the same way by GORTTW (proof of Thm 1.2) and Kim–Lee (Thm 4) |
| Obreschkoff, Acta Math. Acad. Sci. Hungar. 12 (1961) | finite-degree version behind Thm 171 | classical |

## Certificates (Riemann/Deep/Jensen9/verify.py; 29 s, peak 77 MB, one thread)

| part | check | result |
|---|---|---|
| A | exact rational computation of D, the old bound T², ratio 16/9, conditional and Ξ-family variants | D = 16001870271185050232604444 |
| B | sector geometry: edge zero 3/8 + iT, 2000 random admissible zeros inside the sector | pass |
| C | FLINT/Arb N(3,000,175,332,800) by Turing's method | 12,363,153,437,138 = Platt–Trudgian's count (10 s) |
| D | sector lemma spot checks with exact integer polynomials and certified roots: sharpness (pair −5 ± i: J^26 real-rooted, J^27 not), 600 Jensen polynomials of random sector polynomials, Gauss–Lucas | pass |
| E | rigorous γ(n) (Arb power series of ξ, 4400 bits) agree with .scratch gamma.txt; J^{9,n} certified for 0 ≤ n ≤ 150 and J^{d,0} for 2 ≤ d ≤ 40 by sign alternation at exact dyadic points | pass (consistency only; the theorem covers them) |
| F | FLINT isolates zero #12,363,153,437,139 = 1/2 + i·[3000175332800.28652 ± 4.3e-6] on the line | D' = 16001870271188106411210124 |

Parts C–F need python-flint (pip install python-flint; we used 0.9.0 with FLINT 3.6.0); A–B run on plain Python.

## Side results (same argument)

- **Ξ-family** (Jensen polynomials of Ξ^(n), Ξ(t) = ξ(1/2+it)): Kim–Lee Thm 4 with strip half-width 3/8 instead of 1/2
  gives d ≤ 1 + 64T²/9 = 64,007,481,084,740,200,930,417,778 ≈ 6.4·10^25 (was 1 + 4T² ≈ 3.6·10^25, Kim–Lee).
- **Conditional half-planes:** Liu, arXiv 2610.12234 (0.874957…): D ≈ 1.60055·10^25; Chang–Wang, arXiv 2611.00924
  (0.87446205769353, conditional on cited inputs): D ≈ 1.60479·10^25.

## The original d = 9 plan

Moot: d = 9 lies far inside the all-n range. What the plan would have needed: GORTTW's c and Kim–Lee's n₁ are not
computed, so an explicit N(9) would have to be built from scratch (sharp Laplace-method bounds for γ(n+j)/γ(n)).
The ball-arithmetic pipeline it asked for is in verify.py part E (γ(n) for n ≤ 160, J^{9,n} for n ≤ 150, 3 s).

## Gap above the record

- d = D + 1: no shift n is proved. Large-n results (GORTTW c·e^{d/2}; Kim–Lee O(d^{1/2+ε}); Holland
  n³log²(n+2) ≥ Kd^5) all have unknown constants, so no explicit N(D+1) exists.
- The only working route is a narrower sector: either a higher verified RH height (each extra unit of height adds
  (32/9)T ≈ 1.07·10^13 to D; part F adds one zero) or a better zero-free strip (D ≈ T²/(4a²) for |Re ρ − 1/2| ≤ a).
  The lemma is sharp for its class (part D), so the constant 1/sin²δ cannot be improved without more information on
  the zeros above T.
- Heuristic only (not tested): the many real zeros below T are unlikely to help, since in the heat-flow picture real
  zeros push a nearby complex zero away from the axis. I did not pursue a sharper lemma.

## Review and novelty

- Own verification: verify.py (above). Reviewer (subagent, verdict only): "CORRECT — conditional on I1–I4 …
  routine corollary … increasing the degree range by approximately 16/9; no new proof mechanism."
- Novelty: arXiv listings math.NT, math.CV, math.CA for 2026-10 and 2026-11 (titles: no Jensen paper using 7/8);
  Zenodo "Jensen polynomials" (102 records; all-n statements use 9·10^24 or less); GXL evidence API, 2 claims (hits:
  GORTTW, Kim–Lee only); our repo (R103–R107 list 7/8 consequences, none on Jensen). Not found elsewhere.
- Classification: routine deduction from the 7/8 half-plane → COROLLARIES.md, not DISCOVERIES.md. The bigger
  practical change is the corrected N6 record.

## Suggested repo edits (orchestrator)

- TREE N6 record: "d ≤ 1.6002·10^25 for all n (ours, Chasse/Kim–Lee sector lemma + Platt–Trudgian + 7/8;
  before 9.0·10^24, Kim–Lee 2021 / GORTTW 2022; GORZ's d ≤ 8 was a sketch); large n: n ≫ d^{1/2+ε}
  (Kim–Lee, ineffective)". Barrier: "the all-n range is ≈ T²/(4a²); no explicit N(d) beyond it."
- COROLLARIES.md: one line for D = (16T² − 4)/9 and the Ξ-family 6.4·10^25.

## Files

- Riemann/Deep/Jensen9/SOTA.md — state of the art with citations, tables, search log.
- Riemann/Deep/Jensen9/verify.py — certificates A–F (29 s).
- .scratch/jensen9/STATE.md — agent state; .scratch/jensen9/gamma.txt — the orchestrator's γ(0..40) table.
