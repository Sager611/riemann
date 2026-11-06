NO IDEA: mean-periodicity (Fesenko–Ricotta–Suzuki) encodes only meromorphic continuation plus functional equation, so Davenport–Heilbronn and Epstein functions pass it; its one RH-sensitive add-on, the single-sign property of the boundary function, is Landau's theorem with a double real pole: already in the literature, blind on every computable window (DH's off-line zero flips the sign only below x ≈ 10^-275) and false in the literal F_q analogue although RH holds there; total positivity is Gröchenig's PF criterion (N6) and Nyman–Beurling is spectral synthesis in H², where synthesis holds (no Malliavin barrier) (nodes N13, N1, N5, N7)

# explore14 — mean-periodic functions, spectral synthesis, Pompeiu, integral geometry → RH tree

Orca exploration, 6 Nov 2026, about an hour of reading and compute. All numerics in this directory (`runs.txt` has the verbatim output).

## Survey (one page)

- **Mean-periodic functions.** f is X-mean-periodic if its translates do not span X, that is, f ∗ μ = 0 for some μ ≠ 0 in X*. Delsarte 1935; Schwartz 1947 (Ann. Math. 48): in C(R) every mean-periodic f is a limit of exponential polynomials whose exponents are its spectrum, which equals the poles of its Carleman transform; Kahane 1959 (Tata lectures): Fourier–Carleman series. Synthesis fails in C(R²) (Gurevich 1975) and in L¹(G) for every non-compact LCA group G, R included (Malliavin 1959). It holds in C^∞_exp(R) (Gilbert 1970, cited by FRS).
- **Fesenko–Ricotta–Suzuki** (arXiv 0803.2821, Ann. Inst. Fourier 2012). For Z in a class 𝓕 of quotients of completed zeta functions, the boundary function h(x) = f_Z(x) − ε x^{-1} f_Z(1/x), with f_Z the inverse Mellin transform, is C^∞_poly(R₊^×)- and S(R₊^×)*-mean-periodic. This happens if and only if Z has continuation and functional equation ("extension of Hecke–Weil"), and the spectrum of h is then the pole set of Z. For Λ itself the boundary function is trivial: h = x^{-1} − 1, with spectrum {0, 1} (Poisson summation, Tate's boundary term). Zeros of zeta functions enter only through quotients: (a) §5.1, Z_K(s) = Λ_K(2s)Λ_K(2s−1)/Λ_K(s); (b) §5.3, Epstein families Ê(τ,s)/Λ(s) and Λ(2s)Λ(2s−1)/Ê(τ,s); (c) arithmetic surfaces, ζ_ℰ(s)² = (ζ(s)ζ(s−1)/L(E,s))²·(…).
  - Prop. "mono" (Landau, via Widder II §5): if h has one sign near 0 and no real poles in (1/2, ∞), then all poles of Z lie on Re s = 1/2.
  - Prop. "ss": RH, plus simple zeros, plus Σ_{γ<T} |ζ_K(1/2+iγ)|^{-1} ≪ e^{AT} imply that h_K has one sign near 0. The hypothesis is printed with ζ_K; ζ_K′ is meant.
  - Their closing remark: Epstein functions have off-line zeros (Davenport–Heilbronn), and comparing the single-sign behaviour of the two families is left open.
- **Fesenko** ("Analysis on arithmetic schemes II", 2006, J. K-Theory 2010, Thm 52; Moscow Math. J. 2008): one sign of the fourth logarithmic derivative Z_ℰ of the two-dimensional boundary function h_ℰ, plus no real zeros of L(E,s) in (1,2), implies GRH for L(E,s).
- **Suzuki** (arXiv math/0703052): K-Bessel expansion of the boundary term. Assume GRH, simple zeros and Σ|L′(ℰ,ρ)|^{-2} = O(T). Then Z_ℰ(x) = −C x (log 1/x)^{2r+2J+1}(1+o(1)), so the sign holds if r + J ≥ 1, where r is the rank and J counts the reducible fibres of the model. For r = J = 0 one gets −(C + v(x)) x log(1/x) with v bounded and oscillating, so the sign is a quantitative question.
- **Follow-ups**: Suzuki 0805.4547 (GL(2) cusp forms); Oliver 1307.6706 (automorphicity versus mean-periodicity) and 1311.6964 (adelic mean-periodicity on arithmetic surfaces). Search: arXiv API phrase queries, web_search (arxiv), and three GXL evidence queries. None found any 2013–2026 work on the single-sign property, or any link from spectral synthesis, Pompeiu or Radon range conditions to RH, beyond the items below.
- **Synthesis and Nyman–Beurling**: Nyman 1950 (thesis: closed ideals of L¹(R₊) are standard, plus his criterion); Beurling 1955; Bercovici–Foias 1984; Nikolski 1995 (Ann. Inst. Fourier 45); Burnol math/0001013; Manzur–Noor–Santos 2111.05398 (the closure of the NB span is shift-invariant if and only if RH). Connes 1999 and Meyer 2005 (Duke) realise the zeros as a spectrum of the scaling action, with synthesis in nuclear spaces.
- **Total positivity**: Schoenberg 1951 → Gröchenig 2007.12889: RH ⟺ (1/2π)∫e^{−ixτ}/ξ(1/2+τ)dτ is a Pólya frequency function. Already in N6 (2611.02944 PF order k*; R211, R313, R332).
- **Pompeiu / integral geometry**: Brown–Schreiber–Taylor 1973 (the Pompeiu property holds exactly when the transforms of the moved sets have no common zero). We found no arithmetic link. Helgason–Ludwig was dismissed at scout level (R345, R362); Lax–Phillips is R5.5 and R58.

## Transfer 1 — single sign of the FRS boundary function (N13/N1 → really a Landau criterion)

**Closed form (ours, routine).** Take Z(s) = Λ(2s)Λ(2s−1)/Λ(s), with Λ(s) = π^{-s/2}Γ(s/2)ζ(s). Z(1/2+u) is even in u, so the residue a_{−1} at s = 1/2 vanishes and the residues at the zeros are purely imaginary. On the line, Λ(2ρ)Λ(2ρ−1) = |Λ(1+2iγ)|² > 0 and Λ′(ρ) = −iΞ′(γ), with Ξ(t) = Λ(1/2+it) real. Assuming RH with simple zeros:

  x^{1/2} h(x) = a·t − 2 Σ_{γ>0} b_γ sin(γt),   t = log(1/x),   a = −1/(4Λ(1/2)) = 0.0628619872,   b_γ = |Λ(1+2iγ)|²/Ξ′(γ).

b_γ alternates in sign along simple zeros, and |b_γ| ≈ e^{−3πγ/4}/|ζ′(ρ)| up to powers of γ.

**Test (ii)** (`frs_zeta.py`): the definition (inverse Mellin on Re s = 2) was compared with this Fourier–Carleman expansion at 8 points x ∈ [0.02, 4]. With 10 zeros they agree to ≤ 4·10^{−31}. The zeros' part, about 10^{−13}, is resolved (b_{14.13} = 2.85·10^{−14}). The oscillation is at most 2Σ|b_γ|γ/a = 1.3·10^{−11} of the main term for every t > 0, so h > 0 on (0,1) and h < 0 on (1,∞). This uses the first 10 zeros; later zeros add terms of size e^{−3πγ/4}/|ζ′(ρ)|, below 10^{−40} unless |ζ′(ρ)| is extraordinarily small.

**Adversarial test (iii-a), Davenport–Heilbronn** (`dh.py`, `dh2.py`). Set Λ_f(s) = (5/π)^{s/2}Γ((s+1)/2)f(s); the functional equation checks to 10^{−41}. Template: Z_f = s(1−s)Λ(2s)Λ(2s−1)/Λ_f(s). This is the same shape as before; the factor s(1−s) removes the poles at 0 and 1 that 1/Λ cancels in the ζ case.
- Λ_f ∈ [1.136, 1.161] at 39 sample points of (0,1), so there are no real poles and FRS's hypotheses hold. The main coefficient is a_f = −0.0550.
- A scan with step 0.1 finds 43 on-line zeros in (0, 86); a missed close pair would change nothing visible, since its residues would be below 10^{−80}. Their residues alternate, and the oscillation is at most 0.025 of the main term.
- The definition matches the expansion to 10^{−31} at x = 0.1, 0.5, 3.
- The off-line zero ρ* = 0.808517182456637 + 85.6993484853776i has residue |c*| = 2.0·10^{−84}. It beats the main term only for t > 633.9, that is **x < 10^{−275.3}**. So DH has mean-periodicity (FRS Thms 3.2–3.3 apply: Z_f has the functional equation and convexity bounds) and a single sign on all of [10^{−275}, 1).
- Rule: an off-line zero β+iT appears only below t ≈ (3πT/4)/(β−1/2). For ζ above the verified height 3·10^{12}, with β ≤ 0.8734, that means x < 10^{−8.3·10^{12}}.
- As a zero detector, checking the sign to depth t costs about 0.65t digits of cancellation in the definition. It sees only zeros below height ≈ 4(β−1/2)t/(3π), which is exponentially worse than computing zeros.

**Adversarial test (iii-b), known equivalent.**
- This is Landau's theorem for one-signed Laplace transforms plus Ingham-type dominance: the double real pole beats simple complex poles.
- Single sign ⟹ RH is FRS Prop. mono; RH + simple zeros + a weak Gonek–Hejhal bound ⟹ single sign is FRS Prop. ss.
- Pólya's Liouville conjecture uses the same mechanism and fails because its real pole is only simple.
- The criterion is multiplicity-sensitive. Under RH, a zero of order ≥ 3 kills the sign (a t^{m−1} oscillation beats a·t), while Weil's criterion (N1) ignores multiplicity.
- The Fesenko–Suzuki elliptic version is the same, with real pole order 2 + 2r + 2J.

**Sanity test (i), F_q, where RH is a theorem** (`fq.py`, q = 5; u = q^{−s}; the boundary function becomes the coefficient sequence [uⁿ]):
- *Honest template* G = 1/((1−q^{1/2}u)² Z_C(u)): one sign for large n ⟺ RH for C (simple zeros). Constant sign for n = 40..90 for E: y² = x³+x+1/F_5 (#E = 9, P = 1+3u+5u²) and for an on-circle genus-2 polynomial. Sign changes for the functional-equation-symmetric fakes 1+5u+5u² (real roots off the circle) and roots 3e^{±i}, (5/3)e^{±i}. The test passes, but trivially: these are finite exponential sums.
- *Literal FRS template* Z_C(u²)Z_C(qu²)/Z_C(u) = P(u²)P(qu²)/[(1+u)(1+qu)(1−qu²)²P(u)] has **ghost poles** at u = −1 and u = −1/q, that is s = πi/log q and s = 1 + πi/log q. They lie off the line and come from no zero. The pole at −1/q dominates, so the signs alternate for every curve although RH holds. The number-field template works only because ζ(2s) has a single pole on Re s = 1/2.
- *Fesenko's surface analogue* (rational elliptic surface: Z = 1/((1−u)(1−qu)^ρ(1+qu)^m(1−q²u)), squared, abscissa pole removed). One sign holds when the Picard number ρ over F_q beats the multiplicity m of the eigenvalue −q:
  - split surface (ρ = 10): one sign;
  - constant quadratic twist (ρ = 2, m = 8): alternating;
  - ρ = m = 5: decided by constants.
  RH holds in all three. Over F_q the sign is RH plus a Galois-multiplicity condition, a consequence of Weil/Deligne and never a route to it.

**Side observation (ours, routine) on FRS's open comparison of the Epstein families.** For the first family, Z_Q(τ,s) = Ê(τ,s)/Λ(s) has no real singularity at all: Λ has no real zeros, and 1/Λ cancels Ê's poles at 0 and 1. By the same Landau argument as Prop. mono, h_Q(τ,·) therefore changes sign in every interval (0, x₀) unless Z_Q(τ,·) is entire, in which case h_Q(τ,·) ≡ 0 (shift the contour), as at the CM points of their remark. So the single-sign hypothesis of their Prop. 5.3 for this family is never met. For the second family, Λ(2s)Λ(2s−1)/Ê(τ,s), the double real pole makes the sign hold down to x ≈ e^{−(3πT/4)/(β−1/2)}, where β+iT is the lowest off-line zero of Ê(τ,·); for DH this is the 10^{−275} computed above.

**Verdict.** Known (Fesenko 2006/08, Suzuki 2007, FRS 2008). It is RH-equivalent up to simplicity of zeros, there is no mechanism that would prove the sign, and DH passes on every computable window. Fesenko's hope, a sign coming from the boundary of the 2D adelic object, has seen no progress, and has no F_q precedent.

## Transfer 2 — mean-periodicity in a smaller space (Beurling/Carleman growth; N13, N10)

The spectrum of a mean-periodic function lies in the strip set by its exponential growth rates: the two halves of the Carleman transform converge outside that strip, and bounded functions have real spectrum. Let H(t) = e^{−t/2}h(e^{−t}), which is odd. Then RH ⟺ H(t) = O_ε(e^{ε|t|}) for every ε, that is, H is mean-periodic in the subexponential class (FRS: the exponential class ⟺ the functional equation).
- ⟸: the Laplace transform converges for Re w > 0 and equals Z(1/2+w) minus an entire function.
- ⟹: under RH, 1/ζ(1/2+ε+iy) ≪ |y|^δ bounds the shifted contour integral.
- DH: H_f grows like e^{0.3085|t|}, so it fails, correctly.

**Verdict:** Littlewood's M(x) = O(x^{1/2+ε}) criterion in new dress. Routine, known content.

## Transfer 3 — positivity of the spectrum's generator (Bochner, convexity, total positivity; N1, N6)

- **Positive definiteness** (Bochner) of the even function H′ would need weights −b_γγ ≥ 0. These alternate in sign under RH, so H′ is not positive definite even if RH holds; the requirement is too strong.
- **Convexity** fails: the oscillating part has mean zero.
- **Weil's form:** for Z = Λ′/Λ the residues are the multiplicities, and the FRS function is Weil's distribution Σ_ρ e^{(ρ−1/2)t}. Its spectral synthesis is the explicit formula, and its positive definiteness is Weil's criterion (N1).
- **Total positivity** (Schoenberg's PF_∞) is the only positivity that forces a spectrum onto a line. That is Gröchenig 2007.12889, and its finite orders are N6's PF_k.
- **Annihilator side:** h_Q is annihilated by Riemann's positive kernel Φ, which is not PF_∞ (R7.6).

**Verdict:** known.

## Transfer 4 — Nyman–Beurling as spectral synthesis / Pompeiu in R₊^× (N5)

Nyman–Beurling (NB) asks whether ρ(1/x) is cyclic under the dilation semigroup on L²(0,1). Under the Mellin transform this is the shift semigroup on H²(Re s > 1/2). By Beurling–Lax the closed invariant subspaces are ΘH². The hull is the set of off-line zeros, and the Dirichlet series (→ 1 as Re s → ∞) excludes a singular inner factor. So synthesis holds there (Nyman 1950, Beurling 1955, Bercovici–Foias 1984, Nikolski 1995, Burnol 2001, 2111.05398).
- Malliavin's failure of synthesis needs two-sided L¹. Closed ideals of L¹(R₊) are standard (Nyman 1950), so there is **no Malliavin-type barrier and no new criterion**.
- The "Pompeiu" form is just the annihilator form of NB. The Pompeiu problem itself is trivial in one variable: (b^s − a^s)/s has zeros, so the property always fails. The rotations that make it non-trivial in Rⁿ have no counterpart in R₊^×.
- Two-generator variant (dilates of ρ and of its χ-twist): dense ⟺ ζ and L(s,χ) share no off-line zero. This is weaker than RH but just as out of reach. The operator side was done in explore13.

## Transfer 5 — support theorems / Radon / Helgason–Ludwig (N7)

Support theorems (Paley–Wiener–Schwartz, Helgason) trade support for exponential type. Λ has order 1 and maximal type, so there is nothing to support. Zero-free half-planes are growth statements about 1/ζ (Transfer 2, N10), not support statements. Range conditions were already dismissed (R345, R362), and scattering is R5.5 and R58. **Dead.**

## Ranking (RH relevance × chance × gain)

| Transfer | Node | RH relevance | Chance | Gain if it worked | Status |
|---|---|---|---|---|---|
| 1 single sign (FRS/Fesenko) | N13/N1 | implies RH | ≈ 0: DH passes to 10^{−275}; literal F_q analogue fails | none beyond Landau | known 2006–08 |
| 2 subexponential mean-periodicity | N13/N10 | ⟺ RH | 0 | none (Littlewood) | routine |
| 3 positivity: Bochner / convexity / PF | N1/N6 | Bochner and convexity false under RH; PF ⟺ RH | low | already N6 | known (Gröchenig) |
| 4 NB as synthesis / Pompeiu | N5 | ⟺ RH | 0 | none (synthesis holds) | known since 1950–95 |
| 5 support theorems | N7 | none | 0 | none | dead |

Honest trade-off: mean-periodicity is the right language for "functional equation ⟹ spectrum" (it is DH-blind by design). Every RH-sensitive addition we found (sign, growth, total positivity, synthesis hull) turns out to be an existing criterion: Landau, Littlewood, Gröchenig/N6 or Beurling–Lax.

## Suggested NOGO line

- explore14 (Orca, exploration: mean-periodic functions, spectral synthesis, Pompeiu, integral geometry → N1, N5, N7, N13) NO IDEA: Fesenko–Ricotta–Suzuki mean-periodicity encodes only continuation + functional equation (DH and Epstein pass). Its single-sign criterion is Landau's theorem with a double real pole, already published (Fesenko, Suzuki, FRS). It is blind on every computable window: DH's 0.8085+85.70i flips the sign only below x ≈ 10^-275, and zeta's oscillation is 1.3e-11 of the main term. It is false in the literal F_q analogue although RH holds (ghost pole at s = 1+πi/log q). Smaller-space mean-periodicity is Littlewood, total positivity is Gröchenig/N6, and NB is synthesis in H² with no Malliavin barrier [ExploreMeanPeriodic](Riemann/Deep/ExploreMeanPeriodic/REPORT.md)

## Files

- `REPORT.md`: this report.
- `frs_zeta.py`: FRS Z for ζ; Laurent data, residues, definition vs expansion.
- `dh.py`: DH functional-equation check, real values, off-line zero, residue, sign-flip scale.
- `dh2.py`: DH on-line zeros < 86, residues, definition vs expansion.
- `fq.py`: F_q sanity tests (honest, literal FRS and Fesenko-surface templates).
- `runs.txt`: verbatim outputs, run times 259 s, 166 s, 384 s and < 2 s; single-threaded, mpmath via uv.

RH relevance: N13 and N1 (single-sign/positivity criteria), N5 (synthesis), N7 (support); no node advanced.
