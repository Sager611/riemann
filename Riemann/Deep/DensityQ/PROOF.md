# PROOF: zero density for families of Dirichlet L-functions just below the Lean-checked half-plane

Nodes N7 (zero density) and N11 (L-function analogues). Certificate: [verify.py](verify.py) (exact rationals, standard
library only, ALL PASS). Paper proofs, not refereed, no Lean.

## 0. Setting, notation, sources

B := 0.87338196661320752. (H): every Dirichlet L-function is non-zero in Re s > B (Lean-checked, see TREE.md N7;
Bebel's newer B_2 + 2.2·10^−254 < B only helps).

Families. (F_q) all characters mod a fixed q, X := qT. (F_Q) primitive characters of conductor ≤ Q, X := Q²T. Pairs
(χ, t), |t| ≤ T, T ≥ 2, are *well spaced* if (χ, t) ≠ (χ', t') implies χ ≠ χ' or |t − t'| ≥ 1. For (F_Q) put
λ := log(QT)/log(Q²T) ∈ [1/2, 1] (the largest analytic conductor q(|t|+1) is X^{λ+o(1)}); for (F_q), λ = 1.

Exponents (the family version of Tao–Trudgian–Yang's language, arXiv 2501.16779): a polynomial D(χ, t) =
Σ_{N<n≤2N} a_n χ(n) n^{−it}, |a_n| ≤ N^{o(1)}, has length N, and τ := log X / log N. LV(σ, τ) is the least exponent
with R ≪ N^{LV(σ,τ)+o(1)} for every well-spaced family set of R points with |D| ≥ N^σ. The Montgomery conjecture (MC)
at (σ, τ) is LV(σ, τ) ≤ 2 − 2σ. A(σ) for a family is the least A with Σ N(σ, T, χ) ≪ X^{A(1−σ)+ε}; A*(σ) for (F_Q),
A_q(σ) for (F_q). Smaller A is better.

Inputs (all published):
- [PY] Petrow–Young, arXiv 1908.10346 (Duke 2023), Thm 1.1: L(1/2 + it, χ) ≪_ε (q(1 + |t|))^{1/6+ε} for every
  primitive χ mod q (checked in the TeX source, line 196).
- [J] Jutila, Acta Arith. 32 (1977) 55–62 (open scan, matwbn.icm.edu.pl/ksiazki/aa/aa32/aa3216.pdf), Theorem, (1.2)
  for (F_q) and (1.3) for (F_Q): for fixed k ≥ 1, R ≪ (GNV^{−2} + X(G⁴N²V^{−8})^k + (XG²V^{−4})^k)X^ε with
  G = Σ|a_n|² ≤ N^{1+o(1)}. With k = 2, V = N^σ: LV(σ, τ) ≤ LV_J(σ, τ) := max(2 − 2σ, τ + 12 − 16σ, 2τ + 4 − 8σ).
  So MC holds for τ ≤ min(3σ − 1, 14σ − 10) = 3σ − 1 when σ ≥ 9/11. Its Corollary (1.5)–(1.6) gives
  A_q, A* ≤ max(2, ...) = 2 (density hypothesis) for σ ≥ 21/26; nothing below 2.
- [Z] Classical zero detection with a short mollifier M_{X^ε}, the weight e^{−n/Y}, class I/II zeros and Huxley's
  powers (Montgomery, Topics, LNM 227, Ch. 12; Jutila [J] §4; written out for (F_q) in Chen–Gupta–Li, arXiv
  2507.08296v2, §5).
- [L] Conductor-uniform Lindelöf beyond B (COROLLARIES, R354): under (H), L(σ + it, χ) ≪_ε (q(|t| + 1))^ε for σ > B
  (Titchmarsh §14.2: Borel–Carathéodory and three circles, with log(q(|t|+3)) in place of log t).

## 1. The reflected bound (under (H))

Lemma 1. Assume (H). For primitive ψ mod k ≥ 2 and fixed η, ε > 0:
L(1 − B − η + iv, ψ) ≪ (k(|v| + 2))^{B+η−1/2+ε}. For k = 1: ζ(1 − B − η + iv) ≪ (|v| + 2)^{B+η−1/2+ε} for |v| ≥ 1.
Proof. Functional equation Λ(s, ψ) = (k/π)^{(s+a)/2} Γ((s + a)/2) L(s, ψ) = ε(ψ) Λ(1 − s, ψ̄), |ε(ψ)| = 1. At
s = B + η − iv, Stirling gives |Γ((s+a)/2)/Γ((1−s+a)/2)| ≍ (|v| + 2)^{B+η−1/2} uniformly (both arguments have real
part ≥ (1 − B − η)/2 > 0). Hence |L(1 − B − η + iv, ψ̄)| ≍ (k(|v|+2))^{B+η−1/2} |L(B + η − iv, ψ)|, and [L] bounds the
last factor by (k(|v|+2))^ε. Replace ψ̄ by ψ. For ζ this is DensityB (C1'). ∎

## 2. Theorem A: a longer Montgomery range for both families (under (H))

τ_B(σ) := (2σ − 2 + B)/(B − 1/2), g_B(τ) := 1 − B + (B − 1/2)τ; g_B(τ_B(σ)) = 2σ − 1 (verify A1).

Theorem A. Assume (H). Let 3/4 < σ < (1 + B)/2, (a_n) 1-bounded on (N, 2N], W a well-spaced set of pairs from (F_q)
or (F_Q) with |D(χ, t)| ≥ N^σ on W, and X ≤ N^τ with τ < τ_B(σ). Then |W| ≪ N^{2−2σ+ε}.

Proof. Take b smooth, b ≥ 1 on [1, 2], supp b ⊂ [1/2, 3], with Mellin transform b̃ decaying faster than any power in
vertical strips. By Cauchy–Schwarz (Halász–Montgomery), with unimodular c_r,
(R N^σ)² ≤ (Σ_r |D_r|)² ≤ (Σ_n |a_n|²/b(n/N)) Σ_{r,s} |K(r, s)|, K(r, s) := Σ_n b(n/N) χ_r(n) χ̄_s(n) n^{−i(t_r−t_s)}.
(i) r = s: |K| ≤ 3N. (ii) χ_r = χ_s, r ≠ s, u := t_r − t_s, |u| ≥ 1: K is the Mellin integral of
ζ(w + iu) Π_{p|q_r}(1 − p^{−w−iu}) N^w b̃(w); moving to Re w = 1 − B − η passes the pole w = 1 − iu, whose residue is
≪_A N(1 + |u|)^{−A}; the rest is ≪ N^{1−B} X^{B−1/2+3η} by Lemma 1 and |Π| ≤ 2^{ω(q_r)}. Summed over s, the residues
give O(N) by the spacing. (iii) χ_r ≠ χ_s: ψ := χ_r χ̄_s is non-principal (two distinct primitive characters, or two
distinct characters mod q, never agree on all n coprime to the moduli). Let ψ* mod k* induce ψ; k*(|u| + 2) ≤ 4X in
both families. K is the Mellin integral of L(w + iu, ψ*) Π_{p | q_r q_s}(1 − ψ*(p)p^{−w−iu}) N^w b̃(w); no pole; on
Re w = 1 − B − η Lemma 1 and |Π| ≤ 2^{ω(q_r q_s)} = X^{o(1)} give |K| ≪ N^{1−B} X^{B−1/2+3η}.
Hence R² N^{2σ} ≪ N(NR + R² N^{1−B} X^{B−1/2+3η}) and N^{2−B} X^{B−1/2+3η} ≤ N^{1+g_B(τ)+3τη}. As g_B is increasing and
g_B(τ) < g_B(τ_B) = 2σ − 1, a small η makes the second term at most N^{2σ}R²/2 for large N, so R ≪ N^{2−2σ}. ∎

Comparison (verify A3–A9): the best published family ranges are Jutila's 3σ − 1 ([J], k = 2) and Huxley's 4σ − 2;
the Petrow–Young line gives τ < 12σ − 9. τ_B(σ) beats all of them exactly on (σ*, s_PY), σ* = (5 − 4B)/(7 − 6B) =
0.856092, s_PY = (10B − 13/2)/(12B − 8) = 0.900522. (Under (H) the full off-diagonal exponent is the lower hull of the
vertices (0, 1/2), (1 − B, B − 1/2), (1/2, 1/6), (B, 0) of the conductor-aspect μ; the (H) vertex is the best one for
1 ≤ τ ≤ (B − 1/2)/(B − 2/3) = 1.806, verify A8.)

## 3. Lemma Z: zero density from large values in the families

For a large-value bound LV (as a function of τ) put F(τ) := min_{k≥1} k·LV(σ, τ/k)/τ (Huxley's powers) and
G(a) := sup_{a≤τ≤2a} F(τ). Since F(τ) ≤ F(τ/m) for every integer m ≥ 1, sup_{τ≥a} F = G(a).

Lemma Z. For 1/2 < σ < 1 and every a ≤ cut := (6σ − 3)/λ: A(σ)(1 − σ) ≤ max(G(a), A_ζ(σ)(1 − σ)), where A_ζ(σ) is any
admissible zeta exponent (N(σ, T) ≪ T^{A_ζ(σ)(1−σ)+ε}).
Proof (standard, [Z]). Take the mollifier length X^ε and Y := (max conductor·T)^{1/(6σ−3)+η} = X^{λ/(6σ−3)+λη}.
A class-II zero gives Y^{1/2−σ} |L(1/2 + iγ', χ)| (qT)^ε ≫ 1 for some |γ' − γ| ≤ log²T; by [PY] |L(1/2 + iγ', χ)| ≪
X^{λ/6+ε}, so for η > 0 fixed and ε small there are no class-II zeros. A class-I zero yields a large value
|Σ_{N<n≤2N} a_n χ(n) n^{−it}| ≫ N^{σ−o(1)} with |a_n| ≤ N^{o(1)}, well-spaced points, and some N ∈ [X^ε, Y log²Y],
that is τ ≥ cut − O(η). Each dyadic N contributes ≤ X^{F(τ)+o(1)} (a large value of D gives one of D^k), and
sup_{τ ≥ cut − O(η)} F = G(cut − O(η)) → G(cut) as η → 0 (F is continuous); for a < cut, G(a) ≥ sup_{τ≥cut} F anyway.
Characters (repaired after review). Imprimitive characters in (F_q) reduce to primitive ones of conductor q* | q at the
cost d(q) = q^{o(1)}: L(s, χ) and L(s, χ*) differ by finitely many Euler factors, which vanish only on Re s = 0. The
detection above is run only for primitive non-principal χ, where L(s, χ) is entire and the Mellin identity has no
residue term (ε(χ) = 0 in Chen–Gupta–Li (5.1)). The one remaining primitive character is the trivial one of conductor
1, i.e. ζ itself, which induces every principal character; its zeros are all zeros of ζ (not O(log²) of them) and are
counted by a zeta density bound: N(σ, T) ≪ T^{A_ζ(σ)(1−σ)+ε} ≤ X^{A_ζ(σ)(1−σ)+ε} since T ≤ X, at most d(q) times in
(F_q) and once in (F_Q). ∎

Remark. Without [PY] (only the fourth moment for class II, as in [J] and Chen–Gupta–Li §5) the scheme stops at A = 2 on
[0.84, B): the class-II count X^{1+ε}Y^{2−4σ} forces Y ≥ X^{1/2} for A ≤ 2, and then polynomials of length just
above X^{1/2} (τ → 2−) cost F(2) = min(LV(σ, 2)/2, 2 − 2σ) = 2 − 2σ, since LV(σ, 2) ≤ 2(2 − 2σ) is the best known
bound there (Montgomery ranges are < 2σ). This is why all earlier family results stop at the density hypothesis.

## 4. Theorem C (unconditional): below the density hypothesis on (5/6, 7/8]

Theorem C. For 5/6 ≤ σ ≤ 7/8, both A_q(σ) and A*(σ) are ≤ 4/(6σ − 3). That is,
Σ_{χ mod q} N(σ, T, χ) ≪ (qT)^{4(1−σ)/(6σ−3)+ε} and Σ_{q≤Q} Σ*_χ N(σ, T, χ) ≪ (Q²T)^{4(1−σ)/(6σ−3)+ε}.
4/(6σ − 3) < 2 and 4/(6σ − 3) < 3/(3σ − 1) for σ > 5/6 (verify E1, E2); at σ = 0.84, 0.86, 0.87, B: 1.9608, 1.8519,
1.8018, 1.7855 (3/(3σ − 1): 1.9737, 1.8987, 1.8634, 1.8517).

Proof. Lemma Z with a := 6σ − 3 ≤ cut (λ ≤ 1) and LV := LV_J from [J] (k = 2), c := 2 − 2σ. On [a, 2a]:
(1) τ ∈ [6σ − 3, 6σ − 2]: τ/2 ≤ 3σ − 1, MC: F ≤ 2c/τ ≤ 2c/a.
(2) τ ∈ [6σ − 2, 9σ − 9/2] (non-empty iff σ ≥ 5/6, verify C1a): k = 2 beyond MC, 2·LV_J(τ/2)/τ =
max(2c/τ, 1 + (24 − 32σ)/τ, 2 + (8 − 16σ)/τ); the first is ≤ 2c/a, the last two increase in τ
(24 − 32σ < 0, 8 − 16σ < 0); at τ = 9σ − 9/2 the second is (39/2 − 23σ)/(9σ − 9/2) ≤ 2/9 (σ ≥ 0.82, C1d) and the
third is exactly 2/9 (C1b); and 2/9 ≤ 2c/(6σ − 3) iff σ ≤ 7/8 (C1c).
(3) τ ∈ [9σ − 9/2, 12σ − 6]: τ/3 ≤ 4σ − 2 ≤ 3σ − 1, MC: F ≤ 3c/τ ≤ 3c/(9σ − 9/2) = 2c/(6σ − 3).
So G(a) ≤ 2c/(6σ − 3) (also checked on 22 rational σ by the generic evaluator, C1).
ζ part: A_ζ(σ) ≤ 3/(2σ) for σ ≥ 3831/4791 = 0.7996 (Ivić, The Riemann Zeta-Function, Thm 11.4; Bourgain, IMRN 2000
no. 3, extends it to σ ≥ 3734/4694; ANTEDB Table "Current best upper bound on A(σ)"), and 3/(2σ) ≤ 4/(6σ − 3) iff
σ ≤ 9/10 (verify F1). Lemma Z gives A ≤ 4/(6σ − 3). ∎
The half-plane is not used. The binding term is the window start a = 6σ − 3: polynomials of length X^{1/(6σ−3)}, i.e.
the L-sums just below the length where [PY] stops excluding large values.

## 5. Proposition B: the class-II barrier, and why (H) gives nothing more for the uniform exponents

Proposition B. Let 0.84 ≤ σ < B, λ = 1, and let the large-value input be MC for τ ≤ τ_M plus Huxley subdivision
(LV ≤ 2 − 2σ + max(0, τ − τ_M)) with τ_M ≤ 8σ − 5. Then G(a) ≥ (4 − 4σ)/(6σ − 3) for every a ≤ 6σ − 3,
with equality at a = 6σ − 3 when τ_M ≥ (11/2)σ − 13/4.
Proof. Every branch satisfies k·LV(σ, τ/k)/τ ≥ kc/τ. For a ∈ [τ_M + c, 6σ − 3]: the k = 1 branch is
(c + a − τ_M)/a ≥ 2c/a, the others are ≥ 2c/a, so F(a) ≥ 2c/a ≥ 2c/(6σ − 3). For a < τ_M + c ≤ 2a the window contains
τ_M + c, where every branch is ≥ 2c/(τ_M + c) ≥ 2c/(6σ − 3) (as τ_M + c ≤ 6σ − 3 iff τ_M ≤ 8σ − 5). For 2a < τ_M + c,
F(a) ≥ c/a > 2c/(τ_M + c). Equality: at a = 6σ − 3 the k = 2 branch is 2c/a (a/2 ≤ τ_M), the local maximum
(3 − 3σ)/(τ_M + 1 − σ) inside the window is ≤ 2c/a iff τ_M ≥ (11/2)σ − 13/4, and F(2a) ≤ 4c/(2a) as 2a/4 ≤ τ_M.
(verify C2a–C2c: minimum over windows at five σ with τ_M = τ_B.) ∎

Consequences.
- Under (H), Theorem A gives τ_M = τ_B with (11/2)σ − 13/4 ≤ τ_B ≤ 8σ − 5 on [0.84, B] (C2a, C3a): the half-plane
  reaches the cap 4/(6σ − 3) without Jutila's theorem, and cannot go below it in this scheme.
- (H) does not move the cut: under (H) and [PY] the conductor-aspect exponent of an L-sum of length N at X = N^τ is
  min(τ/2, g_B(τ), 1/2 + τ/6, B), equal to 1/2 + τ/6 on [2, 6(B − 1/2)] ⊃ [2, 6σ − 3] (B1, B2); the interpolated
  lines μ*(u) ≤ (B − u)/(6(B − 1/2)) give their best cut at u = 1/2, namely 6σ − 3 again (B3, B4); interpolated moments
  between the fourth moment at 1/2 and Lindelöf at B are weaker than the fourth moment itself (DensityB Prop. 2(b),
  same algebra).
- Without t-subdivision (T tiny compared with q) the (H) route gives max(4/(6σ − 3), 3/τ_B) (C3b), and Theorem C
  (which needs no subdivision) is better.
- What would break the barrier: an L-sum large-value estimate at lengths ≈ X^{1/(6σ−3)}, e.g. a hybrid twelfth moment
  Σ_χ ∫_{−T}^{T} |L(1/2 + it, χ)|^{12} dt ≪ X^{2+ε} (known only in the q-aspect for smooth or prime-power moduli:
  Nunes arXiv 1802.05180, Milićević–White arXiv 1908.04833), or a hybrid μ*(1/2) < 1/6.

## 6. Theorem D (under (H)): the large-sieve family when T is a moderate power of Q

τ0(σ) := τ_B(σ) + 1 − σ; λ_b(σ) := 3(6σ − 3)/(4τ0); λ_1(σ) := (1 + (1 − σ)/τ0)/2;
λ_a(σ) := max((6σ − 3)/(4σ), λ_1(σ)).

Theorem D. Assume (H). For σ* ≤ σ < B and λ_1(σ) ≤ log(QT)/log(Q²T) ≤ λ_b(σ),
Σ_{q≤Q} Σ*_χ N(σ, T, χ) ≪ (Q²T)^{3(1−σ)/τ0(σ) + ε}, 3/τ0 = 3(B − 1/2)/((5/2 − B)σ − (5/2 − 2B)).
In T (T = Q^{(2λ−1)/(1−λ)}): the bound holds for Q^{0.18} ≤ T ≤ Q^{13.8} at σ = 0.86 and Q^{0.15} ≤ T ≤ Q^{14.7} at
σ → B; the comparison below is claimed for λ ≥ λ_a, i.e. T ≥ Q^{0.69} (σ = 0.86) to T ≥ Q^{0.79} (σ → B) (verify table D).
Proof. Lemma Z with a = 4τ0/3 ≤ cut (this is λ ≤ λ_b). MC holds for τ < τ_B (Theorem A); t-subdivision gives
LV ≤ max(2 − 2σ, 3 − 3σ + τ − τ0) on [2τ0/3, τ0], which needs T ≥ X^{(1−σ)/τ0}, i.e. λ ≥ (1 + (1 − σ)/τ0)/2; this is
≤ (3 − 3σ)τ/τ0 there (τ0 ≥ 3 − 3σ), and powers k = 2, 3 cover [4τ0/3, 8τ0/3] (Tao–Trudgian–Yang Cor. zero-large-cor3),
so G(4τ0/3) ≤ (3 − 3σ)/τ0. ζ part: DensityB Theorem 2 (same hypothesis (H)) gives A_ζ(σ) ≤ 3/τ0 on [σ*, B), and
T ≤ X. ∎
Same method without (H) (MC up to 3σ − 1 from [J], with or without subdivision): for λ ≥ (6σ − 3)/(4σ) every window
contains τ = 4σ (where F = (3 − 3σ)/(2σ)) or a worse point, so the best is exactly 3/(2σ); the comparison is claimed
for λ ∈ [λ_a, λ_b] only (verify D2; for smaller λ both methods improve and are not compared). Gain: 3/(2σ) → 3/τ0, from 0 at σ* to 2.3% at B
(1.7442 → 1.7349 at 0.86; 1.7175 → 1.6783 at B−). Both are below the published 2 ([J]); the comparison bound 3/(2σ)
in this regime is itself unpublished (it needs [PY]).

## 7. Self-review

Pass 1 (hypotheses, ranges, directions).
- Theorem A: Lemma 1 needs only Re s = 1 − B − η > 0 Euler factors and Stirling; the pole of ζ appears only for
  χ_r = χ_s and is summable by the spacing; conductors of χ_r χ̄_s are ≤ q (F_q) or ≤ Q² (F_Q); |u| ≤ 2T. Range
  3/4 < σ < (1 + B)/2 keeps 1 ≤ τ_B ≤ 2 (A2), though only σ ∈ [σ*, B) is used for density.
- Lemma Z: the class-II exclusion uses [PY] at heights |γ'| ≤ T + log²T for primitive characters; q = 1 is the
  classical Weyl bound. The class-I lengths start at τ = cut − O(η); F is continuous, so η → 0 is harmless. The
  coefficients after zero detection are divisor-bounded, as allowed by [J] (G ≤ N^{1+o(1)}).
- Theorem C uses [J] exactly as stated (k = 2, any T ≥ 2, both families), no subdivision.
- Theorem D uses subdivision in t only (sub-families keep the same conductor bound), hence the lower end λ_1.
- Direction: smaller A is better; 4/(6σ − 3) < 2 iff σ > 5/6, < 3/(3σ − 1) iff σ > 5/6 (E1, E2).
Pass 2 (newer records, unproved inputs, what the half-plane really adds).
- Old record (after review, conservative): Chen–Gupta–Li (arXiv 2507.08296v2, TeX l. 137–148, eq. zeeshmh) state that
  the Halász–Montgomery–Huxley estimate (Huxley 1975) gives Σ_{χ mod q} N(σ, T, χ) ≤ (qT)^{3(1−σ)/(3σ−1)+o(1)}, below 2
  for σ > 5/6; we take 3/(3σ − 1) as the old record for both families. (Jutila [J] p. 57 lists only density-hypothesis
  results for families, σ ≥ 5/6 and his σ ≥ 21/26; in the scheme of Lemma Z a bound below 2 needs a hybrid Weyl-type
  input, see the remark there; we do not rely on this.) Theorem C beats 3/(3σ − 1) exactly on (5/6, 7/8] (E2).
- Principal characters (review gap 1): repaired in Lemma Z; the ζ part is covered by Ivić/Bourgain's 3/(2σ) (Thm C)
  and by DensityB's 3/τ0 (Thm D).
- Theorem C is unconditional; its only recent input is [PY] (2019/2023). The half-plane contributes nothing to the
  uniform exponents (Proposition B); its contributions are Theorem A (large values) and Theorem D (a sub-regime).
- Novelty searches (REPORT.md) found no use of [PY] in family zero-density estimates.
