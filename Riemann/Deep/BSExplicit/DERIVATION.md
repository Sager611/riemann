# bsexpl: derivation (Orca agent bsexpl, 2026-11-06)

Goal: an explicit c(k) with (1 − β) log d_F ≥ c(k) for every quadratic F/k and every real zero β ∈ (0,1) of
L(s, χ_{F/k}) = ζ_F(s)/ζ_k(s), k ≠ Q; then an explicit effective Brauer–Siegel constant over k.
Notation: n = [k:Q]; d = |discriminant|; M = normal closure of F over Q; k̃ = normal closure of k over Q;
T* = B_2 + 2.2·10⁻²⁵⁴, B_2 = (777 − √1081)/852 = 0.8733819666… (zero-free half-plane Re s > T* for every Dirichlet
L-function, Lean-kernel-checked, Bebel 1 Nov 2026, TREE.md N7); fallback 7/8 (OpenAI, Lean port Riemann/Deep/QRH).

## 1. Banks's chain (arXiv 2610.31402v2, TeX `BrauerSiegel-CURRENT.tex`) and what is not explicit

c(k) of cor:siegel (l. 266, proof l. 785–822) = min( 1/(16 c_1 C), min_{small F} (1 − β_F) log d_F ), C = max(C_0(1/8), 8c_1, 8).

| # | input (TeX line) | status in the paper |
|---|---|---|
| 1 | lem:count (l. 1510): A, W, E_p with errors O(D^{s−1}), O(D^s) | implied constants not given |
| 2 | prop:main (l. 1553): "N large in terms of n_K, H", error O_{n_K,H}(1) | not explicit |
| 3 | thm:precise (l. 381): threshold X_0(n_K, H), error O_{n_K,H}(1) | not explicit |
| 4 | thm:quadratic (l. 244, proof l. 648–703): ε = 1/8, C_0 = 32, prime ℓ with ℓ − 1 ≥ 64n, H = ⌈128 n(ℓ − 1)⌉, n_K ≤ 2n(ℓ − 1); "d_F large in terms of k, ℓ, H, C"; O_{k,ℓ,H,C}(1); C_k of lem:generator | thresholds and O(1) not explicit; for k = Q(i): ℓ = 131, H = 33 280, n_K ≤ 520 |
| 5 | lem:bias (l. 725): absolute c_1 (Stark (8)–(10), "O([F:Q])", Minkowski) | not stated; easy |
| 6 | Mertens for k (l. 770), O_k(1) | explicit in Garcia–Lee (Ramanujan J. 2022) |
| 7 | small fields (l. 806–822): Hermite + mean value theorem + Rademacher's convexity bound | effective; gives (1 − β_F) ≳ d_F^{−A} at the threshold |
| 8 | lem:residue (l. 859) = Stark Lemma 4: c_2 = 1/c_4 | explicit: c_4 = 43.17 (§5 below) |
| 9 | η_k (l. 923), Stark Lemma 3 + Thm 3 (l. 837) | explicit / qualitative |

Items 1–4 and 7 make Banks's c(k) astronomically small even if every O(1) were computed: the threshold sits above
X_0(520, 33 280) and every field below it is handled by a mean-value bound of size d_F^{−O(1)}. We did not evaluate it:
§§2–4 give a far larger explicit c(k) by a classical route that bypasses items 1–7.

## 2. Lemma D (a double real zero is far from 1)

Lemma D. Let M ≠ Q and let β ∈ (0,1) be a zero of ζ_M of multiplicity ≥ 2. Then (1 − β) log d_M > 6 − 4√2 = 0.34314…

Proof. Stark 1974, p. 140 (from his (9), (10); all Γ-terms are negative for 1 < σ < 2): for 1 < σ < 2 and any
conjugation-closed sub-multiset of zeros, Σ' 1/(σ − ρ) < 1/(σ − 1) + ½ log d_M. Keep the double zero:
2/(σ − β) < 1/(σ − 1) + L/2, L = log d_M ≥ log 3. Put b = (1 − β)L and suppose b ≤ 6 − 4√2. Take σ = 1 + a/L with
a = (1 + √2)b ≤ 0.829 < log 3 (so σ < 2). Then 2/(σ − β) − 1/(σ − 1) = L(2/(a + b) − 1/a) = L(3 − 2√2)/b ≥ L/2,
contradiction. ∎ (Stark's Lemma 3 gives the constant 1/4 for "at most one zero"; the double real zero gives 6 − 4√2.)

## 3. Theorem A (explicit c(k))

Inputs (Stark, Invent. Math. 23 (1974), read from the GDZ scan; quotes verbatim):
- Lemma 3 (p. 139): "if M ≠ Q, ζ_M(s) has at most one zero in the region S given by all s with σ ≥ 1 − (4 log|D_M|)^{−1},
  |t| ≤ (4 log|D_M|)^{−1}; if such a zero exists, it is real and simple."
- p. 140: "For s = σ with 1 < σ < 2, all the terms on the right of (9) after ½ log|D_M| are negative and thus by (9) and (10),
  for 1 < σ < 2, Σ′_ρ 1/(σ − ρ) < 1/(σ − 1) + ½ log|D_M|" (Σ′: any subset of zeros closed under ρ ↦ ρ̄).
- Theorem 3 (p. 144): "Let B be a normal extension of A and let χ be a character of G(B/A). Suppose ρ is a simple zero of
  ζ_B(s). Then L(s, χ, B/A) is analytic at s = ρ. Further, there is a field F, cyclic over A and contained in B such that for
  any field E between A and B, ζ_E(ρ) = 0 if and only if F ⊂ E. If ρ is real then F is either A itself or quadratic over A."
- Lemma 4 (p. 140): "Let M be a number field other than Q and let κ be the residue of ζ_M(s) at s = 1. We let β_0 be the
  exceptional zero (if it exists) of ζ_M(s) so that β_0 ≥ 1 − (4 log|D_M|)^{−1}. There is a constant c_4 > 0 (effectively
  computable) such that κ > c_4^{−1}(1 − β_1) where β_1 = β_0 if β_0 exists and β_1 = 1 − (4 log|D_M|)^{−1} otherwise."
(S3) Stark 1974 Thm 3, p. 144: B/A normal, ρ a simple zero of ζ_B ⇒ there is F′ cyclic over A, A ⊆ F′ ⊆ B,
with ζ_E(ρ) = 0 ⟺ F′ ⊆ E for every A ⊆ E ⊆ B; ρ real ⇒ F′ = A or [F′:A] = 2. (AB) Aramata–Brauer: E ⊆ B with B/E
Galois ⇒ ζ_B/ζ_E entire. (To) Toyama 1955: d_{K_1K_2} | d_{K_1}^{[K_1K_2:K_1]} d_{K_2}^{[K_1K_2:K_2]} (proof: for θ ∈ O_{K_2} with K_2 = Q(θ), minimal
polynomials g | f over K_1 and Q, 𝔇_{K_1K_2/K_1} | g′(θ) | f′(θ); the f′(θ) generate 𝔇_{K_2/Q}; take norms in 𝔇_{K/Q} = 𝔇_{K/K_1}𝔇_{K_1/Q}).

Theorem A. Let [k:Q] = n ≥ 2, F/k quadratic, β ∈ (0,1) with L(β, χ_{F/k}) = 0. Then one of:
 (a) L(β, χ_D) = 0 for the quadratic character of a quadratic field Q(√D) ⊆ F; then β ≤ T* and
     (1 − β) log d_F ≥ (1 − T*) log d_F ≥ (1 − T*) log 9 = 0.27820…; with a k = Q constant c_Q instead
     ((1 − β) log|D| ≥ c_Q): (1 − β) log d_F ≥ n c_Q, since |D|^n ≤ d_F;
 (b) (1 − β) log d_M > 6 − 4√2, and log d_M ≤ λ log d_F with λ = r[M:F] ≤ 2^{n−1} n!.
Hence (1 − β) log d_F ≥ c(k) := min{0.2782, (6 − 4√2)/λ_n} = (6 − 4√2)/(2^{n−1} n!) (half-plane version), and
≥ min{n c_Q, (6 − 4√2)/(2^{n−1} n!)} (c_Q version).

Proof. ζ_F = ζ_k L(s, χ_{F/k}) vanishes at β; by (AB) so does ζ_M. If β is a multiple zero of ζ_M, Lemma D gives (b).
If it is simple, (S3) with A = Q, B = M gives F′; E = F shows F′ ⊆ F; F′ ≠ Q because ζ(β) ≠ 0 on (0,1); so F′ = Q(√D)
and ζ_{F′}(β) = ζ(β)L(β, χ_D) = 0: case (a). Bounds in (b): M is the compositum of the r distinct conjugate fields of F;
σ(F) = σ(k)(√σ(α)) depends only on σ|_k, so r ≤ n; iterating (To), log d_M ≤ r[M:F] log d_F; M = k̃(√α_1, …, √α_n), so
[M:Q] ≤ 2^n n! and λ = r[M:Q]/(2n) ≤ 2^{n−1} n!. In (a): d_F = d_{F′}^{n} N(𝔡_{F/F′}) and d_F ≥ d_k² ≥ 9. ∎

Refinement, n = 2 (every quadratic k, in particular Q(i), Q(√−3)); this case needs neither (S3), (AB) nor (To).
F/Q is V4, C4 (only for real k: the quadratic subfield of a cyclic quartic field is real) or D4.
- V4 / C4: F/Q abelian, Ind_k^Q χ_{F/k} = ψ_1 ⊕ ψ_2 with Dirichlet characters ψ_i, so L(s, χ_{F/k}) = L(s, ψ_1)L(s, ψ_2) and
  β ≤ T*: (1 − β) log d_F ≥ (1 − T*) log 9 = 0.2782. (c_Q version: ψ_i quadratic ⇒ |cond ψ_i|² ≤ d_F ⇒ ≥ 2c_Q;
  ψ_i quartic ⇒ L(β, ψ̄_i) = 0 too, a double zero of ζ_F ⇒ Lemma D with M = F ⇒ ≥ 0.3431.)
- D4: M = F F̄ is V4 over k with quadratic subfields F, F̄ = τF, F_3 (τ the non-trivial automorphism of k);
  L(s, χ_{F̄/k}) = L(s, χ_{F/k}) identically (reindex 𝔭 ↦ τ𝔭), so β is a zero of ζ_M = ζ_k L(χ_F) L(χ_F̄) L(χ_{F_3}) of
  multiplicity ≥ 2 (all factors entire near β): case (b) without (S3) or (AB). Conductor–discriminant over k and
  𝔣(χ_{F_3}) | 𝔣(χ_F)𝔣(χ_F̄) give d_M = d_k⁴ N𝔣(χ_F)² N𝔣(χ_{F_3}) ≤ d_k⁴ (d_F/d_k²)⁴ = (d_F/d_k)⁴, so λ ≤ 4.
So c(k) = (6 − 4√2)/4 = (3 − 2√2)/2 = 0.085786… for every quadratic k (c_Q version: min{2c_Q, 0.0857}, e.g. 1/50 with
Guo's c_Q = 1/100). The D4 case uses no k = Q input at all; the V4 case is pure k = Q input.

## 4. Theorem B (explicit effective Brauer–Siegel over k)

Theorem B. Let [k:Q] = n ≥ 2 and K/k Galois. (i) Every real zero β ∈ (0,1) of ζ_K has (1 − β) log d_K ≥ c_1(k) :=
min{6 − 4√2, c(k), η_k}, η_k ≥ min{(1 − T*) log 3, (6 − 4√2)/n!} (η_k = +∞ if ζ_k has no real zero in (0,1), e.g.
k = Q(i), Q(√−3)). (ii) h_K R_K > (2 min{1/4, c_1(k)}/c_4) · d_K^{1/2}/((2π)^{[K:Q]/2} log d_K), c_4 = 43.1621…
For n = 2: c_1 = (3 − 2√2)/2 and the constant is 2(3 − 2√2)/(2c_4) = 3.9750·10⁻³ (Banks's normalization
d_K^{1/2}/((2π)^{[K:Q]} log d_K) is weaker, so the same constant works there).

Proof. (i) Multiple zero: Lemma D with M = K. Simple: (S3) with A = k, B = K gives F′, k ⊆ F′ ⊆ K, [F′:k] ≤ 2. F′ = k:
ζ_k(β) = 0 and (1 − β) log d_K ≥ (1 − β) log d_k ≥ η_k (d_k ≤ d_K). [F′:k] = 2: ζ_k(β) ≠ 0 (F′ ⊄ k), so L(β, χ_{F′/k}) = 0,
and (1 − β) log d_K ≥ (1 − β) log d_{F′} ≥ c(k). η_k: a real zero β′ of ζ_k is a zero of ζ_{k̃} (AB); multiple ⇒ Lemma D and
log d_{k̃} ≤ r_k[k̃:k] log d_k ≤ n! log d_k (To); simple ⇒ (S3) with A = Q, B = k̃ gives a quadratic F″ ⊆ k with
L(β′, χ_{F″}) = 0, so β′ ≤ T* and (1 − β′) log d_k ≥ (1 − T*) log 3.
(ii) Stark Lemma 4 (§5): κ_K > (1 − β_1)/c_4 with β_1 the zero of ζ_K in S = {σ ≥ 1 − 1/(4 log d_K), |t| ≤ 1/(4 log d_K)}
(real, simple: Stark Lemma 3) or 1 − 1/(4 log d_K) if none. By (i), 1 − β_1 ≥ min{1/4, c_1}/log d_K. Class number formula
h_K R_K = w_K d_K^{1/2} κ_K/(2^{r_1}(2π)^{r_2}), w_K ≥ 2, 2^{r_1}(2π)^{r_2} ≤ (2π)^{[K:Q]/2} (as 2 ≤ √(2π)). ∎

## 5. Stark's Lemma 4 is explicit: c_4 = 2 exp(21/8 + (γ + 2 log 2)/(4 log 3)) = 43.1621…

Stark 1974, Lemma 2 (p. 138): f entire of order 1 as in his Lemma 1, (i) f′/f(σ) ≤ 1/σ + 1/(σ − 1) + c_1 log α on (1, 5/4],
(ii) d ≥ 3, d ≥ α, log d ≥ n/c_2, (iii) at most one zero in S, real ⇒ f(1) ≥ c_3^{−1} E f(σ_0), σ_0 = 1 + 1/(4 log d),
E = (1 − β_0)/(σ_0 − β_0) or 1, c_3 = exp(21/8 + c_1/2 − (c_2/8)ψ(1/2)) ((7), p. 139). Lemma 4 (pp. 140–141): f = s(s − 1)ζ_M,
α = d_M/(2^{2r_2}π^n), d = d_M, c_1 = 0 (ζ′/ζ < 0 on σ > 1), c_2 = 2/log 3, then κ = f(1) > (1 − β_1)/(2c_3).
Checked here: the arithmetic of (7) from the last display of p. 138 (2(σ_0 − 1)/σ_0 < 1/2, 2, (2c_1 + ½)/4, (c_2/8)|ψ(½)|);
σ_0 ≤ 1 + 1/(4 log 3) < 5/4; the h(x) ≤ 2h(x_1) case split for zeros outside S; E(σ_0 − 1) ≥ ½(1 − β_0);
c_2 = 2/log 3 ⟺ d_M ≥ 3^{n/2}: n = 2 trivially, n ≥ 3 from Minkowski and n! ≤ e√n (n/e)^n, since
(πe²/(4√3))^n = 3.351^n ≥ e² n for n ≥ 3 (cert.py §3).

## 6. Without the Lean-checked half-plane

Replace (a) by any explicit k = Q constant: Guo 2610.30476 c_Q = 1/100 (audited, Riemann/Deep/SiegelGuo), 2610.21511
c_Q = 5.811·10⁻⁴ (Lean threshold), OpenAI 7/8 half-plane (Lean port, Riemann/Deep/QRH): (1 − 7/8) log 9 = 0.2747.
The D4 part of the n = 2 result and Lemma D use nothing from 2026.

## 7. Comparison with the general-field half-planes (not used above)

Cureton–Thorner (Hexagon 2610.00200v2, Lean: AxiomMath/DedekindZetaZeroFree; our COROLLARIES treats it as Lean-checked):
ζ_F ≠ 0 on Re s > 31/32 gives (1 − β) log d_F ≥ (1/32) log d_F at once: 0.1553 for k = Q(i) (min d_F = 144),
0.1488 for k = Q(√−3) (min d_F = 117), and h_K R_K > (2 min{1/4, (log 3)/32}/c_4) d_K^{1/2}/((2π)^{[K:Q]/2} log d_K)
= 1.59·10⁻³ (…) for every number field K ≠ Q. Liu (arXiv 2610.26835, audited sound): Re s > 0.874887 for every
finite-order Hecke L-function: (1 − β) log d_F ≥ 0.1251 log d_F. Banks (l. 143–152) and Frei–Kai (2610.00466, l. 223)
both note that these preprints would settle the question.
