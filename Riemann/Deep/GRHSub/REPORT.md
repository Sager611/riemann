RESULTS: 3 of 4 confirmed
(Items 1 and 2 in full, item 4 for the full Sato–Tate range I = [0, π] only; item 3's claimed transfer fails, a positive-proportion version survives.)

# grhsub: four GRH-conditional theorems of Sep–Nov 2026 with a zero-free half-plane in place of GRH (R493 candidates, node N11)

Agent grhsub, 20 Nov 2026. Inputs: (Q1) ζ_K(s) ≠ 0 for Re s > 31/32, every number field K (Cureton–Thorner, Lean-checked; no human
referee, see HPAudit); (Q2) L(s, χ) ≠ 0 for Re s > θ = θ₁₄₁ = 0.8733458561067185768, every Dirichlet χ (QuasiRiemannTracker PR 141).
Chebotarev tool (H), citable: Frei–Sofos–Widmer, arXiv 2610.18491, Thm generaldelta (TeX l.194–207): if ζ_L has no zeros with
Re s > 1 − δ, δ ∈ (0, 1/2], then |π_C(x, L/K) − (#C/#G) Li(x)| ≤ c₁ (#C/#G) x^{1−δ} log(Δ_L^{1/δ} x^{n_L}) for every Galois L/K,
class C, x > 2 (c₁ absolute; LO77's GRH proof with |x^ρ| ≤ x^{1−δ}, |ρ| ≥ δ); Cor. quantitative_splitting_qrh: a prime of every class
with N𝔭 ≤ c(δ)(log Δ_L)^{1/δ}(log log Δ_L)^{2/δ}. With δ = 1/32 both follow from Q1 for ζ_L alone. Working notes: NOTES.md.
Reviewers: one 'reviewer' subagent per item, statement, GRH lemma and replaced steps pasted inline; verdicts quoted below.

## Ready COROLLARIES.md lines

- Generalized Euclidean pure orders of odd degree without GRH (N11; R493, 20 Nov 2026; routine; reviewer CORRECT, ROUTINE): for every fixed d ≥ 2 only finitely many pure orders Z[m^{1/d}] are generalized Euclidean (Johnson–Queen–Sevilla), given the Lean-checked 31/32. Method: arXiv 2610.25793 (19 Oct 2026, 'Orders in the Euclidean Court') proves this under GRH (Thm maingeneralGRH, TeX main-revised.tex l.91; even d unconditional) using GRH only through Prop. chebotarev (l.290–295, Serre's (20_R) over Q) in Lemma generalGRH (l.197–200, proof l.297–313), to get π_C(|m|^{1/6}) ≥ 2 + ω(p|m|) for E = Q(m^{1/p}, ζ_p), p the least prime of d, C the fixed-point-free elements, |C|/|G| = 1/p; (H) with δ = 1/32 gives |π_C(x) − Li(x)/p| ≤ (32c₁/p) x^{31/32} log(Δ_E x^{p(p−1)}), and log Δ_E ≤ (p(p−1) − 1) log(p|m|) + p(p−1) log(p(p−1)) (ramification only at pm, Serre Prop. 6), so π_C(|m|^{1/6}) ≫_d |m|^{1/6}/log|m| (error/main ≪ |m|^{−1/192}(log|m|)²); the rest of the proof is unconditional. Before: odd d only under GRH (the same paper). Novelty: lst.py (math.NT 2026-10/11 titles: Euclidean, pure order, quasi, Chebotarev) finds only the source and norm-Euclidean or class-number papers (2610.24900 removes GRH from Weinberger's UFD ⇔ Euclidean for totally real abelian fields via the quasi-RH: a different notion); Zenodo none.
- Koblitz–Zywina least prime divisors without GRH (N11; R493, 20 Nov 2026; routine; reviewer CORRECT after a bookkeeping fix, ROUTINE): every E/Q without congruence obstruction (CM or not) has infinitely many good primes p ≥ 5 with c_E(p) > κ log p, for every κ < 1/128 (c_E(p) the least prime divisor of #E_p(F_p)), given the Lean-checked 31/32. Method: arXiv 2609.30573 (24 Sep 2026) proves κ < 1/8 under GRH (Thm Quanticep, main.tex l.262–266) using GRH only through Thm ECDTA (l.1376: Bach–Sorenson's least prime p ∤ m, unramified in K, Frobenius in C, p ≪ (log d_K̃ + [K̃:Q])², K = Q(E[M]), M = Π_{ℓ≤n} ℓ, K̃ = K(√m), m = rad(2N_E M), C = {A: det(I − A) ∈ (Z/M)^×}); (H)'s corollary with δ = 1/32 (or chebotarev-numeric (A)) gives p ≤ c(log d_K̃)^{32}(log log d_K̃)^{64} unramified in K̃, hence p ∤ m (m ≡ 2 mod 4, every prime of m ramifies in Q(√m)), with Frobenius in the class of a lift of some A ∈ C; log d_K̃ ≪_E M^4 log M (l.1417–1436) gives log p ≤ 128 log M + 96 log log M + O_E(1) = (128 + o(1))n while c_E(p) > n. Before: Jones (MR2727656, appendix) unboundedness only; unconditional least primes ≪ Δ^{16} (Kadiri–Ng–Wong) give c_E(p) > (1/4 − ε) log log p; GRH κ < 1/8; Liu's unrefereed B₂ (chebotarev-numeric (B)) would give κ < 1/31.59. Novelty: lst.py (Koblitz, least prime divisor, congruence obstruction, elliptic curve) finds only the source, 2611.01204 (entanglement) and 2610.20010 (CM sieve for y² = x³ − x, which has an obstruction); Zenodo none. CM curves: CM sieves (Iwaniec–Jiménez Urroz) may give more; not checked.
- Divisor moments of F(a_f(p)) without GRH, full Sato–Tate range (N11; R493, 20 Nov 2026; routine; reviewer CORRECT, ROUTINE): for f a non-CM cuspidal eigenform of weight k ≥ 2, level N, trivial nebentypus and trivial inner twist, F = c Π F_i^{e_i} ∈ O_{K_f}[T] and m ≥ 1, Σ_{p ≤ x, p ∤ N, F(a_f(p)) ≠ 0} d_{K_f}(F(a_f(p)))^m ≍_{F,f,m} π(x)(log x)^{K_m}, K_m = Σ_v κ_v((v+1)^m − v^m), given the Lean-checked 31/32. Method: arXiv 2611.10489 (7 Nov 2026) proves this under GRH for Dedekind zetas (Thm application 2, Manuscript.tex l.1106–1120, case I = [0, π]) using GRH only in (a) LO77 Thm 1.1 (l.2616: error x^{1/2}(N𝔬)^4 log(xN N𝔬)) and (b) Serre's §7.2 bound #{p ≤ x: F(a_f(p)) = 0} ≪ x^{7/8}(log x)^{1/2} (l.2634); (H) with δ = 1/32 gives (a) with x^{31/32}(N𝔬)^4 log(xN N𝔬), and Serre's argument with (H) (ℓ split in K_f with image ⊇ SL₂(F_ℓ), trace class of density ≤ 2/ℓ, degree ≤ ℓ^4, ℓ ≍ x^{1/128}(log x)^{−1/2}) gives (b) with x^{127/128}(log x)^{−1/2}; then h(x,𝔬)/h(𝔬) = 1 + O(x^{−1/128}(log x)²(N𝔬)^{8+ε}), and Thm main_geometric (l.926–955, needs only some c ∈ (0, 1)) runs with c < 1/(128(8 + ε)) (the paper: c < 1/(8(10 + ε))). Not for I ≠ [0, π]: the joint Chebotarev–Sato–Tate input (l.2643–2655, LL26 Thm 1.3, Thorner 2021 Thm 1.3) needs zero-free regions for L(s, Sym^m f ⊗ χ) and Hypothesis A, which our half-planes do not give (unconditional Sato–Tate saves only powers of log x, the sieve-dimension estimate needs x^{−δ}). Before: GRH (the same paper; Lee26 for K_f = Q and squarefree F). Novelty: lst.py (divisor moment, vanishing density, Sato–Tate, Chebotarev, modular form) finds only the source and conditional Chebotarev–Sato–Tate papers (2610.36464, 2610.42022); Zenodo none; our Harvest6 (l.56) had filed 2611.10489 as needing hypotheses we lack, true only for I ≠ [0, π].
- Signs of the Alladi–Johnson constants without GRH, three quarters of the classes (N11; R493, 20 Nov 2026; routine, off-tree; reviewers: Markov step CORRECT, density refinement CORRECT, ROUTINE): for every ε > 0 and k ≥ k₀(ε), C(k,l) < 0 for at least (2(2 − θ)/3 − ε)φ(k) reduced classes l mod k, θ = θ₁₄₁, 2(2 − θ)/3 = 0.7511027625955209488 (C(k,l) is the constant in Σ_{n≤X, P⁻(n)≡l (k)} μ(n)ω(n)/n = C(k,l)/log X + O((log log X)⁴/(log X)²)). Method: arXiv 2610.23660 (16 Oct 2026) proves #{l: C(k,l) ≥ 0} ≪ φ(k)/log k under GRH (Thm grh, l.400–510) using GRH only through θ(x,χ) ≪ x^{1/2}(log kx)² (l.379–398) at X = k(log k)³, a Selberg pair sieve covering (k, X]; θ₁₄₁ does not give that conclusion by this route (the tail needs X ≥ k^{3.948} pointwise, and for X = k^A the uncentered pair bound only gives Σ_l F(l)² ≥ (Σ_l F(l))²/φ(k) ≍ (log k)²/φ(k)). Instead: Lemma lem:lem at X = k^A gives C(k,l) = 1_{l prime} w_l + F(l) − ζ(1,X)/φ(k) + T_X(k,l), w_p = ζ(1,p)/p, F(l) = Σ_{k<p≤X, p≡l} w_p ≥ 0; Parseval and the mean value Σ_{χ≠χ₀}|Σ_{p>y} χ(p)/p|² ≪ y^{−η} for y ≥ k^A, A > 3/(2(2 − θ)) (explicit formula, Montgomery's density estimate Σ_{χ mod k} N(σ,T,χ) ≪ (kT)^{3(1−σ)/(2−σ)}(log kT)^B, θ₁₄₁ only to cap β ≤ θ) make the classes with |T_X(k,l)| ≥ ε log k/φ(k) o(φ(k)); Markov with Σ_l F(l) = ζ(1,X) − ζ(1,k) and Mertens bounds the remaining classes with C(k,l) ≥ 0 by φ(k)(1 − (1 − b)/(A − b)), b = 2εe^{−γ}. The pointwise θ₁₄₁ bound alone gives 2(1 − θ) = 0.2533; Montgomery's large-σ estimate (2(1−σ)/σ, σ ≥ 4/5) would give 4/5 and Huxley's 12/5 would give 5/6 (not checked); the density hypothesis on [1/2, θ] would give the paper's 1 − o(1). Before: GRH (the same paper); numerics C(k,l) ≠ 0 for 3 ≤ k ≤ 300; nothing unconditional. Novelty: lst.py (Alladi, smallest prime, Möbius; math.NT 2026-10/11) finds only the source; Zenodo none.

## What fails (with reasons)

- Item 3, the R493 claim 'any θ < 1 gives the same φ(k)(1 − o(1))': false for this proof. Its two steps need incompatible X: the tail
  (eq:GRH) needs square-root cancellation at X = k^{1+o(1)} (θ₁₄₁ pointwise: X ≥ k^{1/(2(1−θ))+o(1)} = k^{3.948}; with zero density:
  X ≥ k^{1.331}), the middle range (eq:BT) needs X/k ≤ (log k)^{O(1)}. Only a positive proportion survives (line above).
- Item 4 for proper intervals I ≠ [0, π]: needs zero-free regions for L(s, Sym^m f ⊗ χ) (Chebotarev–Sato–Tate with power saving) and
  Hypothesis A of LL26; our half-planes cover Dedekind, Hecke and Dirichlet L-functions only.
- Item 2: nothing fails, but the exponent drops from 1/8 (GRH, Bach–Sorenson (log d)²) to 1/128 ((log d)^{32}); this is the method's limit
  for any least-prime input of shape (log d)^{1/δ}.

## GRH uses found (rg -n -i 'GRH|Riemann|zero-free|quasi' in each TeX, each hit read; details NOTES.md)

- 2610.25793: l.91, 182, 197–198, 289–295, 297 (one use: Prop. chebotarev). 2609.30573: l.262, 268, 1376 (one use: Thm ECDTA).
- 2610.23660: l.379–398 (one use: θ(x,χ) bound), l.400–510 (Thm grh); l.294 is the classical zero-free region (unconditional Thm thm).
- 2611.10489: l.1111 (hypotheses), l.2616 (LO77), l.2634 (Serre §7.2), l.2643–2655 (Sato–Tate, I ≠ [0, π] only); l.249, 302 are comments.

## Reviewer verdicts (one 'reviewer' subagent per item, inputs inline)

1: CORRECT, ROUTINE; no missed GRH use. 2: 'GAP as written' (64 log log log d is not O(1): fixed to 128 log M + 96 log log M + O_E(1);
the constant 0.484 needs chebotarev-numeric's own proof, (H)'s corollary suffices), conclusion CORRECT, ROUTINE. 3: obstruction (A) correct
after rephrasing (uncentered second moment), Markov version CORRECT ('not a pure substitution'), coefficient (1 − b)/(A − b); density
refinement (second reviewer): CORRECT given Q2, ZD citation correct to its knowledge (theorem numbers unchecked), ROUTINE. 4: CORRECT,
ROUTINE for I = [0, π]; proper intervals not removed.

## Lead (sketch, reviewer PLAUSIBLE; not a claim): item 2's theorem is weak; the sieve route transfers too

The paper notes (l.268) that sieve methods may beat Quanticep. With (H), Chebotarev in Q(E[d]) (degree ≤ d^4, class of density
δ_d = d^{−1+o(1)}) has error ≪ x^{31/32} d^3 log(dNx), so Σ_{d≤D} |r_d| ≪ x^{31/32} D^4 log x: level of distribution x^{1/128} for the
sequence #E_p(F_p), p ≤ x (GRH: x^{1/8}). The linear sieve (f(s) > 0 for s > 2; sieve dimension 1, no congruence obstruction, small
entangled primes fixed by one Chebotarev class at level A_E) should give ≫ x/(log x)² primes p ≤ x with no prime factor of #E_p below
x^{1/256 − ε}, i.e. c_E(p) > p^{1/256 − ε} infinitely often (vs κ log p above), hence #E_p = P_r with r ≤ 256 (≈ 128 needs a weighted
sieve with square-divisor control), for non-CM E, given 31/32 (GRH analogues: Miri–Murty 2001, Steuding–Weng, David–Wu 2012). Reviewer:
PLAUSIBLE, 'the load-bearing advance is the half-plane, not the sieve'; not known unconditionally to it. Worth one transfer round on the
GRH almost-prime Koblitz papers; not checked here.

RH relevance: N11
