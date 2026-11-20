RESULTS: 2 confirmed of 4 A-items (27 ids triaged; 2 A-items give no gain)

# GRHHarv5 — GRH-conditional arXiv papers of 1 Jul – 31 Dec 2025 (grhharv5, Orca, 20 Nov 2026, from 10:57 UTC)

Harvest: `.scratch/axq.py` with the brief's query (abs: GRH / generalized / generalised / extended Riemann hypothesis, submittedDate
[202507010000 TO 202512312359]) returned 40 ids (< 100, no split; `harvest.txt`). Skipped as SEEN: R504's ids (2509.00635, 2507.06913,
2510.14306, 2507.06609, 2508.17701, 2512.01779, 2512.00337, 2511.16452, 2509.21518, 2512.19761) and 2511.23381, 2509.24119, 2508.15715,
2508.14888, 2508.11125. Triaged: R504's four A-candidates, the 21 NEW ids, and two GRH papers of the period the query missed, found in the
references of 2509.20192 (2509.07651: abstract misspells "Riemnn"; 2508.16967: abstract says only "conditional"). Classes: A = GRH only
through a power-saving prime/Chebotarev/character-sum count or a least-prime bound; B = needs the critical line, GRH for L-functions outside
our half-planes, or a GRH-conditional machine computation; C = applies a half-plane itself (none); 0 = no GRH-conditional result of its own.
Inputs: L(s, χ) ≠ 0 for Re s > θ₁₄₁ = 0.87334585610671857680… (every Dirichlet χ), ζ_K(s) ≠ 0 for Re s > 31/32 (every K), both Lean-checked;
(H) = Frei–Sofos–Widmer 2610.18491 Thm generaldelta. Working notes: `NOTES.md`; reviewer replies: one 'reviewer' subagent per item.

| id | title (short) | class | GRH step (TeX line) |
|---|---|---|---|
| 2510.02309 | Effective Brauer–Siegel for Artin L-functions | A | Prop grh-bound l.262–271; proof l.1010–1048: only LO Chebotarev x^{1/2} log(D_K x), x ≥ (log D_K)² (l.1012–1019) |
| 2509.20192 | Large quadratic character sums, multiplicative coefficients | A | Thms 1.2 (l.246), 1.3 (l.254) via Darbar–Maiti Lemma 2.2 (l.265–271), used l.322–330, 395–416 |
| 2509.07651 | Large quadratic character sums (Dong–Zhang) | A | Thms 1.1–1.3 (l.242, 260, 268) via DM Lemma 2.2 (l.294), used l.363, 401 |
| 2512.24147 | Large values of quadratic character sums revisited | A | Thm 1.1 l.239–244 via DM Lemma 2.2 l.265–276, used l.316–329 |
| 2510.09005 | Note on large quadratic character sums | A | Thm 1.2 l.259–264 via DM Lemma 2.2 l.282–288, used l.336–352 |
| 2508.16967 | Quadratic character sums, multiplicative coefficients | A | Thm 1.1 l.294–298 via DM Lemma l.308–316, used l.435–439 |
| 2510.14167 | Maximal order of ω*(n) = #{p : p − 1 ∣ n} | A | π(x; d, 1) ≫ x/(φ(d) log x) for d ≤ x^{1/2}/(log x)³ (l.203–205, 248–252); principle l.372–383 |
| 2511.13659 | Average hardness of SIVP, module lattices | B | ERH for all trivial-modulus Hecke L (introduction.tex l.216, 438; self-reducibility-bulk.tex l.18, 579–592): Größencharaktere of Pic⁰_K (preliminaries.tex l.856); only unit rank 0 is covered, where n = 2r is fixed |
| 2512.01588 | Class/unit groups in provable subexponential time | B | ERH for characters of Arakelov ray class groups (introduction.tex l.41–42, 148; dBDPW random walks), infinite order |
| 2512.23302 | Chebyshev's bias without LI | B | zeros on the line in the explicit formula (l.100–122) |
| 2512.20066 | One-level density, Γ₁(q) | B | GRH for the family's zeros (l.183–185) |
| 2511.10594 | Binary correlations of cusp-form coefficients at primes | B | Littlewood zero-free region for L(s, f, χ), GRH for L(s, f) (l.127–131, 199–236) |
| 2511.08783 | Keating–Snaith for cubic Hecke L over Q(ω) | B | GRH for L(s, χ_𝔣), 1-level density over zeros 1/2 + iγ (l.183, 511–515) |
| 2511.02344 | Lower bounds, high moments of twisted coefficients | B | Cor l.199–205 adds the GRH upper bound (critical-line moments); lower bounds unconditional |
| 2510.12092 | x¹³ + y¹³ = zⁿ | B | GRH-conditional Magma class groups for J₇'s 2-Selmer rank (l.658–671) |
| 2509.20085 | High moments, quadratic twists of modular coefficients | B | Soundararajan–Harper moments under GRH (l.144–231) |
| 2509.01602 | Joint equidistribution of newforms | B | GRH for newform L, short Euler products at s = 1/2 (l.159–203, 1015) |
| 2508.15964 | Decorrelation of SO(5)×SO(2) Bessel periods | B | GRH for L(s, sym³f × χ_d) etc. (l.191, 556–570) |
| 2508.14534 | Shifted moments, cubic/quartic Dirichlet L | B | moments on the line under GRH (l.140–171) |
| 2508.13301 | Argument of Dirichlet L-functions | B | Selberg-type S(T, χ) on the line (l.94–107) |
| 2508.07734 | Selberg orthogonality, half-integral weight | B | GRH + GRC (l.143–146, 236) |
| 2507.14566 | Effective non-vanishing of Hecke–Maass L | B | GRH for L(s, f) (l.202, 532) |
| 2507.05862 | Norm-Euclidean cyclic cubic fields | 0 | Thm T:big is unconditional (l.267–275); GRH only in McGown's prior result (l.264–266, 372) |
| 2510.27573 | Sárközy's theorem in F₂[x] | 0 | Weil's RH for function fields (l.87, 106) |
| 2510.25630 | Average rank over function fields | 0 | Weil's RH (l.454); GRH only in Young's result over Q (l.160, 195) |
| 2508.08996 | Artin-type problems over number fields | 0 | unconditional under series conditions (l.84, 103, 177) |
| 2510.05501 | Class and unit groups in Magma (survey) | 0 | describes GRH-conditional regimes (l.112, 134–154); Bach's bound with (c log Δ)^{32} is Harvest4's 'E (minor)' |

## A1 — 2510.02309 (effective Brauer–Siegel for Artin L-functions, Prop grh-bound): CONFIRMED

- Source: K/k Galois, G = Gal(K/k), χ any character of G, κ(χ) = lim_{s→1}(s − 1)^{⟨χ,1_G⟩}L(s, χ) (l.95–97), q(χ) the conductor, χ̃ the
  induction of χ to Q through the Galois closure, μ(χ) = min_g Re χ(g) (l.99–102). Prop grh-bound (l.262–271): GRH for ζ_K ⇒
  |κ(χ)| ≪ (log log eD_k)^{χ̃(1)−χ(1)}(log log q(χ))^{χ(1)−⟨χ,1⟩} and |κ(χ)| ≫ (log log eD_k)^{μ(χ̃)−μ(χ)}(log log q(χ))^{μ(χ)−⟨χ,1⟩}, constants
  depending on [k:Q], |G|, χ(1). Unconditional Thm cond (l.131–145): log in place of log log, lower bound times ε(χ) (potentially exceptional
  quadratic characters, l.140–143). The authors: Prop grh-bound "indicates how Thm cond is near the limit of existing unconditional methods" (l.272–275).
- GRH uses (rg 'GRH|Riemann|zero-free': l.65, 260–275, 1010–1048; the rest is Stark/Landau–Siegel material of the unconditional theorems):
  only the proof l.1010–1048: (i) LO Chebotarev Σ_{N𝔭≤x} 1_C(𝔭) = (|C|/|G|)Li(x) + O(x^{1/2} log(D_K x)), x ≥ (log D_K)² (l.1012–1015), and
  by orthogonality Σ χ(𝔭) = ⟨χ,1⟩Li(x) + O(x^{1/2} log(D_K x)) (l.1016–1019); (ii) Mertens over k for x ≥ y ≥ log D_k (l.1020–1023; primes in
  (log D_k, (log D_k)²] discarded trivially); (iii) tail product with {1 + O(y^{−1/2} log(D_K y))}^{χ(1)} (l.1025–1029); (iv) short Euler
  product for T ≥ (log D_K)² (l.1029–1033); (v) prop:residue applied to faithful characters on subextensions with T(ψ) ≥ log q(ψ), T(1) =
  log eD_k, T(ψ) = log q(χ) (l.1033–1040); (vi) prop:small-primes with the GRH Mertens (l.1041–1047).
- Replacement: (H) with δ = 1/32 for K/k needs only ζ_K ≠ 0 on Re s > 31/32: |π_C(x) − (|C|/|G|)Li(x)| ≤ c₁(|C|/|G|)x^{31/32}(32 log D_K +
  n_K log x) for all x > 2. (i') Σ χ(𝔭) = ⟨χ,1⟩Li(x) + O(x^{31/32} log(D_K x)). (iii') tail error ∫_y^∞ |E(t)|t^{−2}dt + |E(y)|/y ≪
  y^{−1/32} log(D_K y), plus prime powers O(χ(1)[k:Q]/y): below 1/2 once y ≥ (C log D_K)^{32}, C = C([k:Q], |G|, χ(1)). (ii') the prime ideal
  theorem for k from (H) gives Mertens beyond (C log D_k)^{32}; on (log D_k, (C log D_k)^{32}] Σ 1/N𝔭 ≤ [k:Q] Σ_p 1/p + O(1) ≤ [k:Q] log 33
  + O(1) (as lem:mertens l.659–668 discards its range), so (ii) holds for x ≥ y ≥ log D_k. (iv') T ≥ (C log D_K)^{32}. (v') T(ψ) ≥
  (C log q(ψ))^{32} (log D_{K_ψ} ≪ log q(ψ) for the fixed field K_ψ of ker ψ, the paper's own reduction); the primes in (log q(ψ), (C log q(ψ))^{32}]
  are discarded: each Euler factor has modulus in [(1 + 1/N𝔭)^{−ψ(1)}, (1 − 1/N𝔭)^{−ψ(1)}] and their Σ 1/N𝔭 = O_{[k:Q]}(1), so that partial
  product is ≍ 1 in both directions and T(ψ) = log q(χ) stands. (vi) unchanged. No exceptional zero: real zeros β ≤ 31/32 < 1 − 1/(4 log D_K)
  for D_K > e^8 (smaller fields of bounded degree go into the constants). The half-plane holds for K, k and every subfield.
- R504's worry ('power saving only from x ≥ D_K^c') was wrong: (H) saves x^{−1/32} log D_K, a power saving from x ≥ (C log D_K)^{32} on.
- Conclusion (given the Lean-checked 31/32): Prop grh-bound holds unconditionally, with constants depending on [k:Q], |G|, χ(1).
- Reviewer (subagent, round 1): CORRECT, ROUTINE ("controls Euler-product tails; does not give relative Chebotarev asymptotics at that
  cutoff"; no further GRH use; the discard is valid for moduli both ways; small discriminants/conductors absorbed into constants).


## A2 — long quadratic character sums with the GCD resonator (2509.20192 Thm 1.3, 2509.07651 Thm 1.3, 2512.24147 Thm 1.1): CONFIRMED

- Sources (d over fundamental discriminants F, χ_d = (d/·)): 2509.20192 Thm 1.3 (l.254–258): GRH, exp((log X)^{1/2+ε}) < x ≤ X^{1/2}, every
  completely multiplicative f with |f| = 1, Re f(n)f̄(m) ≥ 0 (l.241–244): max_{X<|d|≤2X} |Σ_{n≤x} f(n)χ_d(n)| ≥ √x exp((1+o(1))√(log(√X/x)
  log_3(√X/x)/log_2(√X/x))); 2509.07651 Thm 1.3 (l.268–272): f = 1, X < d ≤ 2X; 2512.24147 Thm 1.1 (l.239–244): same range,
  max |Σ_{n≤|d|/x} χ_d(n)| ≥ √(X/x) exp(same). Method: resonator R_d = Σ_{m∈M} χ_d(m) over de la Bretèche–Tenenbaum's GCD-extremal set
  (|M| = N, Σ_{m,n∈M} √((m,n)/[m,n]) = N exp((2+o(1))√(log N log_3 N/log_2 N)), max P₊(m) ≤ (log N)^{1+o(1)}, max M ≤ 2 min M; elements
  up to exp((log N)^{1+o(1)})); 2512.24147 first applies Pólya's expansion (z = √(|d|x) log|d|, C_d(z) = Σ_{1≤|k|≤z} a_k χ_d(k), a_k =
  (1 − cos 2πk/x)/k, l.296–308).
- GRH use (rg 'GRH|Riemann|zero-free|quasi' each, every hit read): only Darbar–Maiti's Lemma (Math. Ann. 2025; Lemma 2.2: 2512.24147
  l.265–276, 2509.20192 l.265–271, 2509.07651 l.294–300): Σ_{|d|≤X, d∈F} χ_d(n) = (X/ζ(2))∏_{p|n} p/(p+1)·1_{n=□} + O(X^{1/2+ε}f(n₀)g(n₁)),
  n = n₀n₁², f(n₀) = exp((log n₀)^{1−ε}), g(n₁) = Σ_{d|n₁} μ²(d)d^{−1/2−ε}; used in M₁ = Σ_d R_d², M₂ = Σ_d R_d²|S_d|² (2512.24147 l.316–329;
  2509.20192 l.395–426; 2509.07651 l.401ff.), with N = X^{1/2−δ}/x (2512.24147 l.310) or X^{1/2−κ}/x (2509.20192 l.370, 424).
- Lemma DM_θ (θ = θ₁₄₁): for each δ > 0 there is η = η(δ) > 0 with error O_δ(X^{θ+δ} exp((log n₀)^{1−η}) G(n₁)), G(n₁) = ∏_{p|n₁}(1 + 3p^{−θ}).
  Proof: for n₀ > 1, d ↦ χ_d(n₀) is a non-principal character ψ mod 4n₀ or 8n₀; the series Σ_{d∈F, (d,n₁)=1} ψ(d)d^{−s} (both signs, the
  2-adic conditions through characters mod 8) is a finite combination of L(s, ψξ)/L(2s, ψ²ξ²) times local factors at p | 2n₁ (for n₀ = 2 the
  principal combinations' residues cancel: n is no square). Perron with T = X, shift to Re s = θ + δ: there |log L(s, χ)| ≪_δ
  (log q(|t| + 2))^{1−η} (Borel–Carathéodory and Hadamard's three circles, using only zero-freeness on Re s > θ and convexity; GRHHarv4 A3),
  1/L(2s) is bounded (2θ > 1), the local factors are ≤ G(n₁). n = □: unconditional squarefree count.
- Sizes: log n ≤ 2(log N)^{1+o(1)} + 2 log z = (log X)^{1+o(1)}, so exp((log n₀)^{1−η}) = X^{o(1)} and log G(n₁) ≤ 3Σ_{p|n₁} p^{−θ} ≪ (log 3n₁)^{1−θ}
  = o(log X) (sum over the first ω(n₁) primes); the paper's f(n₀) ≤ exp(P₊(n)^{1−ε}) plays the same role.
- Replacement, short sums (2509.20192, 2509.07651): M₁ ≤ XN/ζ(2) + O(X^{θ+δ+o(1)}N²) ≪ XN; M₂ error O(X^{θ+δ+o(1)}N²x²) (a, b ≤ x, |f| = 1);
  N = X^{1−θ−κ}/x, δ < κ/2: error/M₁ ≪ X^{δ−κ+o(1)}x, against the main term (1/N)Σ_{m,n∈M} #{a, b ≤ x : am = bn} ≫ x exp((2+o(1))√(log N log_3 N
  /log_2 N)) (2509.20192 l.426–463, unconditional). Pólya sums (2512.24147): (Σ_{0<|k|≤z} |a_k|)² ≪ (log z)², so the M₂ error is
  O(X^{θ+δ+o(1)}N²(log z)²); N = X^{1−θ−κ}/x: error/M₁ ≪ x^{−1}X^{δ−κ+o(1)}, against the main term x^{−1}(log X)^{−c} exp((2+o(1))√(…)) (l.329–370).
  The GCD-tail step (l.351–360: x^{−2η}exp(y_M^{1/2+η}) ≪ 1 for x > exp((log X)^{1/2+δ})) and the combinatorics are unconditional. Take δ, κ → 0
  slowly; x ≤ X^{1−θ−ε} gives log(X^{1−θ−κ}/x) = (1 + O(κ/ε)) log(X^{1−θ}/x).
- Erratum (2512.24147 l.322–329, 363–370): the written M₂ error X^{1/2+ε}N²x² (Σ_{k,ℓ≤x} 1) leaves the ratio X^{−δ+ε}x, which beats the main
  term x^{−1}exp(…) only for x ≤ X^{δ/2}; the (log z)² bound above repairs it, and 'x/N', 'x exp' at l.363–370 should read x^{−1}/N, x^{−1}exp.
- Conclusion (given θ₁₄₁): for fixed ε > 0 and exp((log X)^{1/2+ε}) < x ≤ X^{1−θ₁₄₁−ε}: (i) max_{X<|d|≤2X} |Σ_{n≤x} f(n)χ_d(n)| ≥ √x exp((1+o(1))
  √(log(X^{1−θ₁₄₁}/x) log_3(X^{1−θ₁₄₁}/x)/log_2(X^{1−θ₁₄₁}/x))) for every f as above (f = 1, d > 0 likewise); (ii) max_{X<|d|≤2X}
  |Σ_{n≤|d|/x} χ_d(n)| ≥ √(X/x) exp(same).
- Why the half-plane is needed: Granville–Soundararajan's unconditional Lemma 4.1 (error X^{1/2}n^{1/4} log n) needs n ≤ X^{O(1)}, i.e. a
  resonator on [1, X^{O(1)}]; such resonators give only exp(c√(log X/log_2 X)) (Soundararajan's r(p) = L/(√p log p), as in 2510.09005), the
  log_3 gain being the long-resonator phenomenon (Bondarenko–Seip; de la Bretèche–Tenenbaum 1804.01629, who for characters mod q reduce M
  modulo q, r(h)² = #{m ∈ M : m ≡ h (q)}, l.1057, a device with no analogue in the d-family). Not a proof that no polynomial-support route exists.
- Reviewers: round 1 GAP (slow δ, κ; '4^{ω(n)} = X^{o(1)}' unsupported; δ < κ; GS comparison open), all but the last repaired above;
  round 2: points 1, 2, 5 OK, point 3 'GAP' read the short-sum case with the Pólya main term x^{−1}exp (for short sums the main term is
  x exp, 2509.20192 l.426), point 4 GAP: the GS comparison is heuristic; round 3 (narrow, on point 3): CORRECT, ROUTINE (x cancels in
  both comparisons, relative error X^{−κ/2+o(1)}(log X)^c). Novelty relative to the literature found, with the comparison caveat stated.
- GRHHarv3 filed 2608.15773 (Aug 2026, same DM Lemma 2.2, l.546–551) as B ('Lindelöf on Re s = 1/2'); DM_θ shows the lemma transfers, so its
  long-range Omega results are A-candidates for a later round.


## A items that give no gain

- A3 — polynomial-resonator Omega results for quadratic character sums (2509.07651 Thms 1.1–1.2 l.242–264, 2508.16967 Thm 1.1 l.294–298,
  2509.20192 Thm 1.2 l.246–251, 2510.09005 Thm 1.2 l.259–264): NO GAIN. GRH enters only through the same DM Lemma, with resonators supported on
  n ≤ y = X^{1/2−δ}/(2 log z)², X^{1/2−α}/x² etc. (2510.09005 l.328, 2509.20192 l.307); DM_θ gives y ≤ X^{1−θ₁₄₁−δ}, so constants like √2/2 become
  √(1 − θ₁₄₁) = 0.35588501498838276206… (e.g. 2510.09005: √(X/x) exp((√(1−θ₁₄₁) + o(1))√(log X/log_2 X))). But the unconditional
  Granville–Soundararajan Lemma 4.1 (error X^{1/2}n^{1/4} log n; cited at 2509.07651 l.279–291, 2509.20192 l.264) already allows y ≤ X^{1/3−δ}/x
  for short sums (M₂ error ≈ X^{1/2}y^{3/2}x^{5/2}Σr², main ≈ X x Σr²) and y ≤ X^{1/6−δ}x^{−5/6} for the Pólya-expanded sums (Σ_{|k|≤z}|a_k|k^{1/4}
  ≪ z^{1/4}), i.e. constants √(1/3), √(1/6) > √(1 − θ₁₄₁); for small y the GS error is the smaller one (X^{1/2}n^{1/4} < X^{θ} for n < X^{4θ−2}).
  Reviewer: CORRECT, ROUTINE ("the displayed constants require x = X^{o(1)}"; a comparison of the supplied methods, not an impossibility result).
- A4 — 2510.14167 (ω*(n) > exp((log((1+√5)/2) + o(1)) log n/log log n) i.o. under GRH, l.174–183): NO GAIN. GRH is used only for π(x; d, 1) ≫
  x/(φ(d) log x) at d = x^{1/2−ε+o(1)}, d | k (l.203–205, 248–252); by the paper's own principle (l.372–383) a range d ≤ x^α gives the constant
  max_t f_α(t), f_α(t) = (t log t − (t−α) log(t−α) − α log α)/(t + 1 − α) (l.320), whose maximum is −log(1 − τ), 1 − τ = τ^α (stationarity of
  αH(τ)/(α + (1−α)τ), τ = α/t; prachar.py: α = 1/2 → 0.481211825059603 = log φ, α = 0.4736 → 0.466944945316581 = the paper's 0.4669).
  θ₁₄₁ gives π(x; d, 1) = Li(x)/φ(d) + O(x^{θ₁₄₁} log²(dx)), α = 1 − θ₁₄₁ − ε, constant 0.210472675714423 (0.3036 log 2), below the paper's
  unconditional 0.6736 log 2 (Harman's smooth moduli d < x^{0.4736}, l.296–312); density routes restricted to σ ≤ θ₁₄₁ keep the exponent
  12/5 near σ = 3/4 (range x^{5/12}). Reviewer: CORRECT, ROUTINE.
- R504's other two A-candidates: 2507.05862 is unconditional already (class 0); 2511.13659 needs ERH for infinite-order Hecke characters of
  the Arakelov class group (class B; R504's "these divide ζ_H" holds only for class-group characters, i.e. unit rank 0, where n = rd = 2r is
  bounded and SIVP is solvable directly).
- Lessons: (1) a GRH family lemma of Lindelöf type in the conductor (Darbar–Maiti: X^{1/2+ε}exp((log n₀)^{1−ε})) transfers to X^{θ+δ}exp((log n₀)^{1−η})
  via three circles (GRHHarv3's B for 2608.15773 is superseded); (2) it is worth something only where no unconditional lemma competes: for
  polynomial supports Granville–Soundararajan beats θ₁₄₁, for the long (GCD) resonator only the half-plane works; (3) a range parameter α
  entering an entropy optimisation (2510.14167) loses to any unconditional range α > 1 − θ₁₄₁.


## Ready COROLLARIES.md lines

- Effective Brauer–Siegel bounds with log log for Artin L-functions without GRH (N11; GRHHarv5, 20 Nov 2026; routine; reviewer CORRECT,
  ROUTINE): for every Galois extension K/k of number fields with group G and every character χ of G, κ(χ) = lim_{s→1}(s − 1)^{⟨χ,1_G⟩}L(s, χ)
  satisfies |κ(χ)| ≪ (log log eD_k)^{χ̃(1)−χ(1)}(log log q(χ))^{χ(1)−⟨χ,1_G⟩} and |κ(χ)| ≫ (log log eD_k)^{μ(χ̃)−μ(χ)}(log log q(χ))^{μ(χ)−⟨χ,1_G⟩}
  (μ(χ) = min_g Re χ(g), χ̃ the induction to Q; constants depending on [k:Q], |G|, χ(1)), given the Lean-checked 31/32 (Cureton–Thorner).
  Method: arXiv 2510.02309 (2 Oct 2025) proves this under GRH for ζ_K (Prop grh-bound, TeX l.262–271), using GRH only in its proof (l.1010–1048)
  through Lagarias–Odlyzko's Chebotarev error x^{1/2} log(D_K x) for x ≥ (log D_K)² (l.1012–1019) and the Mertens and Euler-tail estimates built on
  it; Frei–Sofos–Widmer 2610.18491 Thm generaldelta with δ = 1/32 gives the error x^{31/32}(32 log D_K + n_K log x), so the short Euler products hold
  from T ≥ (C log D_K)^{32}, the primes in (log q(ψ), (C log q(ψ))^{32}] change the products by factors ≍ 1 (Σ 1/N𝔭 ≤ [k:Q] log 33 + O(1)), and real
  zeros ≤ 31/32 leave no exceptional zero. Before: unconditionally log in place of log log and an exceptional factor ε(χ) in the lower bound (the
  paper's Thm cond, l.131–145; the authors expect that to be "near the limit of existing unconditional methods", l.272–275). Novelty: arXiv API
  (abs "Brauer-Siegel", Oct 2025 – Nov 2026): 2610.31402 (effective Siegel zeros over a base field), 2601.18408, 2608.12257, other questions;
  math.NT 2026-10/11 listings: nothing on it; Zenodo "Brauer-Siegel Artin": nothing. [GRHHarv5](Riemann/Deep/GRHHarv5/REPORT.md)
- Long quadratic character sums with the log₃ gain without GRH (N11; GRHHarv5, 20 Nov 2026; routine; reviewers: round 1 GAP (slow δ, κ;
  the G(n₁) bound), repaired; rounds 2–3 CORRECT, ROUTINE; the comparison with Granville–Soundararajan's unconditional route is heuristic): for
  fixed ε > 0, exp((log X)^{1/2+ε}) < x ≤ X^{1−θ₁₄₁−ε} and every completely multiplicative f with |f| = 1 and Re f(n)f̄(m) ≥ 0,
  max_{X<|d|≤2X} |Σ_{n≤x} f(n)χ_d(n)| ≥ √x exp((1+o(1))√(log(X^{1−θ₁₄₁}/x) log₃(X^{1−θ₁₄₁}/x)/log₂(X^{1−θ₁₄₁}/x))) and max_{X<|d|≤2X}
  |Σ_{n≤|d|/x} χ_d(n)| ≥ √(X/x) exp(same), d over fundamental discriminants, given θ₁₄₁ = 0.87334585610671857680…. Method: arXiv 2509.20192
  (Thm 1.3, l.254–258), 2509.07651 (Thm 1.3, l.268–272, f = 1) and 2512.24147 (Thm 1.1, l.239–244) prove these with √X/x and x ≤ X^{1/2}
  under GRH, using it only through Darbar–Maiti's lemma (Σ_{|d|≤X} χ_d(n) = main + O(X^{1/2+ε}exp((log n₀)^{1−ε})g(n₁)), 2512.24147 l.265–276)
  at de la Bretèche–Tenenbaum's GCD resonator (elements up to exp((log N)^{1+o(1)})); the half-plane gives |log L(s, ψ)| ≪ (log q|t|)^{1−η} on
  Re s ≥ θ₁₄₁ + δ (Borel–Carathéodory, three circles), hence the lemma with X^{θ₁₄₁+δ}, and the resonator length X^{1−θ₁₄₁−κ}/x. Erratum: the
  M₂ error of 2512.24147 (l.322–329) is O(X^{1/2+ε}N²(log z)²), not N²x² (as written its proof covers only x ≤ X^{δ/2}). Before: GRH only;
  Granville–Soundararajan's unconditional lemma (error X^{1/2}n^{1/4} log n) supports only polynomially bounded resonators, which give no
  log₃ gain (heuristic comparison, not an impossibility proof). Novelty: arXiv API (abs "quadratic character sum(s)", Sep 2025 – Nov 2026): only
  GRH papers (2509.07651, 2510.09005, 2512.24147, 2608.15773); math.NT 2026-10/11 listings: 2611.11389 (under RH), 2610.24958 (GRH);
  Zenodo "quadratic character sums": nothing. [GRHHarv5](Riemann/Deep/GRHHarv5/REPORT.md)

Files: `harvest.txt` (API output with abstracts), `dl.sh` (e-print download and hit extraction into .scratch/grhharv5), `hits.py`, `abs.py`
(hit and abstract filters), `prachar.py` (the 2510.14167 constants, mpmath), `NOTES.md` (working notes). TeX sources were in .scratch/grhharv5/
and are deleted; reviewer replies were in .scratch/grhharv5/review1.txt, review2.txt (quoted above).

Scope: removes GRH from an Artin/Brauer–Siegel statement and from three long-range quadratic character-sum Omega results with the
Lean-checked half-planes; no new zero-free input.

RH relevance: N11
