RESULTS: 3 confirmed of 8 A-items (21 ids triaged; 3 A-items fail, 2 not pursued)

# GRHHarv3 — GRH-conditional arXiv papers of 1 Jul – 15 Sep 2026 (grhharv3, Orca, 20 Nov 2026, from 00:47 UTC)

Harvest: `.scratch/axq.py` with the brief's query (abs: GRH / generalized / generalised / extended Riemann hypothesis, submittedDate
[202607010000 TO 202609150000]) returned 25 ids (< 100, no split). Skipped as already in our files: 2609.14471 (R126), 2607.14515 (R330),
2609.27808 and 2608.16034 (N8Fam REPORT). Ids listed only in Harvest2/HARVEST.md (as "conditional", no substitution analysis) were triaged.
Classes: A = GRH only through a power-saving prime/Chebotarev count, least prime, character-sum or short-Euler-product bound; B = needs the
critical line; C = applies a half-plane itself; 0 = no GRH-conditional result of its own (GRH only cited or a citation key).

| id | title (short) | class | GRH step (TeX line) |
|---|---|---|---|
| 2610.00109 | PageRank on LPS graphs | A | Thm main-results (ii) l.216–222 via Lemma grh-short-euler l.944–986 (Granville–Soundararajan 2003 Lemma 2.1, σ₀ = 1/2, z = (log Y)^8); proof l.1380–1411 |
| 2607.21259 | Fractional parts of polynomials mod p | A | Thm estimate_GRH l.136–140 via Lemma l.204–209 = Montgomery–Vaughan 1977 Thm 2 (Σ_{n≤x} χ(n) ≪ √p log log p); proof l.399–404 |
| 2609.02392 | Large cubic/quartic L(1/2, χ), prime moduli | A | Thm l.128–136; Lemma PropDirpoly l.339–346 (prime sums of cubic/quartic residue symbols, error X^{1/2+ε} log(c+2)) via L'/L ≪ log for Re s ≥ 1/2+ε (l.413–417, IK Thm 5.17) |
| 2607.20853 | Murmurations of quadratic Hecke L over Q(i) | A | Thm quadraticmainthm l.110–116 (δ ∈ (3/4, 1)); Hecke prime sums ≪ √x log³x log N(n) (l.234–239) and the GRH prime ideal theorem in [yX, yX + X^δ] (l.414–418) |
| 2607.01126 | Selmer ranks in prime cyclic extensions | A | ERH only through effective Chebotarev (Thm effective_chebotarev l.485–490), used at l.1146–1150 and Prop Yu+GRH l.1663–1700 (degree ≤ B·2^{4g²ω(𝔡)}, ω(𝔡) ≤ 2 log log X) |
| 2609.16048 | Waring–Goldbach in short intervals (GRH) | A | ψ(x, χ) ≪ x^{1/2} log² x (l.329–331) and a θ = 1/2 + ε exponential-sum transform (l.339–341) |
| 2608.23500 | Logarithmic Chowla correlations | A | GRH theorem (abstract l.28–33; §16) via the GRH clause of KMT Thm 1.5 (§16 l.106–110) and Lemma grh-twisted-liouville (§16 l.114–125, X^{1/2+2η}) |
| 2608.02494 | Isogeny graphs in characteristic 0 | A | Thm l.5594ff via Banwait–Najman–Padurariu's GRH determination of IsogDeg(Q(√d)) (l.5578–5592; Larson–Vaintrob / Type 2 least-prime bounds) |
| 2607.01703 | Euclidean systems of ray classes | B | Lenstra's theorem = Artin's primitive root conjecture with density for units (l.60, 652–700, Cor 53 l.749–756) |
| 2608.15773 | Large values of quadratic character sums | B | Lemma 2.2 l.546–551: Σ_{d∈F, |d|≤X} χ_d(n) with error X^{1/2+ε}f(n₀)g(n₁) (n^ε-uniformity = Lindelöf on Re s = 1/2) |
| 2609.01101 | Shifted moments of Dedekind zeta | B | analytic, zero-free Re z > 1/2 for Φ-ratios (l.293–297), critical-line moments |
| 2607.21532 | Simultaneous non-vanishing, CLT | B | log L(s, χ) bounds near Re s = 1/2 (l.656–661, GS Lemma 8.2) |
| 2608.05961 | Moments of twisted quadratic characters | B | Soundararajan–Harper upper bounds (l.201–226) |
| 2608.15063 | Zeros near the central point, newforms | B | low-lying zeros on the line (l.234, 281) |
| 2608.06286 | Average analytic rank of y² = x³ − dx | B | RH for L(s, E_d) (l.222, 251); one-level density |
| 2607.12381 | Average rank, hyperelliptic curves | B | RH for L(s, C_a/K) (l.243, 2317) |
| 2607.09110 | Diagonal Elliott–Halberstam twisted by μ | B | GRH + simple zeros + Gonek–Hejhal (l.327, 1340) |
| 2609.14714 | Prime + non-square-free integer | 0 | GRH only in Lee–O'Clarey's cited theorem (l.33); paper is unconditional |
| 2608.12709 | Abelian division fields, prime levels | 0 | removes GRH (l.187) |
| 2607.22930 | Unshared zeros of Dirichlet L | 0 | unconditional (l.72, 198) |
| 2609.10769 | Symplecto-rigidity, Anosov | 0 | "[GRH]" is a citation key (l.67) |

## A1 — 2610.00109 (Hua, PageRank on LPS graphs): CONFIRMED

- Source: fix p ≡ 1 mod 4; admissible primes q (q ≡ 1 mod 4, (p/q) = −1, q > 2√p); X_q the LPS graph, N_q = q(q² − 1); Δ_q(α) ≥ 0 the
  excess return probability of the geometric-time stopped walk over the tree value; α_{q,c} = 1 − c/log q, c ∈ I = [c₀, c₁] ⊂ (0, ∞).
  Thm main-results (i) (unconditional): 1/N_q ≪ Δ_q(α_{q,c}) ≪ log log q/N_q; (ii) (l.216–222, GRH for L(s, χ_D), D < 0): Δ_q ≍_{p,I} 1/N_q.
- GRH step: only Lemma grh-short-euler (l.944–986: r₃(m) ≪ √m E_z(m) for m ≤ Y/4, z = (log Y)^8), from Granville–Soundararajan 2003
  Lemma 2.1 with s = 1, σ₀ = 1/2; used once, in the proof of (ii) with Y = q^6 (l.1380–1411: "uses GRH only for |D| ≤ q^6").
- Replacement: GS Lemma 2.1 needs only a zero-free rectangle {σ₀ < Re z ≤ 1, |Im z − t| ≤ y + 3}; take σ₀ = θ = θ₁₄₁ (L(s, χ) ≠ 0 for
  Re s > θ, every Dirichlet χ). For log z ≥ 2/(1 − θ) = 15.79, σ₁ = θ + 1/log z and the error is O(log|D| (log z)² e z^{θ−1}) =
  O((log Y)^{1−8(1−θ)} (log log Y)²), 1 − 8(1 − θ) = −0.0132331511…, uniformly in |D| ≤ Y: the lemma holds with the same z (8 > 1/(1 − θ) =
  7.8955174…; z = (log Y)^9 also works, exponent −0.1398872950…). The rest of (ii) is unconditional.
- Reviewer (subagent 7825853c): CORRECT, ROUTINE (checked σ₀ ≥ 1/2, |t| ≤ 2q, rectangle, absolute constant, H ≥ z^4 and the exceptional
  k-range O_p((log 2z)²) = o(log q)).

## A2 — 2607.21259 (fractional parts of x^m/p): CONFIRMED (relative to Hexagon 2610.00256v3)

- Source: f_m(p) = #{1 ≤ x < p/2 : {x^m/p} > 1/2}. Thm estimate_GRH (l.136–140): for even m ≥ 2 and every odd prime p, under GRH,
  f_m(p) = p/4 + O_m(√p log log p). Thm m=2 (unconditional): |f_2(p) − p/4| ≥ C√p log log p for infinitely many p ≡ 3 mod 4.
- GRH step: only Lemma 'dirichlet character' (l.204–209) = Montgomery–Vaughan 1977 Thm 2 (Σ_{n≤x} χ(n) ≪ √p log log p, χ ≠ χ₀ mod p),
  used in the proof l.399–404 (f_m(p) = (p − 1)/4 + ½ Σ_{χ^m = χ₀, χ ≠ χ₀} Σ_{y∈J} χ(y), J = (p/2, p)).
- Replacement: Hexagon 2610.00256v3 Thm every (Bhattacharjee–Bhattacharya; Harvest1 row 7): from the Dirichlet half-plane,
  M(χ) ≤ C_θ √q log log q for every primitive χ; mod a prime every χ ≠ χ₀ is primitive and Σ_{y∈J} χ(y) = −Σ_{n≤(p−1)/2} χ(n).
  So |f_m(p) − p/4| ≤ 1/4 + ((d − 1)/2) C_θ √p log log p, d = gcd(m, p − 1), unconditionally; sharp in order for m = 2.
- Reviewer (subagent 81d4c73d): CORRECT, ROUTINE (character expansion, interval reduction, primitivity; no other GRH use).

## A items that fail or were not pursued

- 2609.02392 (cubic/quartic L(1/2, χ), resonator): FAILS. Lemma PropDirpoly (l.339–349) does transfer (cubic/quartic residue symbols are
  Hecke characters of Kummer extensions, zero-free for Re s > 31/32 by Cureton–Thorner: error X^{31/32+ε} log(c + 2)), and S₂ survives
  (error X^{31/32+ε}|R| Σψ², |R| ≤ N = X^{o(1)}, l.740–755), but S₁ sums that error over the approximate functional equation, l ≤ (√X/A)^{1+ε}
  (l.575–580): X^{θ+1/4+ε}|R|Σψ² against the main term X·A_N Σψ², so it needs θ < 3/4.
- 2607.20853 (murmurations over Q(i)): FAILS. With a prime-sum exponent θ, the large-l part is ≪ X^{1/2+θ−δ+ε}/Z (l.232, 248–254) and the
  small-l dual part ≪ X^{θ−δ+ε} Z (R₂, R₃ from (k/ϖ) prime sums, l.298–331); the best Z = X^{1/4} (the paper's choice, l.254) needs
  δ > θ + 1/4: GRH gives the paper's δ > 3/4, while θ ≥ 3/4 (31/32, or θ₁₄₁ for the PIT of Q(i)) would need a window longer than X.
- 2608.02494 (isogeny graphs over Q(√d), d ∈ 𝒟): FAILS as a theorem. GRH enters through Banwait–Najman–Padurariu's determination of
  IsogDeg(Q(√d)) (l.5578–5592), i.e. Larson–Vaintrob / Type 2 least-prime bounds; our serrenumk isogeny-prime corollary removes GRH
  from the finiteness step, but its Type 2 bound log p < y*(k) ≤ max(3000n_M, 2 log Δ_M + n_M log 4) = 6000 here (n_M = 2) leaves
  primes p < e^{6000} that BNP's GRH bound had excluded: no unconditional IsogDeg, hence no unconditional isogeny graphs.
- 2609.16048 (Waring–Goldbach in short intervals): not pursued. ψ(x, χ) ≪ x^{1/2} log² x (l.329–331) becomes x^{θ₁₄₁} log²(qx), so the
  scale would need θ > θ₁₄₁ = 0.8733…, while the paper's minor-arc transform is stated only at θ = 1/2 + ε (l.339–341); unconditional
  zero-density methods already reach 19/24-type exponents for more variables: low value.
- 2608.23500 (logarithmic Chowla, GRH all-shift theorem): not pursued (time). GRH enters through the GRH clause of Klurman–Mangerel–
  Teräväinen Thm 1.5 (§16 l.106–110) and a twisted Liouville bound X^{1/2+2η} (§16 l.114–125), whose θ₁₄₁ analogue X^{θ+2η} holds;
  whether KMT's clause needs more than a half-plane is unchecked (5,000-line machine-structured paper). Candidate for a later round.
- Lesson: a GRH power saving that is summed over ≥ X^{1/4} family members or AFE terms (2609.02392, 2607.20853) does not transfer from
  any half-plane θ ≥ 3/4; only pointwise uses (2610.00109, 2607.21259) or uses with polylog/X^{o(1)} multiplicity (2607.01126) do.

## Ready COROLLARIES.md lines

- PageRank on Lubotzky–Phillips–Sarnak graphs without GRH (N11; GRHHarv3, 20 Nov 2026; routine; reviewer CORRECT, ROUTINE): for fixed
  p ≡ 1 mod 4 and I = [c₀, c₁] ⊂ (0, ∞), every sufficiently large admissible prime q has Δ_q(1 − c/log q) ≍_{p,I} 1/N_q uniformly in c ∈ I
  (Δ_q the excess return probability of the geometric-time stopped walk on the LPS graph X_q over the (p + 1)-regular tree, N_q = q(q² − 1)),
  given the Lean-checked θ₁₄₁ = 0.87334585610671857680… (L(s, χ) ≠ 0 for Re s > θ₁₄₁, every Dirichlet χ). Method: arXiv 2610.00109 (Hua,
  9 Sep 2026) proves it under GRH for L(s, χ_D), D < 0 (Thm main-results (ii), l.216–222); GRH enters only through Lemma grh-short-euler
  (l.944–986: Granville–Soundararajan 2003 Lemma 2.1 with σ₀ = 1/2, z = (log Y)^8, Y = q^6, l.1380–1411); with σ₀ = θ₁₄₁ the GS error is
  O((log Y)^{1−8(1−θ)}(log log Y)²), 1 − 8(1 − θ) = −0.0132331511…, uniformly in |D| ≤ Y, so the lemma holds with the same z. Before:
  1/N_q ≪ Δ_q ≪ log log q/N_q unconditionally (Thm (i)). Novelty: arXiv API (abs: PageRank, Jul–Nov 2026) and the 2026-10/11 listings:
  only 2610.00109; Zenodo 'Lubotzky-Phillips-Sarnak': nothing on it. [GRHHarv3](Riemann/Deep/GRHHarv3/REPORT.md)
- Fractional parts of x^m/p without GRH (N11; GRHHarv3, 20 Nov 2026; routine relative to Hexagon 2610.00256v3; reviewer CORRECT,
  ROUTINE): for every even m ≥ 2 and odd prime p, #{1 ≤ x < p/2 : {x^m/p} > 1/2} = p/4 + O_m(√p log log p), sharp in order for m = 2
  (the paper's unconditional Ω result), given θ₁₄₁ = 0.87334585610671857680…. Method: arXiv 2607.21259 (23 Jul 2026, rev. 1 Sep) proves
  it under GRH (Thm estimate_GRH, l.136–140); GRH enters only through Montgomery–Vaughan's Σ_{n≤x} χ(n) ≪ √p log log p (Lemma l.204–209,
  used l.399–404), which Hexagon 2610.00256v3 Thm every gives from the Dirichlet half-plane for every primitive χ (mod p every χ ≠ χ₀):
  |f_m(p) − p/4| ≤ 1/4 + ((gcd(m, p − 1) − 1)/2) C_θ √p log log p. Before: O_m(√p log p) unconditionally (the paper's even-polynomial
  theorem). Novelty: arXiv API (ti: fractional parts / abs: half-interval, Jul–Nov 2026) and 2026-10/11 listings: only 2607.21259; Zenodo
  'fractional parts of polynomials': nothing. [GRHHarv3](Riemann/Deep/GRHHarv3/REPORT.md)

## A3 — 2607.01126 (Selmer ranks in prime cyclic extensions): CONFIRMED after one repaired step

- Source: Thm intro (l.64: rank/Selmer proportions of cyclic p-twists of E/K, Gal(K(E[p])/K) ⊇ SL₂(F_p)), Thm intro_superelliptic (l.158:
  average #C_d(K) ≤ B(p, K) over the twists), Thm intro_hyperelliptic (l.167), all under ERH. ERH enters only through Thm
  effective_chebotarev (l.485–490: error c(|C|/|G|)X^{1/2}(log|Disc F| + n_F log X)), used at l.1146–1150 (fan ratios, M₁ = X^{2^{1−m_l}},
  c₁ log log X ≤ m_l ≤ c₂ log log X) and l.1691–1699 (Prop Yu+GRH, used only in Prop crucial1, l.1716; crucial1 used at l.1883).
- Replacement: Frei–Sofos–Widmer 2610.18491 Thm generaldelta with δ = 1/32 (Cureton–Thorner, every ζ_F): error 32c₁(|C|/|G|)X^{31/32}
  (log Δ_F + n_F log X). (1) ratio error M₁^{−1/32+ε} (0 < ε < 1/32), so the error terms become O(X^{(−1/32+ε)/(log X)^{c₂ log 2}} log log X)
  → 0 (l.1215, 1288, 1298). (2) Prop Yu+GRH: 32c₁ X^{31/32} B[K:Q] 2^{4g²ω(𝔡)}(log Nm 𝔡 + log X). (3) Prop crucial1 (Nm 𝔡 ≤ X, ω(𝔡) ≤
  2 log log X, l.1703): relative error O((log X)^{b} log f/f^{1/32}), b = 8g² log 2 + 1, so its hypothesis f(X) ≫ (log X)^{16g² log 2 + 3}
  becomes f(X) ≥ (log X)^{256g² log 2 + 48} (relative error O((log X)^{−1/2} log log X)); the fan threshold T = (log X)^A allows any A > 1
  (Remark parametrize_fan, l.925–945). Conclusion: Thms intro, intro_superelliptic, intro_hyperelliptic hold given the Lean-checked 31/32.
- Reviewer (subagent e1be53c5): GAP in (2) — log Nm 𝔡 must be (log X)^{O(1)} uniformly; otherwise CORRECT, ROUTINE (ratio bound from (b),
  (1) with ε < 1/32, (3) with A = 32b + 16). Repair: Prop Yu+GRH is used only inside Prop crucial1 (rg: l.1716), whose hypothesis is
  Nm 𝔡 ≤ X (l.1703), so log Nm 𝔡 ≤ log X; the paper's own discriminant bound |Disc F| ≤ (Nm 𝔡)^{[K:Q]B₃2^{4g²ω(𝔡)}} (l.1685–1688) then gives
  log Δ_F ≪ (log X)^{8g² log 2 + 1}. Gap closed (orchestrator re-check of this step advised).
- Ready line: Selmer ranks in prime cyclic extensions without ERH (N11; GRHHarv3, 20 Nov 2026; routine; reviewer CORRECT, ROUTINE after
  one repaired step): the KMR-type Selmer/rank distributions of 2607.01126 (Thms intro l.64, intro_superelliptic l.158 — average #C_d(K) ≤
  B(p, K) over cyclic p-twists of y^p = cubic —, intro_hyperelliptic l.167) hold given the Lean-checked 31/32 (Cureton–Thorner). Method:
  ERH enters only through effective Chebotarev (l.485–490); FSW generaldelta with δ = 1/32 gives X^{31/32}-errors; the fan ratios keep a
  stretched-exponential error, and Prop crucial1's prime threshold rises from (log X)^{16g² log 2 + 3} to (log X)^{256g² log 2 + 48}. Before:
  ERH-conditional (Klagsbrun–Mazur–Rubin's unconditional Markov model covers other orderings only). Novelty: arXiv API (abs: Selmer AND
  cyclic AND Markov/ERH, Jul–Nov 2026) only 2607.01126; 2026-10/11 listings: 2611.11413 (2-Selmer in progressions, other question);
  Zenodo 'prime cyclic extensions': nothing. [GRHHarv3](Riemann/Deep/GRHHarv3/REPORT.md)

RH relevance: N11 (GRH families: removes GRH/ERH from three 2026 statements using the Lean-checked half-planes; no new zero-free input).
