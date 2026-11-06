CANDIDATES: 10 (best: hexagon:2610.00545, node N10, first numeric Linnik bound and explicit uniform power-saving prime number theorem in progressions from the Lean-checked half-plane)

Scope: arXiv listings math.NT/CA/CV/SP for 2026-10 and 2026-11 (1855 entries, 994 math.NT with cross-lists) → 106 papers; Hexagon API (716 works) → 39 papers; 145 TeX sources scanned with the brief's phrase list plus 10 similar phrases (`.scratch/harvest1/scan.py`); 119 had hits (171 strong phrases, 708 "sufficiently large"/"absolute constant"), all kept in `.scratch/harvest1/hits/`. Each hit was checked against TREE, COROLLARIES, DISCOVERIES, NOGO, REFUTATIONS, Riemann/Tree and Riemann/Deep. Harvest agent harvest1, 2026-11-06, 02:12–02:56 UTC.

## Ranked candidates

| # | id (source line) | node | statement in the source | what is missing | inputs (explicit where?) | h | gain | risk |
|---|---|---|---|---|---|---|---|---|
| 1 | hexagon:2610.00545 Thm 1 l.36 (Pollack); hexagon:2610.00018 l.134–144 (Naslund); arXiv 2610.29018 l.267 (Platt–Trudgian); Zenodo 23608576 | N10, N11 | Linnik: the least prime p ≡ a (q) is ≤ q^L for every q, "L absolute, effectively computable"; Naslund: q^{3.99} for q ≥ q_0, "we do not compute q_0 or C"; Platt–Trudgian: the explicit PNT in progressions of Bennett et al. "could be improved almost instantly"; Zenodo: π(x;q,a) = Li(x)/φ(q) + O(x^{7/8} log²(qx)), constant unspecified | any numeric bound: no explicit power-saving ψ(x;q,a) bound uniform in q, no numeric Linnik bound valid for all q | Bebel's Lean-checked half-plane Re s > B_2 + 2.2·10⁻²⁵⁴ for ζ and every Dirichlet L (QuasiRiemannTracker); N(T,χ) of Bennett–Martin–O'Bryant–Rechnitzer (Math. Comp. 2021, arXiv 2005.02989, refereed); smooth explicit formula (our ArtinHalfPlane NOTE §1; arXiv 1802.00085; Fiori–Jaskari 2609.23222 for ζ) | 5 | first numeric \|ψ(x;q,a) − x/φ(q)\| ≤ c x^{0.8734} log²(qx) for all q, x; first numeric Linnik bound p ≤ (c q log q)^{7.90} for every q; bonus: first numeric polylog least non-residue (c log q)^{7.90} (OpenAI QRH 30 Sep Cor 1.2: C(log p)^{32}, C unspecified) | medium: non-explicit forms are public (Zenodo 23608576; OpenAI 5 Oct, x^{11/12}); Platt–Trudgian and Tao's explicit-estimates project named this upgrade |
| 2 | arXiv 2610.18491v2 Thm 1.1 l.194–206, Cor 1.2 (Frei–Sofos–Widmer); hexagon:2610.00216 Cor chebotarev l.4690–4702 (Hu) | N11 (N12) | under a zero-free half-plane Re s > 1 − δ for ζ_L, every Chebotarev class has a prime of norm ≤ c(δ)(log Δ_L)^{1/δ}(log log Δ_L)^{2/δ}; c_1, c(δ) "effectively computable"; Hu: q_𝔭 ≪ 𝒟_L^8, constant "absolute and effectively computable" | the numeric constant; no explicit polylog least-prime-ideal bound exists (explicit exponents only: 40 Zaman arXiv 1508.00287, Ahn–Kwon; 16 Kadiri–Ng–Wong arXiv 1902.08640) | Cureton–Thorner hexagon:2610.00200v2 l.172/179: 31/32 for every ζ_K and Hecke L, Lean AxiomMath/DedekindZetaZeroFree (confirmed by Zenodo 23373237); Liu hexagon:2611.00082 B_2 (unrefereed, our Hexagon2 audit); explicit N_K(T) (Hasanalizade–Shen–Wong, Math. Comp. 2022); our ArtinHalfPlane NOTE §3.3 (proof with unspecified C*) | 6 | first numeric polylog bound: N𝔭 ≤ (c log 3Δ_L)^{32} from a Lean-checked input, (c' log 3Δ_L)^{7.90} with Liu | medium: Thorner–Zaman school has explicit Lagarias–Odlyzko tools; FSW v3 |
| 3 | hexagon:2610.00256v3 Thm every l.119, Prop logderiv l.870–877 (Bhattacharjee–Bhattacharya) | N11 | \|L'/L(1,χ)\| ≤ C_5 log log q for every nonprincipal χ, unconditionally; C_5 = 6/(1−σ) + 429K_1(1 + 1/(1−σ)) + 99 with K_1 = (1+2A_0)/(2δ²), δ = 1/2000, A_0 = 1.5ζ(3/2): C_5 ≈ 1·10¹¹ (my evaluation; never stated) | a usable constant | Bebel's half-plane; our Littlewood78 smoothed explicit formula (PROOF.md step 3 already bounds Re L'/L(1,χ)); Dusart arXiv 1002.0442 | 2 | \|L'/L(1,χ)\| ≤ (7.90 + o(1)) log log q + c (GRH gives 2, Ihara–Murty–Shimura) | low; routine (COROLLARIES level) |
| 4 | arXiv 2610.26841v2 Thm 1.1 l.211–216, Rem l.223 (Lamzouri) | N8 | some c > 0 makes ≥ 0.6552 of zeros c-isolated (hence simple and on the line); explicit only c = 1/6 → 0.50875; "we did not try to optimize" | an explicit c for the 65.5% statement | his Prop IsolatedSMALL l.511–517 (the O(τ²N) term is not tracked); constant line l.816 (2 − C_MT − 0.0806α_0 > 0.655202) | 3 | first explicit c (table c ↦ proportion) | low–medium (v3) |
| 5 | hexagon:2610.00519v1 l.60, Cor dVP l.118 (Vega; AI-generated, unrefereed) | N4 | every superabundant Robin counterexample n has log log n − log p_k < C exp(−c√log p_k), C, c absolute | numeric C, c | Fiori–Kadiri–Swidinsky explicit ψ bound; Dusart θ bounds; our RobinCH range (counterexamples exceed 10^(10^46.79)) | 2 | first explicit constants; small significance | Vega's reduction itself unverified |
| 6 | arXiv 2610.38870v1 Thm 1 l.121 (Ono) | N6 | J_d (shift 0, ξ) has ≥ d − C√(d log d) distinct negative real roots for every d ≥ 2, C absolute | numeric C (saddle lemmas l.593, 622, 746 use unspecified absolute c) | the paper's saddle analysis; ξ Taylor coefficients; our record (hyperbolic for d ≤ 1.61·10²⁵) | 6+ | first explicit C; matters only for d > 1.61·10²⁵ | may overrun 6 h |
| 7 | hexagon:2610.00256v3 Thm every (Bhattacharjee–Bhattacharya) | N11 (N9) | M(χ) ≤ C_θ√q log log q for every primitive χ, unconditionally; constants by formula only (checking programs in Zenodo 23672710) | a usable numeric C_θ | transfer theorem, Appendix A zero density (c_0 = 10¹²), half-plane | 6+ | usable constant; beats explicit Pólya–Vinogradov (2/π²)√q log q only for astronomically large q | high: constant may stay astronomical; authors may print numbers |
| 8 | arXiv 2610.16010v1 l.1616 (Trudgian–Zuniga-Alterman) | N12 | "ζ_K has a zero below height 14 for every K with n_K ≥ 2 … we leave this as a task for future research" | the proof for all K | their linear kernel with Stechkin's trick; Odlyzko's unconditional discriminant bounds; LMFDB complete tables (small degrees only) | 6+ | universal statement | high: degrees 10 to 10³ may not close |
| 9 | arXiv 2610.29018v1 l.269 (Platt–Trudgian) | N7, N11 | Kadiri's zero-free region (R = 5.60) "is only valid for q ≤ 4·10⁵"; "one could hope for a version … that holds for all q" | nothing hard: instant from inputs | Kadiri; their Thm turbo (no real zeros, q ≤ 10⁹ even, 10¹⁰ odd); 2610.21511 Thm gap (1 − β ≥ 5811·10⁻⁷/log q) | 1 | exception-free region for every χ: R = 5.60 for q ≤ 10⁹, real zeros outside 1/(1721 log q) beyond; fold into #1 | low gain |
| 10 | arXiv 2610.28822v1 l.422 (Fouvry–Shparlinski–Xi) | N12 | Möbius sums over F_p[t] polynomials with coefficients in boxes save p^{−δ}; "δ can be made explicit" | numeric δ(d, γ) | explicit Weil/Burgess bounds | 2 | first explicit δ | low RH relevance |

## Mini-briefs (top 5)

**1. Explicit progressions (N10, N11).** Goal: numeric c_1, c_2, c_3 with |ψ(x;q,a) − x/φ(q)| ≤ c_1 x^B log²(qx) for all x ≥ 2, q ≥ 1, (a,q) = 1, B = 0.87338196661320752; least prime p ≡ a (q) ≤ (c_2 q log q)^{1/(1−B)}; least n with χ(n) ∉ {0, 1} ≤ (c_3 log q)^{1/(1−B)}.
Inputs: Bebel's Lean theorem (github.com/GTimG/QuasiRiemannTracker; Riemann/Deep/Hexagon2/QRT-chronology.md); BMOR N(T,χ) Thm 1.1 as coded in Riemann/Deep/DirLi/certify.py; ArtinHalfPlane NOTE §1 (uniform smooth explicit formula); arXiv 1802.00085 and 2609.23222 for the unsmoothing; Dusart/Rosser–Schoenfeld for θ versus ψ.
Method: per primitive χ a smoothed explicit formula; zero sum ≤ x^B Σ|K(ρ)| by partial summation over BMOR counts; low zeros use the symmetry β ≥ 1 − B; unsmooth by monotonicity; reduce imprimitive χ; orthogonality; Linnik by θ(x;q,a) > 0; non-residue by Σ Λ(n)(1 − Re χ(n))k(n/x) > 0.
Certificate: Arb, `uv run --no-project --with python-flint python3 cert.py`; sanity against exact ψ(x;q,a) for q ≤ 100, x ≤ 10⁸ and least primes for q ≤ 10⁴.
Verify: reviewer re-derives every truncation and imprimitivity term; adversarial cases x < q, q ≤ 4, prime powers.
Novelty: lst.py 'progression|Linnik|least prime|non-?residue' (Oct–Nov), Hexagon API abstracts, Zenodo ('"arithmetic progressions" AND "quasi-Riemann"': 23608576 is non-explicit), our NOVELTY files, Thorner–Zaman 1604.01750 and 2108.10878, Zaman 1506.01635 (explicit constants only under a Siegel zero), BMOR 1802.00085, GRH-only 2404.08380.

**2. Explicit polylog Chebotarev (N11).** Goal: numeric c, c' such that every class C of Gal(L/K) contains Frob_𝔭 of an unramified degree-one 𝔭 with N𝔭 ≤ (c log 3Δ_L)^{32} (Cureton–Thorner 31/32, Lean-checked), and ≤ (c' log 3Δ_L)^{7.90} with Liu's B_2.
Inputs: hexagon:2610.00200v2 l.172/179 and the Lean repo AxiomMath/DedekindZetaZeroFree; hexagon:2611.00082 (Riemann/Deep/Hexagon2 audit); Hasanalizade–Shen–Wong explicit N_K(T); Lagarias–Odlyzko 1977 Deuring reduction; Bach–Sorenson's kernel (Math. Comp. 1996); ArtinHalfPlane NOTE §3.3 (l.261–296, the proof with C*).
Method: smoothed sum over primes of E = L^{⟨g⟩}; explicit formula over the zeros of ζ_L; zero sum ≤ x^{1−δ} Σ_ρ|K(ρ)| by the Hadamard identity (as in Littlewood78) or by HSW counts; explicit ramified, prime-power and higher-degree terms; Minkowski for n_L.
Certificate: Arb constants; PARI least primes in classes of Q(ζ_q) and S_3 cubic fields as sanity.
Verify: reviewer on the Deuring step and uniformity in n_L; novelty: FSW v2, Hu Cor chebotarev, Kadiri–Ng–Wong, Zenodo '"Chebotarev" AND explicit', arXiv listings 'Chebotarev|prime ideal'.

**3. Explicit L'/L(1,χ) (N11).** Goal: |L'/L(1,χ)| ≤ a log log q + b for every nonprincipal χ mod q ≥ 3 (primitive first; imprimitive adds Σ_{p|q} log p/(p − 1)), a near 1/(1 − B) = 7.90.
Inputs: Riemann/Deep/Littlewood78 (PROOF.md steps 1–3, cert.py); Bebel's B; Dusart 1002.0442 Thms 6.11–6.12; source constant ≈ 10¹¹ from hexagon:2610.00256 Prop logderiv.
Method: Littlewood78's smoothed explicit formula at s = 1 for −L'/L, real and imaginary parts, zero sum through the Hadamard identity (no zero counting), x ≈ (log q)^{1/(1−B)}.
Certificate: Arb cert.py variant; sanity: PARI lfun values of L'/L(1,χ) for q ≤ 10⁴ against the bound.
Verify: reviewer; novelty: Zenodo 23672710 (formula constants only), Ihara–Murty–Shimura (GRH), listings 'logarithmic derivative|Euler.Kronecker'.

**4. Explicit isolated-zero constant (N8).** Goal: an explicit c_0 with liminf B_0(T,c_0)/N(T) ≥ 0.655, plus the table c ↦ proportion for c ∈ [0.001, 1/6].
Inputs: arXiv 2610.26841v2 (Prop IsolatedSMALL l.511–517, Montgomery–Taylor kernel, l.816 constants C_MT and α_0 and their sources).
Method: bound the O(τ²N) term explicitly (|K(x) − 1| ≤ 2π²x²∫u²η²(u)du) and the O(ε) of the modified kernel ψ_ε (l.814, α_ε = α_0 + O(ε)), check that no other term depends on c, optimise the kernel per c.
Certificate: Arb evaluation of the proportion at each c, including α_0 < 0.21462, which the paper takes from a Mathematica maximisation (l.811).
Verify: reviewer checks the τ–c relation and that no hidden c-dependence remains; novelty: Lamzouri v3, listings 'isolated|gaps between zeros|moderate', our N8 files.

**5. Explicit superabundant constants (N4).** Goal: numeric C, c (and c_2) in Vega's Cor windowL and Cor dVP.
Inputs: hexagon:2610.00519v1 l.60, 118, 424, 502; Fiori–Kadiri–Swidinsky |ψ(x) − x| < 9.22022 x(log x)^{3/2}e^{−0.8476836√log x}; Dusart; Riemann/Deep/RobinCH range.
Method: first re-derive Vega's reduction (AI-generated, unrefereed; a flaw goes to REFUTATIONS), then carry explicit constants and use our range for p_k.
Certificate: short Arb script.
Verify: reviewer verdict on the reduction before any constant; novelty: Vega's other Hexagon items (2610.00702), listings 'superabundant|Robin'.

## Dropped hits
- hexagon:2610.00018 q_0 for q^{3.99}: Heath-Brown/Xylouris machinery with 4453 root cases, far beyond 6 h.
- hexagon:2610.00466 ε_K "in principle effective": done by bsexpl ((6 − 4√2)/(2^{n−1}n!), COROLLARIES); arXiv 2610.31402 c(k): bsexpl; 2610.21511: imqclass.
- hexagon:2610.00661/00688 (Suzuki criterion to 2.15·10²⁹): the limit comes from the zero-density bin 0.80 < β ≤ 0.85 (their l.207), below the half-plane; no gain.
- arXiv 2610.27830: making Baker–Harman–Pintz explicit is a full sieve project.
- GRH or conjecture-conditional: 2610.14257, 2610.36464, 2610.42022, 2610.21655, 2610.32672.
- Already explicit: hexagon:2610.00207 (70), 2611.00080 (1/1370), 2610.00702, 2610.00634 (41.007).
- Off-tree or research-level: 2610.32012 (Nowak/Lao ideal counts), 2611.03356 (Beurling systems), 2610.04045, 2610.11974, 2610.18491 l.263 (log log removal), 2610.26575, 2611.01905, 2610.41394.

## Papers scanned
2610.00386 none · 2610.00599 hit · 2610.01035 hit · 2610.01524 hit · 2610.01998 none · 2610.02709 hit
2610.03707 hit · 2610.03728 hit · 2610.04045 hit · 2610.04678 hit · 2610.04972 hit · 2610.05227 hit
2610.06136 hit · 2610.06471 none · 2610.08569 hit · 2610.08965 hit · 2610.09567 hit · 2610.10584 none
2610.10595 none · 2610.10609 hit · 2610.11024 none · 2610.11487 none · 2610.11974 hit · 2610.12234 hit
2610.12803 none · 2610.12816 hit · 2610.13136 none · 2610.13565 none · 2610.13864 hit · 2610.14257 hit
2610.14628 none · 2610.15140 hit · 2610.15291 hit · 2610.15316 hit · 2610.15483 hit · 2610.15559 hit
2610.15784 hit · 2610.15989 hit · 2610.16010 hit · 2610.16112 hit · 2610.16953 hit · 2610.17144 none
2610.17349 none · 2610.17720 none · 2610.18065 hit · 2610.18214 hit · 2610.18491 hit · 2610.18923 hit
2610.20058 hit · 2610.20214 hit · 2610.20792 hit · 2610.21099 hit · 2610.21511 hit · 2610.21730 hit
2610.21800 hit · 2610.21952 hit · 2610.22107 hit · 2610.22291 hit · 2610.22927 none · 2610.22939 hit
2610.23201 hit · 2610.23467 hit · 2610.23660 hit · 2610.23985 none · 2610.24612 hit · 2610.24958 hit
2610.26038 hit · 2610.26408 none · 2610.26575 hit · 2610.26705 hit · 2610.26835 hit · 2610.26841 hit
2610.27343 hit · 2610.27493 none · 2610.27591 hit · 2610.27804 hit · 2610.27830 hit · 2610.28527 hit
2610.28822 hit · 2610.29018 hit · 2610.30476 hit · 2610.31402 hit · 2610.31744 hit · 2610.32012 hit
2610.32694 hit · 2610.32704 hit · 2610.36464 hit · 2610.36885 hit · 2610.37371 hit · 2610.38870 hit
2610.39120 hit · 2610.40267 hit · 2610.41115 none · 2610.41394 hit · 2610.42022 hit · 2611.00525 none
2611.00924 hit · 2611.01905 hit · 2611.02058 hit · 2611.02071 hit · 2611.02675 hit · 2611.02933 hit
2611.02944 hit · 2611.03356 hit · 2611.05000 hit · 2611.06245 none
hx:2610.00018 hit · hx:2610.00161 hit · hx:2610.00178 hit · hx:2610.00200 hit · hx:2610.00207 hit · hx:2610.00216 hit
hx:2610.00220 hit · hx:2610.00236 hit · hx:2610.00256 hit · hx:2610.00275 none · hx:2610.00311 hit · hx:2610.00322 hit
hx:2610.00343 none · hx:2610.00346 hit · hx:2610.00375 hit · hx:2610.00385 hit · hx:2610.00409 none · hx:2610.00449 hit
hx:2610.00458 hit · hx:2610.00461 hit · hx:2610.00466 hit · hx:2610.00504 none · hx:2610.00519 hit · hx:2610.00521 hit
hx:2610.00535 hit · hx:2610.00545 hit · hx:2610.00584 hit · hx:2610.00634 hit · hx:2610.00655 hit · hx:2610.00661 hit
hx:2610.00688 hit · hx:2610.00702 hit · hx:2610.00714 none · hx:2610.00716 hit · hx:2610.00721 hit · hx:2611.00059 hit
hx:2611.00080 hit · hx:2611.00082 hit · hx:2611.00092 hit
("hit" = at least one phrase match; "hx:" = Hexagon id, not on arXiv.)

RH relevance: N10, N11, N8, N4, N6, N12, N7
