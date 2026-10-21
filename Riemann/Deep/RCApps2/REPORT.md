DISCOVERY: for every A > 3/2, all but O(X^{1−c}) fundamental discriminants D with |D| ≤ X, negative or positive, have primes in every ideal class (narrow if D > 0) in every interval (x, x + x|D|^{−κ}] with x ≥ |D|^A; so every primitive binary quadratic form of discriminant D, definite or indefinite, represents a prime ≤ 2|D|^A. Before: A > 1.6536 for D < 0 (RCApps); in the literature |D|^{20/3+ε} with density-zero exceptions (Ditchen 2013, definite forms only); nothing on average for indefinite forms.

# RCApps2 report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Brief: `.scratch/orca/rcapps2-brief.md`. Paper proofs, not refereed, no Lean. Unconditional except the marked conditional rows.

## Theorems (complete proofs in the NOTE)
- **Prop 1 (tool).** For every family covered by Density2 (ring class characters of any imaginary quadratic order, anticyclotomic characters, narrow ring class characters of real quadratic orders), Σ_χ N(σ;t,χ) ≪ (𝔔T)^{2(1−σ)/(3−2σ)+ε} for 0.6 ≤ σ ≤ 1. The old exponent 2σ(1−σ) (Combos Thm 1, Density2 A–C) is larger for every σ in (½, 1). Method: a mollifier as long as the per-D large sieve threshold 𝔔^{1/2}, and the first-moment zero detector of Humphries–Thorner; the k-th power step disappears.
- **Thms 2′, 3′, Cor 4′ (D < 0).** RCApps Thms 2, 3 and Cor 4 hold for every A > 3/2: short intervals, sectors, primes a² + nb². For A ≥ 1.55: κ = 8·10^{−5}, exceptional set O(X^{1−1.7·10^{−4}}).
- **Thm 4, Cor 5 (D > 0, new).** For every A > 3/2 and all fundamental D ∈ (0, X] outside O(X^{1−c}): every narrow class holds ≫ H/(D^{1/2+δ}log x) prime ideals of norm in (x, x+H] for x ≥ D^A, H ≥ xD^{−κ}. So every primitive indefinite form represents a prime ≤ 2D^A, and the least prime a² − nb² is ≤ 2(4n)^A for almost all n.

## How it works
- Zeros with β < 3/2 − 1/A: Prop 1 for each D (before: β < A/2). The rest: the family density X^{10(1−σ)/(3−2σ)} (RCApps Prop 1; Prop 2 for D > 0) and Markov over D. At the handover the Markov condition reads 2(2−A) < 1, i.e. A > 3/2.
- At A = 3/2, σ = 5/6 three exponents meet: X times the per-D bound, the family bound and the Markov bar all equal X^{5/4}.
- D > 0: narrow class characters give Maass forms of eigenvalue ¼ (tempered), so the Humphries–Thorner detector applies; genus characters and Siegel zeros are handled exactly as for D < 0.

## Comparison
| Statement | Before | Now |
|---|---|---|
| per-D zero density, class group families | \|D\|^{2σ(1−σ)+ε}T³ per window (Combos) | (\|D\|T)^{2(1−σ)/(3−2σ)+ε} |
| least prime, every form, almost all D < 0 | 2\|D\|^A, A > 1.6536 (RCApps); \|D\|^{20/3}, density-zero exceptions (Ditchen) | 2\|D\|^A, any A > 3/2, O(X^{1−c}) exceptions |
| least prime, indefinite forms, almost all D > 0 | none found; every D: D^{O(1)} (Weiss), about D^{455} (Zaman) | 2D^A, any A > 3/2 |
| conditional on Liu's half-plane (L): least prime, every D < 0 | \|D\|^{1.7498} (Combos) | \|D\|^{1.5998}; \|D\|^{8/5} under Hu's 7/8 |
| conditional on (L): ring class prime number theorem | \|D\| ≤ x^{0.5715} | \|D\| ≤ x^{0.6251} |

## Below 3/2: NOGO with current tools, and leads
- Needed: a family bound that beats the sum of the per-D bounds by some X^η near σ = 5/6, i.e. a large sieve over all D ≤ X at length about X^{3/4} that beats X·N. Class II is already sharp. The Humphries–Thorner sieve on k-th powers, intermediate mollifier lengths and Rankin–Selberg bounds give nothing (`hyb.py`).
- Payoffs (`hyp.py`): a saving X^{0.01} gives A₀ ≈ 1.47; a family fourth moment X^{3/2+θ}, θ < 1/4, gives A₀ < 3/2; the optimal family sieve up to length X^{1−η} gives A₀ = 1, the GRH exponent.
- Lead: pairs of ideals in one class come from t² + |D|y² = 4nn′, so summing over D gives Hurwitz class number sums, which the Eichler–Selberg trace formula turns into Maass-form Hecke eigenvalues. Open step: match the main terms with the genus characters.

## Novelty check
34 arXiv API queries and 6 evidence-API claims (files in `.scratch/rcapps2`). Relevant: Ditchen 1312.1502 (definite forms), Zaman 1506.01635 (exceptional-zero case), Sardari 1802.06193 and 1803.03218 (GRH), Gaudet 2408.01610, Pascadi–Thorner 2508.14888, Blomer–Thorner 2410.17158 (GL(m) families, not class groups). No average result for indefinite forms; no per-D density for class group families below 2σ(1−σ).

## Significance (for the reviewer)
A new record exponent 3/2 for both signs of D. The method is a known detector (Humphries–Thorner) fed with this project's per-D large sieves, so it is classical; the real-quadratic statement is new. The conditional exponents (1.5998, 8/5, 0.6251) are routine corollaries for COROLLARIES.md. Supersedes RCApps Thm 2 and Combos Thm 1.
