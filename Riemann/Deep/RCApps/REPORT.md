DISCOVERY: for every A > 1.6536 (real root of A³−5A²+14A−14), all but O(X^{1−c}) negative fundamental discriminants |D| ≤ X have primes in every ring class in every interval (x, x + x|D|^{−κ}] with x ≥ |D|^A; so every form of discriminant D represents a prime ≤ 2|D|^A. The density-one record was |D|^{20/3+ε} with a density-zero exceptional set (Ditchen 2013). c = κ = 10^{−4} for A ≥ 1.7.

# RCApps report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Brief: `.scratch/orca/rcapps-brief.md`. Paper proofs, not refereed, no Lean. No 2026 half-plane is used.

## Theorems (complete proofs in the NOTE)
- **Prop 1 (tool).** Σ over fundamental D, X/2 < |D| ≤ X, of Σ_{χ²≠1} N(σ,T;χ) ≪ (XT)^{10(1−σ)/(3−2σ)+ε}. This is Humphries–Thorner's zero detector with our per-D second moment Σ_χ|L(½,χ)|² ≪ |D|^{1/2+ε} in place of convexity. The general GL₂ bound (Pascadi–Thorner 2025) gives X^{11}.
- **Thm 2 (targets 1 and 3).** For A > A₀ = 1.653595…, all fundamental D outside O(X^{1−c}) satisfy: every class C has ≥ ν_D H/(16 h log 2x) degree-one prime ideals of norm in (x, x+H], for all x ≥ |D|^A and H ≥ x|D|^{−κ}. Here |D|^{−δ/2−o(1)} ≤ ν_D ≤ 1, and ν_D = 1 unless a genus factor of conductor ≤ |D|^δ has a Siegel zero. For A ≥ 1.7: c = κ = 10^{−4}, δ = 1/2600.
- **Thm 3 (target 2).** Same, with the angle of the generator of 𝔭𝔠 in any arc of length ≥ |D|^{−κ}.
- **Cor 4.** All but O(N^{1−10^{−4}}) squarefree n ≤ N, n ≡ 1, 2 mod 4, satisfy: every interval (x, x + x(4n)^{−10^{−4}}] with x ≥ (4n)^{1.7} holds a prime a²+nb² with arg(a+ib√n) in any prescribed arc of length ≥ (4n)^{−10^{−4}}.

## How it works
- Zeros with β < A/2: our ring class density |D|^{2σ(1−σ)+ε}, for each single D ([Combos](../Combos/NOTE.md) Thm 1).
- Zeros of non-real characters with β ≥ A/2: Prop 1 and Markov's inequality over D. This step sets A₀. Kowalski–Michel gives 1.7242 (= [LeastPrimeForm](../LeastPrimeForm/NOTE.md) Thm B); Pascadi–Thorner gives 1.7055.
- Genus characters, L(χ_d)L(χ_{D/d}): factors with |d| ≤ |D|^δ are handled for every D by Thorner–Zaman's explicit log-free density with its Deuring–Heilbronn factor. A Siegel zero only scales the main term by ν_D ≥ |D|^{−δ/2}, which the power-saving errors beat. Factors with large |d| are averaged over D. This makes the exceptional set power-saving; LeastPrimeForm Thm B had density zero.

## Comparison
| Statement | Before | Now |
|---|---|---|
| least prime by every form, almost all D | \|D\|^{20/3+ε}, density-zero exceptions (Ditchen 1312.1502) | 2\|D\|^A, any A > 1.6536, O(X^{1−c}) exceptions |
| short intervals, sectors, uniform in D | none found | H ≥ x\|D\|^{−κ}, arcs ≥ \|D\|^{−κ} |
| least prime a²+nb² | n^{L₂}, L₂ unspecified, all n (Kowalski–Michel 2002) | 2(4n)^{1.7}, almost all n, also in short intervals and sectors |
| every D | \|D\|^{455} (Zaman thesis) | unchanged |

## Dependencies
Our ring class density and second moment ([Combos](../Combos/NOTE.md) Thm 1, Lemma 4; unrefereed, reviewed); Humphries–Thorner (Compos. Math. 2024, §5, Cor 4.4); Thorner–Zaman (arXiv 2208.11123, explicit log-free density); Density2 Thm A (sectors only).

## Novelty check
41 arXiv API queries; OpenAlex (Ditchen has 3 citing works, none improves him); 2 GXL claims; the Pascadi–Thorner and Humphries–Thorner sources; Hexagon 2610.00178; Kowalski–Michel 2002; Duke 1989. Nothing found below 20/3 for almost all D, and nothing with power-saving exceptional sets, short intervals or sectors uniform in D. Almost-all-fields Chebotarev results (Pierce–Turnage-Butterbaugh–Wood; Lemke Oliver–Thorner–Zaman; Lemke Oliver–Smith) need a fixed Galois group, but ring class fields have degree 2h → ∞.

## NOGO
- Every D (or |D| ≤ x^δ), and sectors for every n: these need an explicit log-free density near σ = 1 for each D. Unchanged.
- Shortening the mollifier with the per-D large sieve inside the detector costs a factor X in the tail term.

## Leads
- A GL₂ large sieve over {π_χ : |D| ≤ X} with threshold X^{3/2} instead of Q²|𝒮| = X^{7/2} would give A₀ = 1, the GRH exponent, for almost all D.
- Non-fundamental D. A larger explicit c (now limited by Thorner–Zaman's constants).

Significance: a new statement with a record exponent for almost all D, built from known tools; the reviewer should judge. It supersedes LeastPrimeForm Thm B. Two focused referee passes found no fatal error; the fixable gaps they raised are handled in the NOTE (§2.4(e), §2.5).
