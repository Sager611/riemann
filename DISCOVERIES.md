# Discoveries

Gate (AGENTS.md): someone gets a new theorem, tool or barrier from it; dependents of one input count once. One line each; the linked page has the statement, what was known before, the review and the novelty checks.
See [RECORDS.md](RECORDS.md) (constants, ranges, exponents, proportions), [REFUTATIONS.md](REFUTATIONS.md) (broken claimed proofs), [COROLLARIES.md](COROLLARIES.md) (routine conditional deductions), [FORMALIZATIONS.md](FORMALIZATIONS.md) (machine-checked known mathematics).

## On the RH tree (node in parentheses; node list at the top of TREE.md)

- (N4) Sharp odd Robin criterion: RH ⟺ σ(n) < (e^γ/2) n log log n for every odd n > 277093800961598968947975375; Washington–Yang's gap closed. [details](discoveries/odd-robin.md)
- (N5) Báez-Duarte's natural-number Nyman–Beurling criterion holds in L^p for all 1 < p < 2, also for every Dirichlet L-function; unconditionally true for p ≤ 1.14497. [details](discoveries/nyman-beurling-lp.md)
- (N7, N11) Zero-density estimate for ring class L-functions, sharper in the discriminant than any published bound; the same method gives, for almost all discriminants, [least primes represented by forms ≤ |Δ|^{1+ε}](discoveries/least-primes-forms.md), [Sato–Tate in every ideal class](discoveries/sato-tate-classes.md) and [least primes in CM classes](discoveries/cm-least-primes.md). [details](discoveries/ring-class-zero-density.md)
- (N7 → N10, N11, N12) Consequences of the Lean-checked zero-free half-planes (one external input, counted once; [audited](Riemann/Deep/HPAudit/REPORT.md): statements in Mathlib's definitions, Comparator match, standard axioms, three kernels): [explicit Littlewood bound for every χ](discoveries/littlewood-half-plane.md), least primes [in progressions](discoveries/linnik-numeric.md) and [Chebotarev classes](discoveries/chebotarev-numeric.md), [unit-group covers](discoveries/unit-group-cover.md), [least prime primitive roots](discoveries/prime-primitive-root.md), primes [between seventh powers](discoveries/primes-between-powers.md) and [in (x, x + x^0.856]](discoveries/primes-short-intervals.md), [random polynomials irreducible without RH](discoveries/random-polynomials-without-rh.md), [root-discriminant bounds](discoveries/discriminant-bounds-strip.md), [Kopp–Lagarias Problem 1.7](discoveries/unit-generated-orders.md).
- (N8, N13) Jerby's Gram-point repulsion conjecture (|Z′(g_n)| > 4|Z(g_n)| at isolated bad Gram points) is false: certified counterexamples at n = 1,565,734,083 (μ = 3.968), 2,565,007,238, 3,315,629,581. [details](discoveries/jerby-repulsion.md)
- (N3) The Li-coefficient remainder: mean square N^{3/2}√log N under RH, Ω((n log n)^{1/4}) unconditionally, mean 7/4; RH ⟺ Σ m^δ y_m² < ∞ for some δ < 1/2. [details](discoveries/li-remainder.md)
- (N8) Lamzouri's small-gap criterion has a floor: it can never prove a positive proportion of gaps ≤ λ for λ ≤ 0.775 (certified LP dual); our record 0.7851 sits 0.0101 above. [details](discoveries/small-gap-floor.md)

## Side results (off the RH tree)

- Thirteen published conjectures of Z.-W. Sun fail in small cases, each certified by exhaustive search or exact arithmetic. [details](discoveries/sun-counterexamples.md)
- The Sun–Hou conjecture holds for every finite field F_q with q > 7 (prime powers up to 18,888,871 were open). [details](discoveries/sun-hou.md)
- A graph's evolution algebra is isomorphic to its random-walk version iff the graph is regular, or bipartite biregular with degree ratio a cube (open problem of arXiv 2610.22248). [details](discoveries/evolution-algebras.md)
- The open question of arXiv 2610.10095 answered: exactly 11 of the 24 cyclic orders of +, ×, −, ÷ make every number 1..N for all large N. [details](discoveries/operation-orders.md)
- Exactly which triples of prefix reversals generate S_n below the diagonal (m + k < n): a complete classification, refuting the conjecture of Blanco et al. [details](discoveries/prefix-reversals.md)
- Theorem: the Fasi–Higham positivity conjecture holds for m = 3 and every k (all Taylor coefficients of 1/q_{k,3} are positive). [details](discoveries/fasi-higham-m3.md)
- Gale's conjecture fails in R⁷ and R⁸ (the last open dimensions), hence in every dimension n ≥ 3; an explicit family with a hand proof for n ≥ 5. [details](discoveries/gale-r7-r8.md)
- Theorem: the diagonal Fasi–Higham positivity holds for every odd m (all Taylor coefficients of 1/q_{m,m} positive; computer-assisted). [details](discoveries/fasi-higham-diagonal.md)
- Johnston–Leong–Tudzi's conjecture holds: the normalised limsup s_a of Σ_{n≤x}(−a)^{Ω(n)} tends to 1 as a → ∞ (arXiv 2408.04143). [details](discoveries/omega-power-sums.md)
