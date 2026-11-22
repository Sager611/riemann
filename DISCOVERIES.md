# Discoveries

Gate (AGENTS.md): someone gets a new theorem, tool or barrier from it; dependents of one input count once. One line each; the linked page has the statement, what was known before, the review and the novelty checks.
See [RECORDS.md](RECORDS.md) (constants, ranges, exponents, proportions), [REFUTATIONS.md](REFUTATIONS.md) (broken claimed proofs), [COROLLARIES.md](COROLLARIES.md) (routine conditional deductions), [FORMALIZATIONS.md](FORMALIZATIONS.md) (machine-checked known mathematics).

Status (21 Nov 2026): assumes Alpöge's RH + Dirichlet GRH (Lean, 20 Nov 2026). ★ = RH-conditional before, unconditional now. [Masullo-risk] = weakens if Masullo's unverified Hecke-GRH is completed. All audits are internal agent audits, not human referees.

## Important, still alive after Alpöge (ranked; node in parentheses, list at the top of TREE.md)

1. (N5) ★ The continuous Nyman–Beurling distance has the conjectured exact rate D(λ)² log(1/λ) → Σ_ρ m_ρ²/|ρ|² (Báez-Duarte–Balazard–Landreau–Saias, Burnol). [details](discoveries/nyman-beurling-continuous-rate.md)
2. Gale's conjecture fails in R⁷ and R⁸ (the last open dimensions), hence in every dimension n ≥ 3; an explicit family with a hand proof for n ≥ 5. [details](discoveries/gale-r7-r8.md)
3. (N3) ★ The Li-coefficient remainder: mean square N^{3/2}√log N, Ω((n log n)^{1/4}) unconditionally, mean 7/4 (Voros's heuristic, proved); RH ⟺ Σ m^δ y_m² < ∞ for some δ < 1/2. [details](discoveries/li-remainder.md)
4. Exactly which triples of prefix reversals generate S_n below the diagonal (m + k < n): a complete classification, refuting the conjecture of Blanco et al. [details](discoveries/prefix-reversals.md)
5. (N8, N13) Jerby's Gram-point repulsion conjecture (|Z′(g_n)| > 4|Z(g_n)| at isolated bad Gram points) is false: certified counterexamples at n = 1,565,734,083 (μ = 3.968), 2,565,007,238, 3,315,629,581. [details](discoveries/jerby-repulsion.md)
6. Johnston–Leong–Tudzi's conjecture holds: the normalised limsup s_a of Σ_{n≤x}(−a)^{Ω(n)} tends to 1 as a → ∞ (arXiv 2408.04143). [details](discoveries/omega-power-sums.md)
7. Theorem: the Fasi–Higham positivity conjecture holds for m = 3 and every k, and on the diagonal for every odd m (all Taylor coefficients of 1/q_{k,m} positive; diagonal computer-assisted). [m = 3](discoveries/fasi-higham-m3.md), [diagonal](discoveries/fasi-higham-diagonal.md)
8. (N7, N11) Zero-density estimate for ring class L-functions, sharper in the discriminant than any published bound; the same method gives, for almost all discriminants, [least primes represented by forms ≤ |Δ|^{1+ε}](discoveries/least-primes-forms.md), [Sato–Tate in every ideal class](discoveries/sato-tate-classes.md) and [least primes in CM classes](discoveries/cm-least-primes.md). [details](discoveries/ring-class-zero-density.md) [Masullo-risk]
9. (N8, N12) Unconditional proportions of critical-line zeros of linear combinations of Dirichlet L-functions: 34.7% for two, 2.06% for three, 6.8·10⁻⁴ for four, ≥ 3.4·10⁻³/N for N ≥ 11; N = 2 is only 0.2% above Bombieri–Hejhal's GRH value 34.5%, now unconditional. [two](discoveries/linear-combinations-line.md), [three](discoveries/linear-combinations-three.md), [N ≥ 4](discoveries/selberg-linear-combinations.md)
10. (N8) Below √T: for some θ < 1/2, a positive proportion of zeros in every window (T, T + T^θ] are simple and critical, the first crossing of √T (was θ ≥ 0.527; thresholds in RECORDS.md). [details](discoveries/simple-zeros-short-windows.md) [Masullo-risk: he claims simple non-central zeros, unverified]
11. (N8, N12) Every GL(2) newform L-function has ≥ 7.005% of its zeros on the line, and the one-piece Levinson–Conrey method can never exceed 7.005633% (our record sits 10⁻⁷ below). [record](discoveries/gl2-critical-line.md), [ceiling](discoveries/gl2-levinson-ceiling.md) [Masullo-risk: CM and dihedral cases]
12. Thirteen published conjectures of Z.-W. Sun fail in small cases, each certified by exhaustive search or exact arithmetic. [details](discoveries/sun-counterexamples.md)
13. The open question of arXiv 2610.10095 answered: exactly 11 of the 24 cyclic orders of +, ×, −, ÷ make every number 1..N for all large N. [details](discoveries/operation-orders.md)
14. (N7 → N11, N12) Via the Cureton–Thorner Dedekind-ζ strip (31/32; input unrefereed), not covered by Alpöge: [random polynomials irreducible without RH](discoveries/random-polynomials-without-rh.md), [root-discriminant bounds](discoveries/discriminant-bounds-strip.md) [Masullo-risk]
15. A graph's evolution algebra is isomorphic to its random-walk version iff the graph is regular, or bipartite biregular with degree ratio a cube (open problem of arXiv 2610.22248). [details](discoveries/evolution-algebras.md) Concurrent: Liu, Zenodo 23608087 (uploaded 24 Oct 2026; ours 22 Oct), same criterion.

## Superseded by Alpöge (follow from RH/GRH; kept as records, see [SUPERSEDED.md](SUPERSEDED.md))

- (N4) Sharp odd Robin criterion: RH ⟺ σ(n) < (e^γ/2) n log log n for every odd n > 277093800961598968947975375; Washington–Yang's gap closed. [details](discoveries/odd-robin.md)
- (N5) Báez-Duarte's natural-number Nyman–Beurling criterion holds in L^p, 1 < p < 2, up to the exact endpoint (ζ ≠ 0 on Re s > 1/p suffices), also for Dirichlet L-functions. [details](discoveries/nyman-beurling-lp.md) (the non-abelian ζ_K version, NBLpK in COROLLARIES, is not covered)
- (N8) Lamzouri's small-gap criterion has a floor: it can never prove a positive proportion of gaps ≤ λ for λ ≤ 0.775 (certified LP dual); smallgap’s 0.7851 sat 0.0101 above. [details](discoveries/small-gap-floor.md)
- (N8) A hard-core inequality on consecutive gaps crosses that floor: summed along the gap-free zero set, it proves a positive proportion of gaps ≤ 0.768, unconditionally. [details](discoveries/small-gap-hardcore.md)
- (N10, N11) Zero-free half-plane consequences: [unit-group covers](discoveries/unit-group-cover.md), [Kopp–Lagarias Problem 1.7](discoveries/unit-generated-orders.md), [least prime primitive roots](discoveries/prime-primitive-root.md) (constants in RECORDS.md).

## Published first by others (not new)

- The Sun–Hou conjecture holds for every finite field F_q with q > 7. [details](discoveries/sun-hou.md) Published first: Zhao, Zenodo 23202607, 7 Oct 2026 (ours 22 Oct); drop under the gate once confirmed.
