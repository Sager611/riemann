# Superseded by Alpöge's proof (21 Nov 2026)

Assumes Levent Alpöge's RH + GRH for all Dirichlet L-functions (Lean 4, 20 Nov 2026; checked by Sutherland, the Lean/Mathlib communities) is correct. These results now follow from it, or have better versions under it. Pages are kept as records; nothing is deleted. Not covered by it (still alive): Hecke/Dedekind of non-abelian fields, GL(2), linear combinations, simple zeros, side results.

Sources of the comparisons: Riemann/Deep/PostRH/CONDITIONAL.md §C and the 21 Nov 2026 ranking.

| Page | Group | Why superseded |
|---|---|---|
| [small-gap-floor](discoveries/small-gap-floor.md) | Small gaps | RH gives a positive proportion of gaps < (1/2+ε) of the mean (Durkan–Pearce-Crump, arXiv 2610.28527); method floor kept as a diagnostic. |
| [small-gap-hardcore](discoveries/small-gap-hardcore.md) | Small gaps | same: RH gives λ > 1/2 (2610.28527). |
| [small-gap-seven-elevenths](discoveries/small-gap-seven-elevenths.md) | Small gaps | same: RH gives λ > 1/2 (2610.28527). |
| [small-gap-two-thirds](discoveries/small-gap-two-thirds.md) | Small gaps | same: RH gives λ > 1/2 (2610.28527). |
| [large-gaps-higher-order](discoveries/large-gaps-higher-order.md) | Large gaps | RH gives Λ > 3.63846723 (2610.28527) vs our 3.367. |
| [robin-nicolas-ranges](discoveries/robin-nicolas-ranges.md) | Robin/Nicolas | RH ⇒ Robin's inequality for all n > 5040, Nicolas's for all primorials. |
| [odd-robin](discoveries/odd-robin.md) | Robin/Nicolas | criterion moot; the odd inequality is now a consequence of RH (a GitHub draft stated the threshold on 4 Oct). |
| [jensen-kim-lee-wedge](discoveries/jensen-kim-lee-wedge.md) | Jensen | RH ⇔ every Jensen polynomial of ξ is hyperbolic (Pólya). |
| [jensen-two-thirds-wedge](discoveries/jensen-two-thirds-wedge.md) | Jensen | same. |
| [jensen-wedge-constant](discoveries/jensen-wedge-constant.md) | Jensen | same. |
| [mertens-power-saving](discoveries/mertens-power-saving.md) | Mertens | RH ⇒ M(x) ≪ x^{1/2+ε}; our explicit x^0.9382 is a record only. |
| [nyman-beurling-lp](discoveries/nyman-beurling-lp.md) | Nyman–Beurling | RH ⇒ the L^p criterion holds for every 1<p<2; the non-abelian ζ_K version (Riemann/Deep/NBLpK) is not covered. |
| [weil-positivity](discoveries/weil-positivity.md) | Weil | RH ⇒ Weil's form is positive on all test functions. |
| [linnik-numeric](discoveries/linnik-numeric.md) | Dirichlet half-plane | GRH ⇒ Linnik exponent 2 (explicit: Lamzouri–Li–Soundararajan). |
| [littlewood-half-plane](discoveries/littlewood-half-plane.md) | Dirichlet half-plane | GRH ⇒ L(1,χ) ≫ 1/log log q. |
| [primes-between-powers](discoveries/primes-between-powers.md) | Dirichlet half-plane | RH ⇒ primes in (x, x + c√x log x). |
| [primes-short-intervals](discoveries/primes-short-intervals.md) | Dirichlet half-plane | RH ⇒ primes in (x, x + c√x log x) (Carneiro–Milinovich–Soundararajan). |
| [prime-primitive-root](discoveries/prime-primitive-root.md) | Dirichlet half-plane | GRH ⇒ g*(p) ≪ (log p)^c with small c. |
| [class-number-5000](discoveries/class-number-5000.md) | Dirichlet half-plane | census stays true; GRH makes larger ranges available. |
| [class-number-10000](discoveries/class-number-10000.md) | Dirichlet half-plane | same. |
| [unit-group-cover](discoveries/unit-group-cover.md) | Dirichlet half-plane | GRH ⇒ O(log Q) small primes cover every modulus with better exponents. |
| [unit-generated-orders](discoveries/unit-generated-orders.md) | Dirichlet half-plane | classification stays true; its zero-free input is now a theorem. |

Also moot as methods against proving RH: NOGO.md barriers and REFUTATIONS.md (refuted claimed RH proofs); kept as audit record.
Not superseded but scooped: sun-hou (Zhao, Zenodo 23202607, 7 Oct), weil-simple-ground-state (Botha, part), evolution-algebras (concurrent, Liu, Zenodo 23608087).
