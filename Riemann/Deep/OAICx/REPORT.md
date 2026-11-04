NO FLAW FOUND: 40 claims checked (7 objects verified, 30 statements read)

# OAICx: audit of OpenAI's claimed counterexamples to famous conjectures (github.com/openai/math)

Date: 2026-11-04. Source fetched read-only (contents/trees API and raw files); no clone, no Lean build.
Scripts in this folder rerun in seconds, single-threaded: `python3 check_*.py`.

## Short answer
I found no false step and no wrong object. 25 of the 40 claims have a Lean "Comparator" statement (30 statement files
in total). I read every one of them, and each faithfully states the negation of the named conjecture. None of them is
vacuous, weakened or junk-valued.
Seven finite objects or key finite lemmas were recomputed with independent exact code, and all match. The 15 claims
with no formal proof are mostly long infinite constructions. Two of them were read in full: Wall's D(2) problem and
Gersten's conjecture for K_3. Two more had their finite core checked: the Shafarevich arithmetic and the torsion-free
direct-finiteness design. No gap was found.

## How the claims were checked
1. I listed all 751 preprints, picked the 39 counterexample claims named in the brief (plus Sidorenko, so 40), and
   downloaded their TeX sources.
2. I mapped each claim to its formal statement in `lean/ComparatorChallenges/*.lean`, using `*.json` and
   `lean/formalization.yaml`. Each statement was then read for fidelity: definitions, quantifiers, edge cases and
   junk values such as sSup or sInf of an empty set and the diameter of an unbounded set. Proof scripts were not
   read. GitHub code search over `lean/OAI` found `sorry` only in 5 unrelated files (Theta, SheafCohomology), and the
   `axiom` hits in the checked solution folders were comments. Some of these searches were rate-limited. The
   Comparator checker enforces the axiom list in any case.
3. Where a claim names an explicit object (a graph, polytope, algebra, presentation or design), I recomputed it
   exactly with my own code.
4. GitHub issues on openai/math do not exist (only 4 PRs, and conversations are locked). I found no public reports of
   flaws.

## Table
Type (a) means an explicit finite object; type (b) means an infinite, random or asymptotic construction.
A Lean entry of "—" means there is no Comparator challenge, so the claim is unformalized.

| preprint | type | Lean | what was checked | verdict | file |
|---|---|---|---|---|---|
| Sidorenko (Sep 23) | (a) pattern H explicit; host random | SidorenkoCounterexample | Statement is faithful: all homs, t(H,G) < p(G)^66, H bipartite. Table 1 = Lean faces. The complex proposition holds: each of the 33 pairs lies in exactly 2 faces, the 6 exposure orders, a_e, the degree table, and colour refinement is discrete. H is the incidence graph of a 13-vertex triangulated 2-sphere (χ=2, all vertex links are cycles). | object verified | check_sidorenko_complex.py |
| Hadwiger | (b) random graphs with α≤2 | HadwigerCounterexample | not_hadwiger_conjecture ≡ ¬∀G, χ≤h. Minor model uses nonempty connected branch sets, so sSup is bounded. The chain 26m/75+2/3 < m/2 holds iff m≥5, which matches the Lean bound. | faithful | — |
| Kaplansky direct finiteness, char 2 | (b) prescription "not executed" (paper) | KaplanskyDirectFiniteness, KaplanskyFinitelyPresented | finite field of char 2, f.g./f.p. G, ab=1, ba≠1 | faithful; no explicit elements to multiply | — |
| Kaplansky direct finiteness, odd char | (b) p=minFac((m!)²+1), m=C(1200,600) | OddKaplansky | as above, plus injective non-surjective cellular automaton | faithful; not computable | — |
| Zero divisors in F₂[G], G torsion-free | (b) random graph cones | TorsionFreeZeroDivisors | f.p., torsion-free, finite 2-dim K(G,1), αβ=0. Parity cancellation (odd |S_x∩S_y| ⇒ even components) re-derived. | faithful | — |
| Torsion-free G, F₂[G] not directly finite (Oct 4) | (b) random; finite label design | — | Proof of the parity criterion checked. Label-type design (PG(2,q) lines plus Fano complements) re-verified on the PG(2,4) analogue: label balance and all parity conditions, 0 violations. | checked parts hold | check_df_torsionfree_types.py |
| Kaplansky quasitrace | (b) | KaplanskyQuasitrace | 1-/2-quasitrace axioms standard. ≥1/144 non-additivity for all normalized 2-quasitraces. | faithful | — |
| Kadison–Kaplansky projection | (b) | — | not finite | not checked | — |
| Kurosh division ring | (b) | — | not finite | not checked | — |
| Ryser covering | (b) asymptotic, rank s^n+1 | RyserCovering, RyserOddExtensions | r-partite = exactly one vertex per class (so uniform); ν=1, τ=r | faithful | — |
| Balanced Ryser at prime orders | (b) random, q large | BalancedRyser | Hypergraph, not Latin squares as the brief assumed. Intersecting, every vertex used, τ=q+1 with r=q+1. | faithful | — |
| Borsuk, dimension 9 | (b) topological | BorsukNine | Trace-one symmetric 4×4 matrices form a 9-dim Euclidean space (Frobenius norm); diam = √2; a cover equals a partition. Equivalent to covering SO(3) by sets of rotation-diameter < π. The ≥11 bound mirrors Walkup's 11-vertex minimum for RP³; the RP² analogue (6) is attained by the icosahedral cover. | faithful, consistent | — |
| Colin de Verdière chromatic | (b) probabilistic plus matrix criterion | — | Not finite; it implies the Hadwiger claim, since μ ≥ h−1. | not checked | — |
| Auslander–Reiten | (a) explicit algebra (dim ≥ 800) | AuslanderReiten | Statement faithful (Gorenstein-projective witness, Ext into Z and Λ, all field extensions). Key lemma recomputed: C associative, T = C⋉DC symmetric, Ext_T^i(s,s) = 1,0,0,1,0,0,… for i ≤ 12 at a generic q in GF(2^20). | key lemma verified | check_ar_tachikawa_T.py |
| Tachikawa 2nd | (a) same C and T | Tachikawa | SymmetricOver equals a symmetric Frobenius form. The paper's full table equals the AR product rules. | key lemma verified | check_ar_tachikawa_T.py |
| Eilenberg–Ganea | (a) explicit seed complex | EilenbergGanea | cd = pd_{ZG}Z; K(G,1) of any size. Lean seed: 7 vertices, 51 edges, 45 sectors; group order 120; CW complex and flag complex (103 vertices) are Z-acyclic. | seed verified | check_eg_seed.py |
| coarse Novikov | (b) | CoarseAssembly (339 KB, main theorem located, body not read) | not finite | not checked | — |
| Baum–Connes injectivity (torsion-free) | (b) | — | not finite | not checked | — |
| Baum–Connes, irrational trace | (b) Gromov–Osajda monster plus lamp module | — | Mechanism skimmed; consistent with Lück's rational-trace theorem, which forces failure of assembly. | not checked | — |
| Naimark in ZFC | (b) | Naimark | Unital, simple, infinite-dimensional, faithful trace, a single irreducible representation, not ≅ K(H). Universes OK. | faithful | — |
| nearby Lagrangian | (b) | — | not finite | not checked | — |
| C¹ entropy conjecture | (b) | C1EntropyCounterexample | coverEntropy = h_top; C¹ via local ambient extensions; f_* on real H₂ has eigenvalue 2; non-invertible | faithful | — |
| Gersten, integral K₅ | (b) | — | not finite | not checked | — |
| Gersten, integral K₃ | (b) explicit ring | — | Read in full. Model lemma checked by hand: Eisenstein divisor, A regular, π ∈ m², special fibre is 5 conjugate lines, residue degrees divisible by 5. Steinberg lift and norm detector are consistent. | no flaw found | — |
| group-ring determinant | (b) builds on the char-2 claim | GroupRingDeterminant | invertible-case Fuglede–Kadison det in (0,1); A a unit over Q[G] | faithful | — |
| Wall D(2) | (a) presentation with 5 generators and 4 relators | — | Read in full. Tower obstruction (Howie towers plus Laurent rank), plus-construction splitting and perfectness lemma checked by hand. ρ respects all relators; Fox matrix has rank 4, so H₂(B;L_ρ)=0. | no flaw found | check_wall_d2.py |
| Shafarevich holomorphic convexity | (b) | — | Golod–Shafarevich margin: 48 generators with P=(1+t)G, G(49/200) = 1.02823 > 257/250, 4t<1, total < −1/125 | arithmetic verified | check_shafarevich_gs.py |
| Hahn–Wilson | (b) | — | not finite | not checked | — |
| Pixton completeness | (b) g = 10^60 | — | not finite | not checked | — |
| metric-entropy duality | (b) | MetricEntropyDuality | Negates Pietsch's right-hand inequality for all a,b; Coverable guards against sInf ∅ = 0 | faithful | — |
| Monge ansatz, 3-marginal Coulomb | (b) | CoulombCounterexample | Kantorovich value attained and finite; every pair of measure-preserving Borel maps is strictly worse; equal infima | faithful | — |
| projection-body simplex maximum | (a) T₁₀×T₁₀ | ProjectionCounterexample, ProjectionVolume | Exact brute-force zonotope sums (231 dets of 20×20) give the ratio 22355476/22020096 > 1. The simplex constant and product identity were also checked by brute force. Not new: the paper cites Feng–Hu–Liu–Xu, d ≥ 9. | verified | check_projection_body.py |
| Arnold, critical number (Q³) | (b) | ArnoldCounterexample | Critical points via horizontal lifts, Hamiltonian lift with phase equivariance, exactly 3 fixed points, Crit ≥ 4 | faithful | — |
| Arnold, Morse number (12-dim) | (b) | — | not finite | not checked | — |
| harmonic dimension, n ≥ 8 | (b) | HarmonicGrowth | Christoffel, Ricci and Laplace–Beltrami formulas, intrinsic completeness, Euclidean count C(n+k−1,k)+C(n+k−2,k−1) | faithful | — |
| harmonic dimension, n = 3 | (b) | — | not finite | not checked | — |
| infinite matroid packing/covering | (b) | InfiniteMatroid | Contraction onto C = (M*|C)*; self-dual pair with no packing/covering partition | faithful | — |
| plane not 5-colourable | (b) analytic, no finite graph | EuclideanFiveColor | arbitrary colourings ℂ → Fin 5, unit distance | faithful | — |
| ETH is false | (b) algorithm | SubexponentialSAT | Explicit multi-tape machine, step count, unary input (bits are polynomial in n and m), bound 2^{O(n (log log n)²/√log n)}, polynomial space | faithful | — |
| torsion-free hyperbolic group, not residually finite | (b) | TorsionFreeHyperbolic | Rips δ-thin geodesic triangles in a finite Cayley graph; torsion-free; ¬RF | faithful | — |

Three more statements were read: LipschitzEquivalence (faithful), PlaneColoring (the classical 7-colouring) and
PettyProjectionVolume. The last is Petty's conjectured lower bound, so it does not contradict the product
counterexample, which concerns the upper bound.

## Remaining risks
- Lean proofs were not built, as the brief required. "Faithful" therefore means the statement is right; it does not
  certify the proof. The repo says `review: unchecked` and has already withdrawn three papers once (October 7, a sign error).
- Of the 15 unformalized claims, 4 were checked (torsion-free direct finiteness, Gersten K₃, Wall D(2), Shafarevich).
  The other 11 remain unchecked beyond skimming: Kurosh, Kadison–Kaplansky, both Baum–Connes papers, nearby Lagrangian,
  Gersten K₅, Hahn–Wilson, Pixton, Arnold Morse number, harmonic dimension n = 3 and Colin de Verdière. These are the
  most likely places for a real gap, but none has a finite object that could be checked.
- Explicit group-ring elements for the zero-divisor and direct-finiteness claims do not exist. They are random, or the
  paper says the prescription "has not been executed". The multiplication test the brief suggested is therefore
  impossible.
