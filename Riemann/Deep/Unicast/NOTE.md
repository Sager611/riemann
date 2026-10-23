# Audit notes: arXiv 2610.10108, "The Multiple Unicast Conjecture is False" (Braverman, He)

## Claim
Undirected network: the point–line incidence graph H of PG(2, 9) (182 vertices, 910 edges, girth 6)
minus 175 listed edges R. Kept graph G: 735 edges. Sessions: B = R \ (C ∪ J), 157 deleted edges;
source = endpoint earlier in a listed vertex order π, sink = the later one. Coding rate ≥ 1;
routing ≤ 735/(157·5) = 147/157.

## Conventions (where the earlier scout went wrong)
- Section 4 does not restate the encoding but says "the same ... graph as in the previous section".
  Section 3 fixes it: F₉ = F₃[θ]/(θ²+1), field order η_r = (r mod 3) + ⌊r/3⌋θ (labels 0..8 = a + 3b),
  points/lines v₀ = (0,0,1), v_{1+r} = (0,1,η_r), v_{10+9r+s} = (1,η_r,η_s); P_i → label i,
  L_j → 91 + j; edge e_{10i+r} = {i, 91 + j_{i,r}} with the line indices increasing.
- Repository README.md and CONSTRUCTION.md state the same; data/182-vertex-code.json = paper lists
  (order = π, backward = R, C, J identical).
- The local rule does not depend on the choice of basis u, w of z^⊥ (linear changes of (λ, μ) keep
  the space of quartics; at lines, cubing is additive in characteristic 3). Changing a vertex's
  representative only flips signs of its edges, an equivalence. So only the labelling matters.

Degree check d⁻ = d⁺ = 5 under each tried encoding (vertices failing, of 182):

| modulus     | field order η_r = (r mod 3) + ⌊r/3⌋θ (paper) | η_r = ⌊r/3⌋ + (r mod 3)θ |
|-------------|----------------------------------------------|---------------------------|
| θ² + 1      | 0                                            | 150                       |
| θ² + θ + 2  | 122                                          | 149                       |
| θ² + 2θ + 2 | 126                                          | 158                       |

## Results under the paper's encoding (python3 verify.py)
- d⁻ = d⁺ = 5 and d_C⁺ = d_J⁻ at all vertices; |E(G)| = 735; |B| = 157.
- dim W = 166 (exact rank over F₉); the 166 monomials X^aY^b (not divisible by X₀Y₀³ or X₁³Y₁)
  give labellings h(p⁽³⁾, ℓ) in W, independent: the appendix bound holds with equality.
- Causal code with symbolic inputs: exactly 5 known values at each vertex, 735 transmissions, one per
  kept edge, all earlier → later; decoded value at every sink = its own message (identity on F₉¹⁵⁷);
  J-edge values agree at both ends; C-edges stay 0; the global labelling is in W; each message is in
  the span of its sink's incoming symbols and own messages (decoder-independent check).
- Distances: 5 for all 157 pairs, sum 785.

## Why this refutes the conjecture
Model (Li–Li 2004; Harvey–Kleinberg–Lehman 2004; paper §2): each undirected edge has capacity 1 split
between its two directions; coding may be arbitrary. Here one symbol per edge in one direction of an
acyclic orientation, one symbol per message: rate 1. Flow: each unit of session q needs ≥ 5 units of
edge capacity, total 735, so common rate ≤ 147/157 and total throughput ≤ 147 < 157. With pendant
terminals (unit edges, all 314 distinct) the gap persists: flow ≤ 1049/1099.

## Not checked
Braverman–Garg–Schvartzman amplification (Ω((log n)^ε) gap corollary); Section 3 nondeterministic
example; the concurrent OpenAI Fourier-circuit route mentioned in the introduction.

## Files
- verify.py: full check (--alt: encodings table; --control: wrong line rule, everything fails).
- instance.json: R, π, C, J copied from the e-print by parse_eprint.py.
- Cost: ~21 s, 37 MB peak, single thread.
