CORRECT: the 182-vertex PG(2, 9) network of arXiv 2610.10108 (Braverman–He) carries all 157 sessions at coding rate 1, while any fractional flow has rate ≤ 147/157; the undirected multiple unicast conjecture is false.

Checked by our own script (verify.py, exact F₉ arithmetic, 21 s, 37 MB), not the authors' code:
1. The encoding is fixed. The paper's Section 3 ("The candidate graph", 3-certificate.tex) sets
   F₉ = F₃[θ]/(θ²+1), η_r = (r mod 3) + ⌊r/3⌋θ, and the order v₀ = (0,0,1), v_{1+r} = (0,1,η_r),
   v_{10+9r+s} = (1,η_r,η_s); Section 4 says it reuses that graph. The linked repository
   (github.com/zhongtianhe/multiple-unicast-counterexample) states the same rule, and its data file
   holds exactly the paper's R, π, C, J. Our round-179 "under-specified" verdict was wrong: the scout
   ranked a + bθ by a first; with that order 150 of 182 vertices break d⁻ = d⁺ = 5, and the two
   other irreducible quadratics break 122–158 vertices.
2. With the paper's encoding: 735 kept edges, 157 sessions, d⁻ = d⁺ = 5 and d_C⁺ = d_J⁻ at all
   182 vertices.
3. The space W of edge labellings obeying all local rules has dimension exactly 166 (rank 744 of
   the 910 × 910 system); the appendix's 166 explicit labellings lie in W and are independent.
4. The causal code, run with symbolic inputs in the order π: every vertex knows exactly 5 values
   when its turn comes (incoming edges, its own messages, public zeros), interpolates, and sends one
   symbol on each outgoing edge. All 157 sinks output exactly their own message; values computed
   twice agree; the final labelling lies in W. Separately, each message lies in the span of what
   its sink holds. Each kept edge carries one F₉ symbol in one direction of an acyclic orientation,
   so with capacity shared between the two directions the coding rate is 1 (no pipelining needed).
5. Routing: every session's endpoints are exactly 5 apart in the kept graph (≥ 5 also follows from
   girth 6). A flow of common rate r uses ≥ 785r of the total capacity 735, so r ≤ 147/157 < 1;
   total throughput ≤ 147 < 157. Still true with all 314 terminals on separate pendant vertices
   joined by unit edges: flow ≤ 1049/1099 < 1, coding still 1.
6. Negative control (verify.py --control, plain quartic rule at line vertices): dimension drops to
   131 and all 157 sinks fail, so the checks can fail. The authors' own verifier also passes.
Not checked: the amplification corollary (cites Braverman–Garg–Schvartzman 2017) and the Section 3
nondeterministic example; neither is needed for the counterexample.
Records: others' result (Oct 2026); round 179's "under-specified" line should read "correct".
