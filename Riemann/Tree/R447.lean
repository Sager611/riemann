/-!
# Round 447 — node N8 (critical-line zeros, large gaps), round 17, check of the inputs and escape routes of our record Λ > 3.324 (gapbar: Bui–Milinovich's criterion with the amplifier at the midpoint, M(s + iκπ/(2ℒ))): both reductions and the derivative step survive the shifted amplifier; no twisted fourth moment beyond ϑ = 1/4; the midpoint shift also applies to BM's h₂ criterion, but BM's own deduction with the 2026 simple-zero proportion already gives Λ* ≥ 6.1263 for gaps between multiple zeros (was 4.05; routine); mixtures of different r decouple, real r untried; Bredberg and Bui–Hall–Subirà Jorge offer no better route; flag algebras, quantum ergodicity and van der Corput sets transfer no tool

* Node: N8 (17th round; R444 before). Round type: [check] the new inputs and the escape routes of a fresh record of ours
  (gapbar, 12 Nov 2026, Riemann/Deep/GapBar: certified Λ > 3.324 with F = e^{ivtℒ}ζ(s)ζ(s+iκπ/ℒ)M(s+iκπ/(2ℒ)), r = 2,
  degree-80 real P, ϑ = 2499/10000; float ceiling 3.32658). BM's Lemma main with shifts α₃ = −β₃ (l.364-384) was checked
  in the TeX by the orchestrator and by gapbar's second reviewer before this round.
* [check] Reductions with the shifted amplifier (BM 1410.3635v3 l.137-233, l.385-410): complex Wirtinger needs only endpoint
  zeros; the zeros of F still contain t_n and t_n − κπ/ℒ (l.182: consecutive zeros of F differ by ≤ (1 + o(1))κπ/ℒ under
  Λ ≤ κ; zeros of M only add points), multiple zeros of ζ stay multiple zeros of F (Theorem gap2, l.212-224), |F| is free
  of v, and 'mimic a real-valued function' (l.174) is descriptive and unused: HOLDS.
* [check] Twisted fourth moment: no 2019-2026 extension of Bettin–Bui–Li–Radziwiłł's ϑ < 1/4 (general, mollifier or divisor
  coefficients, conditional or not), no product A(s)B(s+iα) version; Bui–Hall–Subirà Jorge 2511.14415 use ϑ < 1/8: HOLDS.
* [escape] BM's second criterion h₂ = 4c₁/(κ²c₂) (l.222, l.329; Λ* > 4.05 with r = 1, ϑ = 1/4, P = 1 − 5.2x + 5.5x², l.341-345)
  uses the same smoothed c_j (Theorem smoothed l.298-313), so the midpoint shift applies. But Λ* (gaps between consecutive
  multiple zeros on the line, l.73-78) obeys Λ* ≥ 2/(1 − C) whenever N^s(T) ≥ (C + o(1))N(T) (BM l.78-81), and
  2610.08965 Thm thm:line gives C = s₀ = 1669159/2478195: Λ* ≥ 2478195/404518 ≈ 6.1263 (routine, COROLLARIES), so the
  h₂ route is moot. No later Λ* bound in arXiv titles 2026-10/11, Zenodo or BM's 22 OpenAlex citers.
* [escape] Amplifier class: BM fix d_r(h)P with r ∈ ℕ (l.244-246); no mixed-divisor amplifier or cross-term main term in print
  (Hall, Bredberg, BMN, Feng–Wu, BHSJ, BBLR). Mixtures of different r give no gain: the d_{r₁} × d_{r₂} cross terms have order
  (log T)^{(r₁+2)(r₂+2)} against the blocks' (log T)^{((r₁+2)² + (r₂+2)²)/2}, so the normalised blocks decouple and h₁ is a ratio
  of sums, at most the best block. Real r (Γ in place of the factorials; BM's lemma l.441 is stated for integers) is untried; a
  quadratic through the float ceilings 3.29751, 3.32658, 3.23819 (r = 1, 2, 3) peaks near r ≈ 1.75 at ≈ 3.330 (crude).
* [escape] Bredberg 1101.3197 uses the same product e^{viθ}M(½+it)ζ(½+it)ζ(½+it+iκ/log T) (l.60-64; Λ > 2.766), i.e. BM's
  family; Bui–Hall–Subirà Jorge 2511.14415 prove Λ > 2.64 with f = ζ(½+it)A(½+it) and only propose
  Z(t)Z(t+2πκ/log T)A(t)A(t+2πκ/log T) (l.492-500), without numbers, at ϑ < 1/8: no better route.
* Exploration (tool transfers, tested): flag algebras (Goodman's triangle inequality needs graph limits; finite zero samples give
  no statement about ζ), quantitative quantum ergodicity (needs a self-adjoint quantisation and observables, absent for zeros),
  van der Corput's inequality on Σ e^{itγ_n} (holds trivially, 3.0835 ≤ 6.4850 at N = 3, no RH information): no transfer.
-/

namespace Fayed.R447
end Fayed.R447
