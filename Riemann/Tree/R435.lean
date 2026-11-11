/-!
# Round 435 — node N12 (function-field and Selberg-class analogues), round 13, falsification of our root-discriminant bounds and new claims: the Hilbert class fields of Q(√−d) for every class number h ≤ 400 lie above odlyzko31's certified bounds (closest d = 23, 47, 71); the LMFDB API now answers with a reCAPTCHA and the Jones–Roberts tables are offline; Booker–Lee–Zubrilina's LLM-written proofs are stated as Lean-formalized; no testable new N12 claim; Bhattacharjee's class-number counts agree with ours (N11); celestial mechanics, lens design and matching markets transfer no tool

* Node: N12 Function-field and Selberg-class analogues (thirteenth round; R416 before). Round type: falsify our own
  discovery odlyzko31 (root-discriminant bounds from Cureton–Thorner's strip) and new claims.
* [falsify] orchestrator (PARI): the Hilbert class field of K = Q(√−d) is totally complex of degree 2h(K) with
  rd = √d; for every h ≤ 400 the smallest such d (fundamental discriminants up to 4·10^5) gives a field above our
  certified bound for its degree (formal r1 = 0, table Riemann/Deep/Odlyzko31/tables.csv): no violation; closest
  d = 23 (h = 3, n = 6: 4.7958 vs 4.557), 47 (h = 5, n = 10: 6.8557 vs 6.606), 71 (h = 7, n = 14: 8.4261 vs 8.139);
  for h ≥ 64 every rd exceeds the asymptotic limit 23.3878: HOLDS.
* [falsify] LMFDB minimal fields, degrees 2–47: the API now returns a Google reCAPTCHA page for scripted queries (both
  scouts; the agent's 6 Nov check of degree ≤ 10 stands); Jones–Roberts small-discriminant tables: LowDisc 404, NFDB
  500: NOT TESTABLE.
* [check] Booker–Lee–Zubrilina 2611.07921 v1: the omitted proofs of Thms 1.2–1.3 were LLM-written and are stated to be
  Lean-formalized; the spanning proof had LLM assistance (TeX l. 127–139, 148–153); one divisibility check of the proof
  (N, s, q) = (36, 6, 5) recomputed: HOLDS.
* [falsify] new N12 items since 4 Nov (arXiv 2611.03535, 2611.02675, 2611.01745, the October function-field and Ihara
  papers; Zenodo 24111210, 24020618, 23624161, 23564899, 23250932, 23185896): none with a finitely checkable claim
  (2611.03535 Thm 1 is asymptotic in q): NOT TESTABLE.
* N11 side check (orchestrator): Bhattacharjee, Zenodo 24124876 (11 Nov), counts and largest d of imaginary quadratic
  fields for every h ≤ 1000 agree row by row with ours (Riemann/Deep/ImqClass2/perh_h10000.txt): HOLDS.
* Exploration (tool transfers, tested): Hill stability Δ > 2√3 as a gap condition (fails at the first pair,
  0.8887), the Abbe sine condition (no ray invariant; optics tried in R407), Gale–Shapley stable matching (no canonical
  preferences on zeros): no transfer.
-/

namespace Fayed.R435
end Fayed.R435
