/-!
# Round 283 — node N4 (Robin/Lagarias): record re-checked (10^{10^{13.11485}}, 21-free); margins computed; no proof claims to refute

* Node: N4 Robin / Lagarias.
* State of the art: Morrill-Platt (arXiv 2110.13478, Thm 1 and Lemma 3) verify Robin for all 5041 <= n <= 10^{10^{13.11485}}
  and prove it unconditionally for all 21-free n; no later extension found (2503.03159 reuses it). TREE.md N4 updated.
* Computation: Robin margins positive for colossally abundant n with log n up to 1.23*10^7 (margin * sqrt(log n) about 9.0*10^4
  there); Lagarias' form: smallest relative margin for n <= 10^7 at n = 12 (0.0114), and among superabundant n beyond, the
  smallest found 0.0293 at n = 2,021,649,740,510,400 (2602.15905 reduces a Lagarias counterexample to superabundant n).
* Adversarial: Zenodo 20931450 proves only restricted families (no all-n claim); nothing to refute. One scout (21-free lemma
  check) timed out.
* Exploration ideas: LLL search for near-relations sum a_p log p ~ 0 with Robin weights to enumerate near-violators (N4);
  entropy-decrement bounds for Liouville correlations (2610.12869) toward explicit psi(x) - x errors (N10, bridge unproved).
-/

namespace Fayed.R283
end Fayed.R283
