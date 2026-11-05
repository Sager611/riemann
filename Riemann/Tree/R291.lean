/-!
# Round 291 — node N7 (zero-free half-planes), second round: Hexagon endpoints and twist step check out; method limit 13/15

* Node: N7 zero-free half-planes (second round; Orca hexaudit is auditing the Hexagon paper in full).
* Endpoint arithmetic (Liu's Hexagon TeX, number_field_new_v1-4.tex, Eqs. tw:eq:Bpoly etc.): B_tw = 0.874058503791859005,
  B_ex = 0.874072847666827908, B' = 0.874221250492902083 are the unique roots in their brackets (exact rational interval
  checks, 80-digit residuals < 1e-73). The twist step's bound alpha_max = (4p - 6)/(p(p + 24)) on p in [112/25, 9/2] checks
  (grid + interval + exact). The arXiv identifier of the Hexagon paper is unclear (a scout found 2610.00655v1 to be an
  unrelated paper; abstracts it found list only B' and B_ex): to be settled by hexaudit.
* Method limit (Liu's own section on proportional widths): B -> 13/15 = 0.8667 as the return width C -> 1/4; reaching
  B = 0.87 needs C > 0.2171 (now C_tw = 0.0952); so this method family cannot pass 13/15.
* Formal status: OpenAI's Lean proof covers zeta and Dirichlet L (7/8) and finite-order Hecke L over Q(sqrt -3) (only
  propext, Quot.sound, Classical.choice); Liu's general-field inputs are not formalised.
* Exploration: Alpoge-Furman (arXiv 2608.13637) machine-verified > 2/3 of zeros simple and on the line (TREE.md N8 updated);
  suggested formalisation target: the kappa >= 13/18 input of Liu 2610.26835 (flagged in R131), est. 3-5 weeks; cubic-theta
  Whittaker coefficients as local positivity detectors (N7); AI-produced additive results (2610.19470, 2609.39314) as audit
  candidates (off-tree).
-/

namespace Fayed.R291
end Fayed.R291
