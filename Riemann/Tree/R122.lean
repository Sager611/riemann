/-!
# Round 122 — quick audits of 7 brand-new arXiv papers (8 scouts, run 4 at a time: at 8 they did not start, machine load 80)

* 2610.21099 Kerr, simpler proof of NO Landau–Siegel zeros (ChatGPT-assisted; follows OpenAI 1 Oct 2026):
  - one slip: "K = Q(√q, √2) has degree four" is false for q = 8. Harmless, since it is a single modulus.
  - no other false step at triage depth. The claim is huge. → Orca siegel: line-by-line audit.
* 2610.23467 Cobeli–Zaharescu, > 8/9 of zeros simple or on the line. The logic is fine (union counted correctly),
  but everything hangs on an uncertified numerical minimum: 0.04400000929 vs the needed 0.044 (margin 9.3e−9).
  → Orca cz89: certify or refute with Arb.
* 2610.24612 Arias de Reyna, exposition of Claude's 0.6725 theorem (simple + on the line): it says the AI source
  "doesn't work" as written and gives a repaired proof.
* 2610.21511 Kimura: 16,453,493 imaginary quadratic fields with h ≤ 2000. Completeness rests on an explicit
  no-Siegel-zero bound, (1 − β) log q ≥ 5.811·10^{-4}, plus sieve logs that are not certified.
* 2610.22939 Punch, explicit Burgess bound for composite moduli: the final constant 4.07325 rechecked.
* 2610.26841 Lamzouri, isolated zeros:
  - the scout's "false functional equation" was a pdftotext artifact (it drops overlines; the TeX has
    Λ(1 − s, χ̄));
  - e^{iθ}Λ(1/2 + it, χ) is real for complex χ too, so there is no flaw.
  Rule added: check notation in the TeX source.
-/

namespace Fayed.R122
end Fayed.R122
