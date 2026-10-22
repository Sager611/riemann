/-!
# Round 121 — quick audits of 6 recent arXiv papers (8 scouts; arXiv PDFs + pdftotext)

* 2610.27830 Visser (20 Oct 2026, primes in short intervals assuming RH), §3.2 proof of Theorem 2: FALSE STEP.
  It says π(n²+n) > π(n²) "is equivalent to" π(x+√x) > π(x). Only the case x = n² holds. The upper-Oppermann
  data (n ≤ 7·10^13) give a prime in (x, x+3√x+2] only, so Theorems 2–3 are unproved for
  1.8·10^19 < x ≤ 4.9·10^27 once x^{1/m} < 3 (m ≳ 40). Verified in the PDF; reviewer pending.
  Theorem 1 is just the crossover with Carneiro–Milinovich–Soundararajan, not a record.
* 2609.26787 Mishra–Sarkar, "finite" ω-form of Robin's inequality ≡ RH:
  - one tail constant is missing a factor e^γ (repairable; the corrected value still fits under 0.095);
  - "finite reduction" is overstated: for each K = ω(n) there are finitely many tuples, but K is unbounded.
  Minor.
* 2509.16240 Kuipers (claimed proof of RH): the contraction inequality is imported from the author's earlier paper
  and assumed. The step Ẽ(X) ≪ X^{1/2} log X is RH-strength, so the argument is circular in effect. No new step.
* 2610.10584 Wunderlich: Eq. (3.1) composes the maps wrongly. With S(z) = i(1+z)/(1−z) and d(w) = 1/2 + iw,
  d(S(z)) = −(1+3z)/(2(1−z)), but the paper writes −(1+z)/(2(1−z)). Later steps use the stated map, so this is a
  normalisation slip.
* 2609.33794 Musin, 2609.02882 Lamzouri (> 67.25% simple zeros on the line, no mollifier; pair correlation +
  Montgomery–Taylor kernel): nothing found at triage depth.
-/

namespace Fayed.R121
end Fayed.R121
