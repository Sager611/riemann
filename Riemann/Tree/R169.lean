/-!
# Round 169 — a Zenodo "proof of RH" refuted (wrong-sign positivity certificate); one open question answered in part

* Refuted: Zenodo 10.5281/zenodo.23555113 (Oct 2026), a claimed proof of RH via Weil positivity. Its certificate
  (PDF p. 8) needs c_a = 2^(1−a) √π Γ((1−a)/2) / Γ(a/2) > 0 at a = 1.4503086547 and states c_a = 2.34; in fact
  Γ((1−a)/2) = Γ(−0.2252) = −5.3018 < 0 and c_a = −5.4587453140 (two scouts and our mpmath check, 40 digits; the PDF
  uses |Γ(−0.225)|). So the kernel K_a is not positive-definite and the step W(Φ) ≥ 0 fails.
* 2610.10095 open question (which orders of +, ×, −, ÷ give perfect sequences for all large N; TeX l. 921–925): exact search
  N ≤ 35 for all 24 orders. Ten orders never work, by the first two moves: the six starting with −, the two starting
  with + −, and × − + ÷, × − ÷ +. Two more work only for N = 2 (÷ − × +, ÷ − + ×). Eleven are perfect for every tested
  N from 8, 9, 10 or 12 up to 35; ÷ × − + is perfect only for N = 2, 4, 5. A full answer is with an Orca agent.
* Checked, nothing new: Zenodo 23462601 (Weil-kernel reductions; the decisive estimate is left open) and 23348776
  (Speiser–Lehmer identity correct; its certificate is conditional); Chen's LIS log-concavity (2610.18235) verified by us
  for n ≤ 30 (Bóna–Lackner–Sagan: n ≤ 50); a series for Catalan's constant conjectured in 2609.04176 (TeX l. 152–155)
  agrees to 10⁻⁹.
* Empty: arXiv NT and CO conjecture finders (API filters); double-comet search 21 ≤ n ≤ 60 timed out.
-/

namespace Fayed.R169
end Fayed.R169
