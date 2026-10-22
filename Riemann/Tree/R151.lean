/-!
# Round 151 — numbering correction; a counterexample in Sun's 1811.10503; reviews; OEIS sweeps

* Numbering correction. arXiv 1309.1679 numbers its conjectures 3.1–3.16 in print (= Nanjing Univ. J. Math.
  Biquarterly 36 (2019) 134–155); the TeX labels differ, and R146–R150 used the labels. Printed numbers of our
  counterexamples: 3.2 (label 3.3), 3.3(ii) (3.4), 3.4(i) (3.5), 3.8 (3.10; arXiv 2610.06289 cites it correctly, so
  R146's "numbering slip" is withdrawn), 3.10(ii) (3.12). Sun–Hou is printed 3.5(i). Fixed in Riemann/Deep/SunPrim
  and in TREE.md / NOGO.md.
* Sun, arXiv 1811.10503, Conjecture 4.9 (|S(A)| ≥ min{p(F), (|A| − n)n(n+1)/2 + n(n² − 1)/6 + 1} for
  S(A) = {Σ k a_k : a_k distinct}): fails for n = 2 and every 3-term progression, e.g. A = {0,1,2}: S(A) = {1,2,4,5},
  4 < 5, over Q and every F_p with p ≥ 5 (Riemann/Deep/SunPrim/verify5.py). Its Conj. 4.1, 4.10 survive Z/m, m ≤ 10.
* 1309.1679 v1 → v9: v1 Conj. 2.2 dropped; 2.4 → 3.3 hypotheses changed ("2 ∤ n or n ∤ |G|" → "|G| odd"; (ii): n < |G| →
  n ∤ |G|; both versions fail for {0,1,2} ⊂ Z/5); 2.5 → 3.4 rewritten, v9 records the (Z/3)^r obstruction itself.
* OEIS sweeps: nothing new ("In general" p. 50–60, "real numbers": A185645(i) survives [−8, 8] and random rational
  sets; "finite additive"); two sweeps timed out. Sun's 1402.6641 and 1211.1588: no failure (tested ranges in NOGO).
* Reviews: STLemP re-review CORRECT (repair covers every AFE block; thresholds unchanged). AllDisc: GAP in
  Proposition 1° (NOTE l. 104–127: the T-dependence is not linear for σ < 7/8; restrict to σ ≥ 7/8, which is all the
  theorem uses) and REPORT line 1 must say "positive definite or indefinite". → Orca fix.
-/

namespace Fayed.R151
end Fayed.R151
