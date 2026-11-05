/-!
# Round 308 — node N13 (Hilbert-Polya), third round: Yakaboylu's similarity lemma fails as an operator identity, but his main eigen-theorems are proved directly; only the multiplicity-m extension relies on it

* Node: N13 Hilbert-Polya (third round; follows R302's possible gap).
* Yakaboylu, arXiv 2408.15135 v17: Lemma lemma_similarity (TeX lines 285-309, stated 'formally' on the core C_c^infty
  of L^2([0, inf))) factorises through omega(T)^{-1} and V = omega(T)^{-2} with omega(t) = t e^t/(1 + e^t)^2 (~ t/4 at
  0, ~ t e^{-t} at infinity) and spec T = [0, inf): both are unbounded, and for generic core vectors f/omega is not in
  L^2. A reviewer confirmed the similarity fails as an operator identity on the stated core.
* Dependence: Theorems theo_eigen_R and theo_spectrum_R_adj are proved by direct ODE/Hankel and weak-adjoint
  constructions (lines 327-381, 405-438): independent of the lemma, no flaw. The simplicity/Jordan theorem is a direct
  ODE argument (691-709). The 'Riemann operator for multiplicity m' theorem replaces omega in R_factorized (line 743)
  and transfers spectra through the formal intertwiner: its proof is invalid as written (minor; not an RH claim).
  The RH-relevant step (W > 0, self-adjointness) is conditional and RH-equivalent, as R302 found.
* State of the art: 2607.19067 (Theorem S-eigenstates) shows the Bender-Brody-Muller eigenstates satisfy
  Delta psi_z = x^{-z}, so |Delta psi_z|^2 = 1/x is not integrable on Re z = 1/2: a published obstruction to BBM.
  Connes-Consani-Moscovici 2511.22755 Thm finmain: self-adjoint finite operators with primes cut at lambda^2; the
  agreement with zeta zeros (e.g. 2.44*10^-55 for the first zero at lambda^2 = 13, 2.04*10^-3 for the 50th) is
  observed, not proved. No new Oct/Nov operator claim.
* Exploration: spectral flow/index counting of zeros (2406.01828) only re-expresses N(T) and is RH-equivalent (none);
  Darboux/Crum insertion of eigenvalues works for finitely many levels (toy: Wronskian 128 e^{3x^2/2}(8x^6 - 12x^4 +
  18x^2 + 9) zero-free), but a limit operator needs uniform control of (log W_n)'' and N(T) ~ T log T (N13, idea);
  Furstenberg systems of Liouville give no rate for M(x) (N10, none).
-/

namespace Fayed.R308
end Fayed.R308
