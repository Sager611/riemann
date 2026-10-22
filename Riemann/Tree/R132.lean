/-!
# Round 132 — independent checks of the new A > 1 result (RCApps5): reviewer + 8 scouts

Claim: for every A > 1, almost all fundamental D (power saving) have primes in every (narrow) class above
|D|^A. This is the GRH exponent for this method.
* Reviewer: the theorems as proved are correct.
  - The family large sieve holds for the detector sequences only.
  - It is FALSE for arbitrary coefficients: b = 1 on primes makes the trivial characters give ≫ Xπ(N)².
    Never advertise it as an arbitrary-coefficient large sieve.
  - Lemma P is valid: blocks are ring class groups with unit factors.
  - The Poisson steps are uniform, and the density and Markov deductions hold.
* Scouts (8/8 ran):
  - checks/lemmaP_check.py: 388 tests pass, plus 257 for D > 0 and 1185 for sectors;
  - poisson_check.py: relative deviation 1.3e−4;
  - §4 exponent arithmetic OK (minimum gap 3.0e−5);
  - §5 Markov deduction OK on a grid A ∈ (1, 1.2] (minimum margin 5.5e−6);
  - §3 application conditions OK; D > 0 narrow units OK; sector parameters and the A = 1 barrier OK.
  - One wording fix, NOTE.md l. 74–76: "an ideal of norm M prime to f_e is (γ)" must read "a PRINCIPAL ideal"
    (the identity-class part, r_{Δ_e}(M) = (w_e/2)·[1](T_M)). It does not affect the lemma.
* Next: adversarial novelty (Orca rcnov5); transfer to CM fields (Orca cm3); possibly lower the
  Frobenius-sign / Sato–Tate thresholds with the same fix.
-/

namespace Fayed.R132
end Fayed.R132
