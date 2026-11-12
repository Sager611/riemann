ISSUES

- **(a) Parameter restriction:** B covers every Lemma T kernel only when **λ ≤ 31/40**, not its full λ ≤ 1 range. A does not fill the remaining interval.
- **Normalization/sign:** Set \(c=\widehat f(0)>0\), replace \(f\) by \(f/c\). All kernel hypotheses persist; \(J_\theta(f/c)>1\). The closed-set sign condition implies Lean’s open-set condition—correct direction.
- **(b) Complex encoding:** Acceptable. A real Schwartz function lifts to a complex Schwartz function with real values. Real even integrable functions have real Fourier transforms, so the real-part inequality is faithful.
- **(c) Generality:** `toDecay` proves inclusion. Strict inclusion holds for λ>0: appropriately scaled sinc² has tent-shaped Fourier transform and satisfies `AdmissibleDecay`, but lacks compact support.
- **(d) Logical direction:** \(J≤1\) contradicts normalized Lemma T premises in the barrier region. This establishes a barrier, **not** Lemma T’s construction of \(F\) or its conclusion \(D(λ)>0\).
- **Certificate qualification:** B remains conditional on `CertB`. Its inhabitance is not evidenced here; compilation alone does not establish an unconditional, nonvacuous certificate-based barrier.
