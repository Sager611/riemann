# Discoveries

- Fayed's RH proof fails only at (55): its Appendix Theorem 1 ("L'Hôpital on a disconnected domain") is false (Lean counterexample on the paper's own D_γ), and every other step is machine-checked from 4 standard analytic inputs (Euler–Maclaurin, truncated Perron, Leibniz rule, convexity bound); `fayed_gap` is axiom-clean. Reviewer: significant. Online: not found. [details](Riemann/Fayed/Theorem1.lean), [proof check](Riemann/Fayed/HARDEN_REPORT.md)
