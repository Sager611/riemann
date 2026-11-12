VERDICT: ISSUES

1. **BLOCKING — Scope:** `ceiling_poly` requires reflection symmetry and `ceiling_bernard` covers Bernard’s special complex family, neither covering arbitrary complex polynomials satisfying only `Q(0)=1` in the supplied Andersen–Thorner formulation; restrict the advertised ceiling to Bernard’s class or establish the unrestricted result.

2. **BLOCKING — Certificate discharge:** `ceiling_bernard_json` remains conditional on coverage, certificates, and three tails, while `CrudeCellCert` and the tails directly assume the required regional lower bounds rather than state the described explicit arithmetic tests; supply proved coverage, verified dual witnesses, and concrete crude/tail inequalities implying these hypotheses before claiming an unconditional ceiling.

3. **MINOR — Closed-form endpoint:** The checker’s formulas agree with `Cl` and `Bl` for positive λ but are undefined at λ=0, which `Cover` includes and Lean handles through `Pl`/`Pld`; explicitly bridge the integral and closed-form definitions and certify the endpoint values `Cl 0 = 1`, `Bl 0 = 1/3`, including interval evaluation touching zero.

4. **MINOR — Source interpretation:** `cC` matches Andersen–Thorner’s modulus square, not Bernard’s literal displayed square for complex coefficients, whereas `cBernard_eq` addresses only real-valued Q; explicitly identify the complex definition as the modulus-squared interpretation supported by Andersen–Thorner rather than a verbatim formalization of Bernard’s display.