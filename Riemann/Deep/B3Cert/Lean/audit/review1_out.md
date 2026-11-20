**ISSUES / INCOMPLETE — literal numerical equalities only.**

1. **REPORT’s displayed decimals are approximations, not exact constants.**
   - Its κ decimal is excluded by the more precise quoted Arb ball.
   - Its `h(371/750; κ) = 3.3546·10⁻⁵` conflicts with §8’s more precise value `≈ 3.354573606042305·10⁻⁵`.
   - Lean correctly retains the exact symbolic κ₂ and h. It proves `h > 3.35·10⁻⁵`, strengthened by `e1p_enclosures` to `h > 3.3545·10⁻⁵`.
   - **Exact fix:** replace REPORT’s decimal equalities with `≈`; retain the exact κ definition and symbolic lower bound, and state the certified inequality `h(371/750; κ) > 3.3545·10⁻⁵`. Do **not** strengthen Lean to `h ≥ 3.3546·10⁻⁵`.

**Otherwise faithful under the disclosed conditional architecture:**

- **H1P/H1G:** correct strict ranges; their ν bound automatically implies `ν < 1/2`.
- **H2P:** matches the explicitly disclosed MomentRange-based extension, not an undisclosed appeal to SL’s printed range.
- **H3P:** matches the quoted θ-free chain on `(0,1)`; optimizing over Wang’s `λ < θ` is represented by the arbitrarily close bounds in `SL51Inputs`.
- **H4P:** its rational interval contains the entire quoted Arb c-ball.
- **HCW:** its broader conditional domain does not assert Gaussian moments outside their proved range. This matches §5’s claim that the Gaussian moment is the transfer’s only analytic input.
- **η arithmetic:** correct. The chosen smaller exponent is the midpoint of θ and the moment edge, hence strictly above that edge and below θ; `θ < 1` supplies its upper bound. Positive η permits the required fixed D.
- **Conclusion:** correct liminf interpretation, θ range, detector parameters, and inequality directions.
- No contradictory or vacuous premise identified. Abstract consistency witnesses establish consistency, not applicability to actual zero counts; that distinction is already disclosed.