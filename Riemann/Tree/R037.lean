/-!
# Round 37 — rigidity and uniqueness: which axioms pin down ζ?

| theorem | says | filter consequence | verdict |
|---------|------|--------------------|---------|
| Hamburger 1921 | Dirichlet series + ζ's FE (conductor 1) + growth ⇒ c·ζ | refines the DH filter: DH kills arguments that also run for conductor-q FEs | CIRCULAR (pins ζ, not its zeros) |
| Kaczorowski–Perelli | Selberg class of degree < 2 classified | impostors excluded only by the Euler product | DEAD |
| Weil / Bochner converse theorems | twisted FEs characterise automorphic L | identifies the function, not zero locations | DEAD |
| Bombieri–Hejhal 1995 | DH-type combinations: 100% of zeros on the line | density one ≠ all | CIRCULAR |
| moving a zero pair of ξ | FE and order survive (cf. Fayed.R11.fe_allows_any_zero) | the Dirichlet coefficients break first | DEAD |
| Katz–Sarnak families | symmetry types | weaker than Euler-product rigidity | DEAD |
| Euler product alone (Beurling, DMV 2006) | counter-models exist | need Euler product + FE + Selberg axioms: open | CIRCULAR |
| Saias–Weingartner, Conrey–Ghosh | combinations have zeros in σ > 1 | multiplicativity necessary | DEAD |

Sharpened filter: an argument must use data that (a) pins ζ (Hamburger) and (b) is not shared by
conductor-q impostors; Fayed's local data fail already by σ₀-blindness (Fayed.R1.fayed_profile_realizable).
-/
