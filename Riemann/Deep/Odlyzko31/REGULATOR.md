NO RECORD: the 31/32 strip, i.e. kernels of sech(15x/32) type, cannot improve the unconditional regulator lower bounds of Zimmert and Friedman. Zimmert proved R/w ≥ 0.04·e^{0.46 r1 + 0.1 r2}; Friedman proved R/w ≥ 0.0031·e^{0.241 n + 0.497 r1}. Their proofs never use the zeros of ζ_K: they run on the partial zeta function of the principal ideal class and get their positivity from nonnegative Dirichlet coefficients and nonnegative Mellin kernels. That function has no Euler product, and it can vanish even in Re s > 1, so a zero-free strip has nowhere to act. No new constant, so nothing to certify.

# odlyzko31 follow-up: regulator lower bounds (6 Nov 2026, 14:19–14:40 UTC)

## Question
Do the explicit-formula regulator bounds use zero positivity, so that the strip |Re ρ − 1/2| ≤ 15/32 could raise the exponents or constants, as it raised the discriminant constants 22.38/60.84 → 23.39/65.64 in [REPORT.md](REPORT.md)? The three candidates are Zimmert (1981), Friedman (1989) and Friedman–Skoruppa (1999).

## Sources read
We read OCR text of the Göttingen digital library (GDZ) scans: Zimmert, Invent. Math. 62 (1981) 367–380, and Friedman, Invent. Math. 98 (1989) 599–622. Both are on GDZ with OCR files, pages `gdzocr/PPN356556735_0062/00000373–386` and `gdzocr/PPN356556735_0098/00000605–628`. The formulas are partly garbled by the OCR; the statements below are taken from the legible text.

**Zimmert**
- **Satz 1, the fundamental inequality.** It applies to Dirichlet series f, g with **nonnegative coefficients**, meromorphic with a simple pole at s = 1 and equal residues κ, linked by a functional equation with factors Γ(s/2)^a Γ((s+1)/2)^b.
  - **Lemma 1:** (1/2πi)∫_{(2)} x^s Q(s) ds ≥ 0 for every x > 0, where Q = R(s)·(a Γ-quotient) and R is a rational kernel.
  - **Proof:** the integral of x^s Q(s) f(s) is ≥ 0 coefficient by coefficient. Moving the line to Re s = −γ picks up the pole at 1. The reflected integral is bounded through |g(1 + γ + it)| ≤ g(1 + γ).
- **Satz 3, regulator.** It applies Satz 1 with A = √d π^{−n/2}, a = r1 + r2, b = r2, and **g the partial zeta function of the principal class**. Its residue is κ = 2^{r1}(2π)^{r2}R/(w√d), and the proof uses only "−(g'/g)(1 + γ) ≥ 0 and g(1 + γ) ≥ 1" (the unit ideal).
- **Korollar (i):** "R/w ≥ 0,04 exp(0,46 r1 + 0,1 r2)", taking γ = 1. Korollar (ii): R > 0.056 for every field.
- **Satz 4, Korollar 1 and 2** (an ideal of small norm in each class C or 𝔇C^{−1}, constants 60.8 and 22.3) use the same scheme with partial zeta functions f, g. Here −f'/f(1 + γ) ≥ log N𝔞₀ is a weighted average of log N𝔞, not an Euler-product log derivative.

**Friedman**
- **Theorem A / Thm 2.1:** R/w = (1/2πi)∫_{(δ)} P(s)Λ_𝔄(s) ds, with P(s) = −P(1 − s) and 𝔄 the principal class (together with the different's class).
  - The text proves it by shifting the contour to Re s = 1/2 and using the functional equation of the partial zeta function Λ_𝔄. It says formula (1.1) "lies essentially no deeper than the functional equation".
- **§4, Zimmert's method.** Lemma 4.1 is the same Mellin positivity. The Corollary on p. 620 is "R/w ≥ 0.0031 exp(0.241n + 0.497r1) ... Proof: set γ = 0.6" in Zimmert's Satz 3.
- **§2, Proposition 2.2.** For small |D| the unit ideal gives the main term. The other terms are positive once all ideals in 𝔄 have norm above √|D|/x₀, which "Odlyzko's estimates prove" for small |D|. The inputs are the discriminant floor d₀(r1, r2) and the Hilbert-class-field threshold d₁ = d₀(3r1, 3r2)^{1/3}: if h ≥ 3, the Hilbert class field has degree 3n and the same root discriminant.

**Friedman–Skoruppa (Invent. Math. 135, 1999)**
We found no open copy: GDZ ends before volume 135, and the paper is not on arXiv. Not checked. Secondary sources (Akhtari–Vaaler, arXiv 1508.01969 l. 434) only cite their relative-regulator inequalities.

## Why the strip cannot act
1. **R is a residue of a partial zeta function, not of ζ_K.** R/w is the residue of the principal-class partial zeta function ζ₀(s) = Σ_{(α)} N(α)^{−s}, which isolates R from h. The full ζ_K has residue proportional to hR/w, so its zeros, and any strip for them, bound hR, not R.
2. **ζ₀ has no Euler product.** There is no explicit formula over zeros of ζ₀, so no term Re Φ(ρ) for a positive-definite kernel to control. In the Zimmert–Friedman proofs positivity comes from the coefficients of ζ₀ (the unit ideal plus nonnegative terms) and from Mellin kernels (Lemma 1 / Lemma 4.1), never from zeros.
3. **No strip exists for ζ₀ in general.** ζ₀ = (1/h)Σ_χ L(s, χ) over class-group characters. For imaginary quadratic fields with h > 1 it is an Epstein zeta function, which has infinitely many zeros with Re s > 1 (Davenport–Heilbronn, 1936). So no zero-free strip can hold for ζ₀, and a sech(wx) kernel has nothing to replace.
4. **Contrast with discriminants.** Zimmert's Satz 4 with f = g = ζ_K, which has an Euler product, reproduces Odlyzko's unconditional constants 60.8/22.3 (his introduction; arXiv 2412.05568 l. 570–597: F_{a,b}(β) → a(γ + log 4 + 1) + b(γ + log 4 − 1)). There the zero-positivity route exists, and the strip improves it ([REPORT.md](REPORT.md)). For a single ideal class no such route exists.

## The only indirect channel
Friedman's §2 tables use unconditional discriminant floors d₀ and d₁. Our strip bounds raise them, for example by about 0.2–2% in the root discriminant for degrees 8–100. That only shifts the discriminant ranges inside his small-signature tables. Those served Theorem B, the complete list of fields with R/w < 1/8, which is already final.

The asymptotic exponents come from Zimmert's Satz 3 for large |D|, where no discriminant input appears. CM fields, whose units come from the totally real subfield, keep Zimmert as the only tool. So no general constant changes.

## Novelty scan (first)
- `.scratch/lst.py 'regulator|units|Zimmert|Friedman' math.NT` over October and November 2026: 4 unrelated hits (prime units, Bloch's formula, explicit cubic units, Lubin–Tate regulators).
- arXiv API, 7 queries (Zimmert; regulator + "lower bound"; Friedman + regulator; regulator + "explicit formula"; regulator + quasi-Riemann; regulators + number fields + lower bounds; Skoruppa + regulator): nothing using a zero-free strip.
- Zenodo, 4 phrases, and HAL, 2: nothing relevant.

## What would be needed
- Either a zero-based lower bound for the residue of the principal-class zeta function itself, which (3) rules out in general,
- or an upper bound for h strong enough to turn residue bounds for ζ_K into regulator bounds. Through h·R/w = √d κ_K/(2^{r1}(2π)^{r2}), the strip does give effective lower bounds for κ_K, Stark-style. That is a Brauer–Siegel statement about hR, outside this question; the project's effective Brauer–Siegel work (Riemann/Deep/BSExplicit) covers that ground.

Files: this note only. Downloads deleted.
