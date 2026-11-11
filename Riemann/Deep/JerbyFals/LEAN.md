FORMALIZED: Jerby's Gram-point repulsion conjecture (arXiv 2406.06548 conj:rep, TeX l. 313) is false in Lean 4 + Mathlib. Main theorem: `JerbyLean.not_jerbyConj : ThetaDerivPos → certAL → certAM → certAR → ¬ JerbyConj`. The named hypotheses are the Arb records of n = 1,565,734,083 (logs/certify_cx.log, lines 2–4 and 6) and the classical fact θ′(t) > 0 for t ≥ 7. The given-neighbours reading is refuted from the Arb records alone (`not_jerbyConjEx`). The same holds for n = 2,565,007,238 and n = 3,315,629,581. Also proved: Z′ is Arb's formula; θ is the continuous branch of arg Γ(1/4 + it/2) − (t/2) log π with θ(0) = 0; Z is real. No sorry, standard axioms only, no native_decide. Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (rev 0653561).

[Discovery](../../../discoveries/jerby-repulsion.md) · [REPORT §1](REPORT.md)

## Main theorems (`Lean/Main.lean`)

```lean
theorem JerbyLean.not_jerbyConj (hθ : ThetaDerivPos) (hL : certAL) (hM : certAM) (hR : certAR) :
    ¬ JerbyConj
theorem JerbyLean.not_jerbyConjEx (hL : certAL) (hM : certAM) (hR : certAR) : ¬ JerbyConjEx
theorem JerbyLean.counterexampleA_full (hθ : ThetaDerivPos) (hL : certAL) (hM : certAM) (hR : certAR) :
    ∃ gl g gr : ℝ, IsGram (1565734083 - 1) gl ∧ IsGram 1565734083 g ∧ IsGram (1565734083 + 1) gr ∧
      gl < g ∧ g < gr ∧ (∀ m t, IsGram m t → gl < t → t < gr → m = 1565734083 ∧ t = g) ∧
      (∀ t, IsGram (1565734083 - 1) t → t = gl) ∧ (∀ t, IsGram (1565734083 + 1) t → t = gr) ∧
      Good (1565734083 - 1) gl ∧ Bad 1565734083 g ∧ Good (1565734083 + 1) gr ∧
      |deriv hardyZ g| < 4 * |hardyZ g|
```

Other theorems:
- `not_jerbyConj_B` and `not_jerbyConj_C` prove the same for n = 2,565,007,238 (`certB*`) and n = 3,315,629,581 (`certC*`).
- `counterexampleA`, `counterexampleB` and `counterexampleC` derive the explicit good–bad–good triple from the certificates alone, including (−1)^n Z(g_n) < 0 and |Z′(g_n)| < 4|Z(g_n)|.

## Definitions (`Lean/Defs.lean`)

```lean
def thetaDeriv (t : ℝ) : ℝ := (digamma (1 / 4 + t / 2 * I)).re / 2 - Real.log π / 2
def theta (t : ℝ) : ℝ := ∫ u in (0 : ℝ)..t, thetaDeriv u
def hardyZ (t : ℝ) : ℝ := (exp (I * theta t) * riemannZeta (1 / 2 + t * I)).re
def hardyZ' (t : ℝ) : ℝ := (exp (I * theta t) * (I * thetaDeriv t * riemannZeta (1 / 2 + t * I) +
    I * deriv riemannZeta (1 / 2 + t * I))).re
def IsGram (n : ℤ) (t : ℝ) : Prop := 7 ≤ t ∧ theta t = n * π
def Good (n : ℤ) (t : ℝ) : Prop := 0 < (-1 : ℝ) ^ n * hardyZ t
def Bad (n : ℤ) (t : ℝ) : Prop := ¬ Good n t
def JerbyConj : Prop := ∀ (n : ℤ) (g : ℝ), IsGram n g → Bad n g →
    ((∃ g', IsGram (n - 1) g') ∧ ∀ g', IsGram (n - 1) g' → Good (n - 1) g') →
    ((∃ g', IsGram (n + 1) g') ∧ ∀ g', IsGram (n + 1) g' → Good (n + 1) g') →
    4 * |hardyZ g| < |deriv hardyZ g|
def JerbyConjEx : Prop := ∀ (n : ℤ) (g gl gr : ℝ), IsGram n g → IsGram (n - 1) gl →
    IsGram (n + 1) gr → Bad n g → Good (n - 1) gl → Good (n + 1) gr → 4 * |hardyZ g| < |deriv hardyZ g|
def ThetaDerivPos : Prop := ∀ t : ℝ, 7 ≤ t → 0 < thetaDeriv t
```

`digamma` is Mathlib's `Complex.digamma` (= logDeriv Γ); `Gamma` and `riemannZeta` are Mathlib's Γ and ζ.

## Named hypotheses

`GramCert k a b zl zu := ∃ t ∈ [a, b], θ(t) = kπ ∧ Z(t) ∈ [zl, zu]`. `GramCertD` adds `hardyZ' t ∈ [dl, du]`.

| Hypothesis | Content (each interval = Arb ball [mid − rad, mid + rad]) | Source |
|---|---|---|
| `certAL` | k = 1565734082; t ∈ [568009639.356005286611718227430573, …227431427]; Z ∈ [0.27634005895678409472, 0.27634005895678410528] | certify_cx.log l. 2 |
| `certAM` | k = 1565734083; t ∈ [568009639.698978187581685032968666, …032969334]; Z ∈ [18.39714778558963574, 18.39714778558964426]; Z′ ∈ [73.00385040031973562, 73.00385040031974438] | l. 3; Z′ l. 6 |
| `certAR` | k = 1565734084; t ∈ [568009640.041951088540347534043632, …534044368]; Z ∈ [4.076959407739047812, 4.076959407739048188] | l. 4 |
| `certBL`, `certBM`, `certBR` | k = 2565007237, 2565007238, 2565007239 (Z < 0 at all three, Z′ > 0) | certify_cx2.log l. 2, 3 (Z′ l. 6), 4 |
| `certCL`, `certCM`, `certCR` | k = 3315629580, 3315629581, 3315629582 (Z > 0 at all three, Z′ < 0) | certify_cx3.log l. 2, 3 (Z′ l. 6), 4 |
| `ThetaDerivPos` | θ′(t) > 0 for t ≥ 7, i.e. Re ψ(1/4 + it/2) > log π | classical (Edwards §6.5); mpmath evidence below |

- `check_certs.py` parses the logs and Main.lean: 9 of 9 MATCH.
- The values of ζ, ζ′, Γ and ψ at heights 5.7·10^8 to 1.2·10^9 remain hypotheses: Mathlib has no Riemann–Siegel formula and no ball arithmetic.
- θ′ > 0 on [7, ∞) remains a hypothesis: it needs Stirling or the series for Re ψ, and Mathlib's `Digamma.lean` has only the functional equations.

## Proved in Lean

- Parity of n and n ± 1 by explicit witnesses; good and bad; (−1)^n Z(g_n) < 0; |Z′(g_n)| < 4|Z(g_n)| from the endpoints (e.g. 73.00385040031974438 < 4 · 18.39714778558963574).
- `hasDerivAt_theta`, `hasDerivAt_hardyZ`: θ′ = `thetaDeriv` by the fundamental theorem of calculus. Γ is holomorphic and zero-free on Re s > 0, so ψ is continuous there. Z′ = `hardyZ'` = `deriv hardyZ`.
- `strictMonoOn_theta`, `isGram_unique`, `consecutive`, `unique_consecutive`: under `ThetaDerivPos`, there is one Gram point per index, ordered by index, and g_n is the only Gram point in (g_{n−1}, g_{n+1}). So {g_{n−1}, g_n, g_{n+1}} is a Gram block with N = 2 and g_n is isolated (l. 769). Conjecture 8.1 (l. 777) falls too.
- `jerbyConjEx_imp_jerbyConj` (no hypothesis) and `jerbyConj_imp_ex` (under `ThetaDerivPos`): the two readings agree.
- `Phase.lean`, `gamma_phase`: Γ(1/4 + it/2) = |Γ(1/4 + it/2)| e^{i(θ(t) + (t/2) log π)} for every real t.
- `Phase.lean`, `theta_unique`: every continuous φ with φ(0) = 0 that satisfies this identity is θ.
- `RealZ.lean`, `hardyZ_eq`: (Z(t) : ℂ) = e^{iθ(t)} ζ(1/2 + it). Proof from `gamma_phase`, Mathlib's `riemannZeta_conj`, `completedRiemannZeta_one_sub` and Gammaℝ.

## Build

Run from the repo root (`.lake` is the real Mathlib directory; see below):

```
lake build Riemann.Deep.JerbyFals.Lean.Main
lake build Riemann.Deep.JerbyFals.Lean.Phase
lake build Riemann.Deep.JerbyFals.Lean.RealZ
lake env lean Riemann/Deep/JerbyFals/Lean/Axioms.lean
python3 Riemann/Deep/JerbyFals/Lean/check_certs.py
uv run --no-project --with mpmath python3 Riemann/Deep/JerbyFals/Lean/thetaprime.py
```

Clean rebuild (my oleans deleted first), 11 Nov 2026, one compile at a time under `taskpolicy -b`. Peak lean memory was 423 MB (`top`).

```
✔ [3219/3220] Built Riemann.Deep.JerbyFals.Lean.Defs (58s)
✔ [3220/3220] Built Riemann.Deep.JerbyFals.Lean.Main (109s)
✔ [3220/3220] Built Riemann.Deep.JerbyFals.Lean.Phase (98s)
✔ [3261/3261] Built Riemann.Deep.JerbyFals.Lean.RealZ (100s)
Build completed successfully.
'JerbyLean.not_jerbyConj' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.not_jerbyConjEx' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.jerbyConjEx_imp_jerbyConj' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.counterexampleA_full' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.not_jerbyConj_B' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.not_jerbyConj_C' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.hasDerivAt_hardyZ' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.gamma_phase' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.theta_unique' depends on axioms: [propext, Classical.choice, Quot.sound]
'JerbyLean.hardyZ_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
```

No warnings. `rg 'sorry|native_decide|admit'` over `Lean/*.lean` finds nothing.

## Fidelity (TeX of 2406.06548 vs Lean)

| TeX | Lean | Note |
|---|---|---|
| l. 182–184: "the real function Z(t) = e^{iθ(t)} ζ(1/2 + it)" | `hardyZ` = Re(e^{iθ}ζ); `hardyZ_eq`: the imaginary part is 0 | The abstract (l. 164) writes e^{−iθ}: a sign typo. e^{−iθ}ζ is not real; l. 184 and Edwards use e^{iθ}. |
| l. 186–188: θ(t) = arg Γ(1/4 + it/2) − (t/2) log t, "see [E, I]" | `theta` = ∫₀ᵗ θ′; `gamma_phase` and `theta_unique`: the continuous branch with θ(0) = 0 | "log t" is a typo for log π, as in the cited Edwards and Ivić. With log t, e^{iθ}ζ is not real. Jerby's first bad Gram point g_126 (l. 743) is the standard one. |
| l. 164, 207: θ(g_n) = πn, "the unique solution" | `IsGram n t := 7 ≤ t ∧ θ t = nπ` | θ < 0 on (0, 7] (minimum −3.531 at 6.2898), so for n ≥ 1 the threshold excludes no positive solution. `ThetaDerivPos` gives the uniqueness the TeX presumes. |
| l. 212, 523: good iff (−1)^n Z(g_n) > 0; "Otherwise, it is considered bad" | `Good`; `Bad := ¬ Good` | The witnesses are strictly bad: (−1)^n Z(g_n) < 0 is proved. |
| l. 313–315: conj:rep, "good consecutive neighbours g_{n±1}" | `JerbyConj`: neighbours exist and every Gram point of index n ± 1 is good. `JerbyConjEx`: given neighbours are good. | `JerbyConjEx ⇒ JerbyConj` is proved, so `JerbyConj` is the weakest reading and `¬ JerbyConj` refutes every reading. `¬ JerbyConjEx` needs no θ′ hypothesis. |
| Z′(g_n) | `deriv hardyZ` | Proved equal to `hardyZ'`. |
| l. 769, 777: isolated bad Gram point, Conjecture 8.1 | `counterexampleA_full` | g_n is the middle point of the Gram block {g_{n−1}, g_n, g_{n+1}} with N = 2. |

How the Arb computation matches the Lean definitions:
- θ: certify_cx.py computes θ as Im acb_lgamma(1/4 + it/2) − (t/2) log π. acb_lgamma is the analytic log Γ: real on (0, ∞) and continuous on Re z > 0. It is not the pointwise principal log of Γ(z). So Arb's θ is continuous, is 0 at t = 0 and satisfies the phase identity, and `theta_unique` makes it Lean's θ.
- θ′: `digamma().real/2 − log(pi)/2` is `thetaDeriv`.
- Z and Z′: in the series ζ(1/2 + it + iε), the coefficient z[1] is iζ′(s). So Arb's Z is `hardyZ` and Arb's dZ is `hardyZ'`.
- Gram points: interval Newton (N ⊂ X) gives a root of θ(t) = kπ in N ⊂ [a, b].

Reviewer (a verdict-only reviewer subagent; definitions, statements, TeX lines and the Arb code were pasted inline) first returned "ISSUES — minor":
- The issue: the ∀-only neighbour reading is vacuous where a neighbour does not exist, so the source need not imply it. The reviewer added: "this does **not** undermine the supplied counterexample".
- Fixed: `JerbyConj` now requires the neighbours to exist, and `jerbyConjEx_imp_jerbyConj` is proved.
- Otherwise the reviewer found "the supplied definitions and certificates are faithful".
- Trust boundary, as the reviewer put it: this is a conditional Lean refutation plus external Arb certification.
- Re-check of the fix (second reviewer run): "RESOLVED — existence clauses remove vacuity; ¬JerbyConj refutes the all-neighbours reading and, by contrapositive of JerbyConjEx → JerbyConj, the exists-reading, provided the stated certificate and ThetaDerivPos hypotheses are discharged."

## Consistency of the hypotheses

- Each certificate follows from the Arb computation. Arb prints balls rounded outward, and the Lean intervals are exactly mid ∓ rad (`check_certs.py`).
- The certificates are certified true statements, so they are jointly consistent.
- They also fit `ThetaDerivPos`: the intervals are ordered like their indices, and all lie above 7.
- The logs give |Im e^{iθ}ζ| ≤ 1.8·10⁻⁶⁵ at every point, as `hardyZ_eq` requires. Log line 8 holds an independent mpmath siegelz recheck.
- `ThetaDerivPos` is classical. Evidence from mpmath (`thetaprime.py`, 35 MB, 5 s):
  - θ′ vanishes at t* = 6.2898359888369, where θ(t*) = −3.53097.
  - θ′(7) = 0.05359, the minimum of θ′ over a grid of t ≥ 7.
  - ∂_y Re ψ(1/4 + iy) > 0 on (0, 20]. For all y > 0 it is a sum of positive terms.
  - ∫₀ᵗ `thetaDeriv` equals mpmath siegeltheta(t) to 20 digits at t = 7 and t = 50.
  - θ′(g_n) = 9.15988594071654 at n = 1,565,734,083.

## Not formalized

- ζ, ζ′, Γ and ψ at these heights: the certificate hypotheses.
- θ′ > 0 on [7, ∞): `ThetaDerivPos`.
- The zero-freeness of the Gram blocks (certify_zerofree.py). The refutation does not need it.

## Files (`Riemann/Deep/JerbyFals/Lean/`, 696 Lean lines)

| File | Lines | Content |
|---|---|---|
| `Defs.lean` | 216 | θ′, θ, Z, Z′, calculus, Gram points, good/bad, `JerbyConj`, `JerbyConjEx`, uniqueness, consecutiveness |
| `Main.lean` | 186 | `GramCert`, `GramCertD`, the 9 certificates, counterexamples A–C, `not_jerbyConj`, `not_jerbyConjEx` (imports Defs) |
| `Phase.lean` | 177 | `gamma_phase`, `theta_unique` (imports Defs) |
| `RealZ.lean` | 78 | `hardyZ_eq` (imports Phase) |
| `Axioms.lean` | 39 | `#check`, `#print`, `#print axioms` |
| `check_certs.py` | 26 | certificates = Arb balls |
| `thetaprime.py` | 17 | mpmath evidence for `ThetaDerivPos` |

`.lake` is left as the real Mathlib directory (7.4 GB). The placeholder file was not restored, on the orchestrator's instruction: lilean builds there. Mathlib has no Riemann–Siegel θ, Hardy Z or Gram points (`rg`).

## FORMALIZATIONS.md entry

- Jerby's Gram-point repulsion conjecture (our discovery [jerby-repulsion](discoveries/jerby-repulsion.md)) refuted in Lean 4 + Mathlib: JerbyLean.not_jerbyConj : ThetaDerivPos → certAL → certAM → certAR → ¬ JerbyConj. 2406.06548 conj:rep is stated over Mathlib's ζ, Γ and digamma. θ = ∫₀ᵗ θ′ is proved to be the continuous branch of arg Γ(1/4 + it/2) − (t/2) log π, Z is proved real, and Z′ is proved to be Arb's formula. Named hypotheses: the Arb records of n = 1,565,734,083 and θ′ > 0 on [7, ∞); the given-neighbours reading falls with the Arb records alone. Also n = 2,565,007,238 and 3,315,629,581. Sorry-free, standard axioms, no native_decide, 696 lines. [LEAN.md](Riemann/Deep/JerbyFals/LEAN.md) Reviewer (fidelity): faithful after one minor fix (neighbour existence).

Build and axiom logs (clean sequential rebuild; #print axioms): logs/lean_build.out, logs/lean_axioms.out.
