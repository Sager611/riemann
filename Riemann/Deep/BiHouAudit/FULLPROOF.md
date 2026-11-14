# Bi–Hou criterion chain, refereed line by line
Source: HAL hal-05838582v1 (PDF dated 23 Oct 2026), pdftotext -layout → `bh.txt` (line numbers below). Repo bzyzz/mertens_bound @ a0b267c.
Notation: α = 1/|ρζ′(ρ)|, ψ = arg ρζ′(ρ), K2(v,T,k) = 2Σ_{0<γ<T} e^{−kγ²} Re(e^{iγv}/(ρζ′(ρ))) = 2Σ α e^{−kγ²} cos(γv − ψ) (l.173–178).

## A. Analytic chain (each step re-derived by us)
| Step | bh.txt | Statement | Our check |
|---|---|---|---|
| M1, M2 | 282–289 | Gaussian averages of M(x)/√x, window \|log x−u\| ≤ 2√(ku) | definitions |
| L4.1 | 305–316 | \|M(x)\| ≤ A√x on open window ⇒ \|M1\| < A | ∫_R e^{w/2−(u−w)²/4k}dw = 2√(πk)e^{u/2+k/4}; strict since mass outside window > 0. OK |
| L4.2 | 319–336 | \|M2−M1\| ≤ 2e^{−u/2+2√(ku)} (u ≥ k) | completed square u/2+3k/4−(d−2k)²/4k OK; erfc tail ≤ ½e^{−a²}, a = √u−√k ≥ 0; they drop the ½ (harmless); f(−r) = e^{−2r}f(r) OK |
| L4.3 | 337–347 | M2 = (1/2πi)∫_(2) e^{k(s²−¼)+u(s−½)}/(sζ(s)) ds | 1/(sζ(s)) = ∫_0^∞ M(e^w)e^{−ws}dw (Re s = 2); (1/2πi)∫e^{ks²+as}ds = e^{−a²/4k}/(2√(πk)); Fubini OK |
| residues | 348–352 | k(s²−¼)+u(s−½) = k(s−½)²+(u+k)(s−½) ⇒ residue e^{−kγ²}e^{iγ(u+k)}/(ρζ′(ρ)) | identity checked; ρ, ρ̄ pair gives 2Re ⇒ K2 at v = u+k. OK |
| path L | 377–389 | L1 σ=¼ (0..T0), L2 t=T0, L3 σ=½+1/u (T0..T1), L4 t=T1, L5 σ=½+(5log t+3)/u (T1..t6), L6 σ=2, t6 = e^{3u/10−3/5} | σ(t6) = 2 exactly OK; σ ≤ 2 on L5 OK |
| P4.4 | 394–426 | Q(u) [β ≤ ½+(5logγ+1)/u ∀γ>0] ⇒ \|M2−K2(u+k,T0,k)\| ≤ (J1+…+J6)/π | orientation: ∫_(2) − ∫_L = 2πiΣRes ⇒ M2 = ΣRes + (1/2πi)∫_L; zeros \|γ\|<T0 inside R (on σ=½ by Platt–Trudgian), T0<γ<T1 left of L3, γ ≥ T1 ≥ 2/u left of L5 by Q(u); s=0, −2n outside (σ<¼); s=1 is a zero of 1/ζ. L6 bound (15/2π²)e^{15k/4+3u/2−kt²}. OK |
| L4.5 B1 | 444,468–473 | 4e^{−u/4}(4e^{150}+T0^{1+c0 log4+11}e^{14}) | exp ≤ e^{−u/4}, 1/\|s\| ≤ 4; t ≤ 4: L3.3; t ≥ 4: L3.2 v=4 (zeros within ¼ have γ < T0+c0 ⇒ on line) OK |
| B2 | 446,474–477 | exp(−kT0²+c0 logT0 log(1/δ)+11logT0+16) | k(σ²−¼) ≤ k(1/u+1/u²) ≤ 1, u(σ−½) ≤ 1, length ≤ 1, L3.2 v=1/δ OK |
| B3 | 448,478–479 | u e^{M+2}Γ(6,kT0²)/(2k⁶) | ∫_{T0}^∞ t¹¹e^{−kt²}dt = Γ(6,kT0²)/(2k⁶) OK; uses (4.2) and 1/\|s\| ≤ 1 |
| B4 | 450,480–483 | exp(−kT1²+16logT1+c0 logT1 log u+18) | exp ≤ e^{4+5logT1−kT1²}; L3.2 v=u (zeros near L4 on line, horizontal dist ≥ 1/u) OK |
| B5 | 452,484–495 | 2e^{18}e^{−kT1²/2}/(kT1) | σ(t)−β ≥ (2−5/(uT1))/u > 1/u under Q(u); \|ds\| ≤ 2dt; needs (k/2)T1² ≥ A log T1 (pdftotext "k2 T12" = (k/2)T1²) and T1 ≥ √(A/k) OK |
| B6 | 454,496–497 | (15/π²)e^{1+3u/2}e^{−kt6²}/(2kt6) | ζ(2)/ζ(4) = 15/π² OK |
| T4.8 | 507–541 | \|K2(u+k,T0,k)\| > 1+E ⇒ (a) x ≤ e^u, \|M(x)\| > √x [¬Q(u), Thm 3.4, Y=e^u] or (b) \|log x−u\| < 2√(ku), \|M(x)\| > (\|K2\|−E)√x | ¬Q(u) ⇒ β0 > ½+5logγ0/u = ½+5logγ0/log Y OK; Q(u): \|M1\| ≥ \|K2\|−ΣB/π−2e^{…} = \|K2\|−E > 1, L4.1 ⇒ (b). Sign case h < −1 covered by \|K2\| (\|M1\| symmetric). Step function ⇒ integer n = ⌊x⌋ ≥ 2 with \|M(n)\| > √n OK |
| T3.4 side cond. | 256–271 | Y > γ0⁵ automatic | β0 < 1 ⇒ 5logγ0 < ½logY ⇒ Y > γ0^{10} OK; (1/Y)∫_0^Y√x dx = (2/3)√Y < √Y OK |

## B. Literature inputs, checked against Pintz's scan (numdam AST_1987__147-148__325_0, pp. 328–330)
- Pintz L4 (p.328): log ζ(s) = Σ_{|γ−t|<c0} log(s−ρ) + θ(11 log t + 14), 0<σ<2, t>4. Bi–Hou L3.1 = its real part. Match.
- Pintz L5 (p.329): |s−ρ| > 1/u ∀ρ ⇒ log|1/ζ(s)| < c0 log t log u + 11 log t + 14. Bi–Hou L3.2 adds v ≥ 1 (safer). Match. Proof sketched only (zero count ≤ c0 log t in |γ−t|<c0 implicit).
- Pintz L7: |ζ(¼+it)| > e^{−150}, |t|<4, "without proof"; Bi–Hou L3.3 replaces it by an Arb check (≥ 0.4556); we recheck below.
- Pintz Thm A (p.329): ζ(β0+iγ0) = 0, β0 > ½, γ0 > 0 ⇒ D(Y) = (1/Y)∫_0^Y|M| > Y^{β0}/γ0⁵ for Y > γ0⁵. Proof (p.330): \|U(Y)\| < (2/3)D(Y), and the residue at ρ0 gives \|U(Y)\| > (2/3)Y^{β0}/γ0⁵. Cor A: β0 > ½ + 5log\|γ0\|/log Y ⇒ D(Y) > √Y. Bi–Hou T3.4 = Thm A + Cor A; side condition Y > γ0⁵ automatic (Y > γ0^{10}). Match.
- Pintz M1 (p.330) identical to Bi–Hou M1 (window e^{u±2√(ku)}).
- Remark 4.9 (StR prints log u for u in (4.1)): irrelevant here, Bi–Hou's (4.1) is exactly what Cor A with Y = e^u needs and what L5's placement needs.
- Platt–Trudgian 2021: RH and simplicity to 3·10¹²; used for \|γ\| < T1 + c0 = 4.3·10⁶ + c0. Fine.

## C. Numerical inputs, recomputed by us
| Input | Bi–Hou (bh.txt) | Ours | Script/log |
|---|---|---|---|
| N(T0), T0 = 1.3·10⁶ | 2 325 572 (623–628) | 2 325 572 rows of Gallego's data below T0, contiguous n, increasing γ | k2_gallego.py → k2_bihou.json |
| γ_N, γ_{N+1} | 1 299 999.5363…, 1 300 000.1829… (627) | 1299999.536327494617836…, 1300000.182865248898197… | same |
| δ | ≥ 0.18286 (628) | 0.182865248898197 | same |
| K2 main n ≤ 6·10⁵ | 1.000305250386499012059 ± 3.2·10⁻²² (639–647) | 1.000305250386499 (float64, exact integer phases) | same |
| K2 tail 6·10⁵ < n ≤ N | ball ±9.43·10⁻¹², Σ\|·\| ≤ 7.39·10⁻¹² (648–656) | value 1.29·10⁻¹⁴, Σ 2e^{−kγ²}\|c_ρ\| = 7.385·10⁻¹² | same |
| K2(v, T0, k) | [1.0003052504 ± 2.30·10⁻¹¹] (626) | 1.000305250386512; rigorous ± 6.9·10⁻⁷ from Gallego's radii (r_γ = 5·10⁻²⁶ ⇒ phase ≤ v r_γ = 1.74·10⁻⁷; Σ2w\|c\| = 3.963) | same |
| B1…B6, tail | log ≤ −8.71e17, −45.39, −23.31, −1495.4, −1360.7, −1.09e18; tail −1.74e18 (669–674) | −8.71243e17, −45.3921, −23.3121, −1495.48, −1360.71, ≤ −e^{2.09e18}; tail −1.74249e18 | e_bound.py → e_bound.log |
| E | ≤ 2.40·10⁻¹¹ | 2.39079·10⁻¹¹ | same |
| hypotheses of T4.8 | (659–668) | all true; A = 82.4935, √(A/k) = 741 590.7 < T1, (k/2)T1² − A log T1 = 126.73 (> 120) | same |
| u ∓ 2√(ku) | …086 669.27, …178 123.80 (676–677) | 3484971332733086669.272949, 3484971332733178123.80372 → ⌈⌉ = …178 124 | same |
| Rem. 7.2 (StR M = 23.05, T1 = 10⁷) | E ≤ 3.2·10⁻¹¹ | 3.17252·10⁻¹¹ | same |
| Table 2 rows k = 3e-9, 1.5e-10, 1e-10 | T0 = 273 596, 1 272 608, 1 566 567 | each is the least integer with E ≤ 10⁻⁶ (E(T0) < 10⁻⁶ < E(T0−1)) | same |

## D. Gallego (Zenodo 23700212 v1.0, 27 Oct 2026) vs Bi–Hou
| Item | Bi–Hou | Gallego |
|---|---|---|
| point | v = 3484971332733132396.5383341489440461 | y = 3484971332733132396.538348536063730810969705971170 (y − v = 1.4387·10⁻⁵) |
| k, T0, N(T0) | 1.5·10⁻¹⁰, 1.3·10⁶, 2 325 572 | 5.7·10⁻¹², 7 999 035.172, 16 622 627 |
| T1, M | 4.3·10⁶, 22.7671 (own, Arb zeros) | 22 494 898, 24.25 (from Platt's LMFDB zeros; min L = −10.2467) |
| sum | K2 = 1.0003052504 ± 2.3·10⁻¹¹ | S = 1.0020700 ± 5·10⁻⁷ |
| E | ≤ 2.40·10⁻¹¹ (tail 2e^{−u/2+2√(ku)}) | < 3.9·10⁻⁴⁶ (StR form, tail 2e^{−u/4}); Bi–Hou form at Gallego's parameters: 1.08·10⁻⁴⁹ (e_bound.py) |
| bound | x₁ < exp(…178 124) | x₁ < exp(B), B = …141 310.42338180385 (smaller window, W = 8 913.9) |
| x₁ definition | min{x > 1 : \|M(x)\| ≥ √x} | least n ≥ 1 with \|M(n)\| > √n (README: same bound for ≥) |
| case ¬Q(u) | Pintz Thm A + Cor A (Y = e^u) | StR Cor D, ε = 1 (\|M\| > e√x), from Pintz Thm A |
| case Q(u) | \|M(x)\| > (\|K2\|−E)√x | one-sided: M(x) > +√x (positive kernel; correct, and also true for Bi–Hou since K2 > 0) |
| zeros | Arb acb_dirichlet_zeta_zeros, 45 digits n ≤ 6·10⁵; cache not published | published (.crho, γ 42 digits, r_γ = 5·10⁻²⁶, c_ρ 30 digits); compared with Platt's LMFDB zeros |
| cites the other? | no (HAL 26 Oct; repo public 15 Oct) | no ("previous best I am aware of: KN"); point not explained ("a given number") |
- Same lattice point (Bi–Hou repo CSV since 15 Oct lists y_base = …396.5383311641); Gallego re-polished it for his k. Not independent searches; independent certifications.
- At Bi–Hou's (k, T0), Gallego's y gives K2 = 0.99976 < 1 (k2_bihou.json, offset gallego_y): each bound needs its own (k, T0); neither paper's point certifies the other's parameters.
