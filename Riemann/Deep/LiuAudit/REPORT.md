AUDIT INCOMPLETE: no false identity, inequality, constant or step found in the Part II inputs of B_2 (970 checks pass, including brute force over Q(ζ12)); open: the family constructions and per-row estimates of the moment inductions (Thms 10.7, 11.2, 12.16) and the theta mean squares behind Prop 13.4, compressed prose that was not refereed line by line.

# Liu, arXiv 2610.26835v2, Part II: audit of the inputs of the human Hecke B_2 (liuaudit, 11 Nov 2026)

Brief: .scratch/orca/liuaudit-brief.md. Source: arXiv e-print number_fields_v4.tex (13,294 lines). "l." = its line, "hx l." = line of Hexagon 2611.00082v1 number_field_new_v2-1.tex. Part II = l.1562–11748 (§§3–16). Numbering parsed from the TeX (amsart, theorem counter per section).

## Short answer
- Nothing false. All 27 Part II items that 2611.00082 imports exist with the hypotheses it uses (§1).
- Verified (§§2–3): the finite Gauss and correlation identities (brute force over Q(ζ12)); the local correlation table C_{t,k}(J), the local Euler factors, principal box and contour region of §13.3; every exponent identity and inequality of the M2 recursion; the comparison formulas E, E_small, R_κ, θ and the exponent steps of the numerator average (Props 14.3, 14.6); the number-field steps (Poisson over unit orbits, compact dual support, Hecke reflection, pointwise prime bound); the continuation argument; B_∞.
- Not verified:
  - that every transformed family returns to the same coefficient family;
  - the positive enlargements (Jensen over the amplifier pool, row-ball enlargement);
  - the per-row raw-volume bounds, i.e. the analytic core of M1 and M2;
  - the theta coefficients and mean squares (Prop 4.1, Lemma 4.2, Prop 5.2, Lemma 13.3) and the detector Sobolev bound (Lemma 13.9).
  
  For M2 the paper itself lists the first three as "described in the text" (Remark 12.15, l.7540–7573) and re-derives the recurrence compactly in App. C.4 (l.12741–13100); we checked every exponent there, not the constructions. The Q(√−3) models of all four items are in OpenAI's sorry-free Lean proof (Riemann/Deep/QRH: 2,924 files). The number-field versions exist only in prose.
- B_2 and B♭ for Hecke L over number fields therefore remain unrefereed human proofs with no flaw found. They are not established.

## 1. Map (part 1)
Verdicts: ✓ correct; ✓n correct and numerically confirmed; G gap (not refereed or prose only); — not tested.

| # | item | l. | hx l. | hypotheses as used | test | verdict |
|---|---|---|---|---|---|---|
| 1 | Lemma 3.4: allowed generators, reciprocity, unit characters | 1908–1963 | 313, 448, 456 | S ⊇ places over 6 and ∞, Cl(R) = 1 | reciprocity factor measured in Lemma 12.10 test | ✓ (standard) |
| 2 | Prop 4.1, Lemma 4.2: exceptional theta, coefficient characters | 2050–2515 | 313 | KP over F ⊇ μ3 | — | G4 |
| 3 | Prop 5.2: radial functional equation | 2912–2943 | 313 | radial tests | — | G4 |
| 4 | Prop 8.1 + unit normalization (8.1), l.3537 | 3494–3600 | 321, 2007, 2110, 2122 | λ of order dividing 6, Φ(0) = 0 | re-derived | ✓ |
| 5 | Lemma 9.1: compact dual support | 3614–3640 | via 4 | compact Fourier support | re-derived | ✓ |
| 6 | Prop 10.2, Lemma 10.6: returned inverse families | 3940, 4977 | 321 | — | — | G3 |
| 7 | Thms 10.7, 11.2: inverse moments M1 | 5062, 5657 | 332, 3171, 3318 | — | — | G3 |
| 8 | Lemma 12.5: extra coprimality | 6119–6148 | 365, 1194 | q_{R*} ≤ Z^C | identity re-derived | ✓ |
| 9 | Lemma 12.6: reflection | 6149–6213 | 1399 | primitive FE | identity re-derived | ✓ |
| 10 | Lemma 12.7: pointwise prime bound | 6214–6258 | 359, 1231, 1320, 2927 | rows ∉ Θ, 13/18 ≤ κ ≤ 1, β_* ≤ (1+κ)/2 | proof checked | ✓ |
| 11 | §12.5 centering, Lemma 12.3, Prop 12.4 | 5917–6064, 6260–6318 | 1433, 2055 | 6κ−1 ≥ 10/3 | exponents_m2 | ✓n |
| 12 | Lemma 12.8: common coefficient, local Gauss table | 6319–6417 | via 14 | p ∤ h | gauss_q12 | ✓n |
| 13 | Lemma 12.9: full Gauss sums | 6435–6472 | 323, 762, 2313, 2350 | — | gauss_q12 | ✓n |
| 14 | Lemma 12.10: full correlation | 6474–6552 | 323, 1520, 2018, 2091 | allowed generators | corr_q12 | ✓n |
| 15 | §12.7: two transforms (displays at l.6628, 6726, 6767, 6850, 6924, 6953) | 6555–7016 | 1331, 1839, 2291, 2396 | — | exponents_m2; bridges re-derived | ✓n identities; G1 enlargements |
| 16 | §12.8: extraction, exceptional rows, centring | 7017–7327 | via 18 | — | exponents_m2; Kummer step by hand | ✓n |
| 17 | Lemma 12.12: finite-depth loss | 7342–7406 | via 18 | — | proof checked | ✓ |
| 18 | Thm 12.16: product moment M2 | 7575–7753 | 328–352, 1151 | z = 0: principal excluded, any bounded lengths; z > 0: rows ∉ Θ, n1+n2+6κz ≤ M, 13/18 ≤ κ ≤ 1, β_* ≤ (1+κ)/2 | rank and parameter arithmetic | ✓ skeleton; G1 |
| 19 | Lemma 13.1: Gram character sums | 8139–8206 | via 23 | nonprincipal | proof checked | ✓ |
| 20 | C_{t,k}(J) table, (13.9)–(13.10), (13.12): local correlations, local factors, correction, principal box | 8552–8838 | 4295 | p ∉ S | ctk_table, local_factors | ✓n |
| 21 | (13.13): principal term f_η | 8882 | 392, 4445, 4467 | — | definition, Mellin pair | ✓ |
| 22 | Lemma 13.7 detector, Lemma 13.9 Sobolev | 8925–9121, 9149–9313 | 372 | β_* ≥ a_f = 51/100 | log-control part checked | ✓ statement; 13.9: G5 |
| 23 | Prop 13.4: direct estimate | 8503–8540 | 384, 3356 | ℓ ≤ 1/5, b ≥ 0, l_y ≥ 11b/6, l_x > ℓ | L(d) step checked | ✓ given Lemma 13.3 (G4) |
| 24 | §13.4 count (13.16), §13.5 assembly | 9226–9269, 9275–9743 | 3458, 3695, 3783 | 13/18 ≤ κ ≤ 1 | comparison_formulas | ✓n formulas; derivation G5 |
| 25 | Lemma 13.10: contour region | 9314–9340 | 378, 3030, 3051, 3067 | 1/2 ≤ a ≤ 11/12, e ≤ 1/1000 | local_factors | ✓n |
| 26 | Lemma 13.11, Thm 13.15: comparison, continuation | 9641–9707, 9904–9933 | 384, 393 | B > 437/500 | argument checked; B_∞ recomputed | ✓ |
| 27 | §14: Lemma 14.2 AFE, Prop 14.3 numerator, Prop 14.6, Lemmas 14.7–14.8, count (14.3), endpoint E (l.10091) | 9994–10841 | 378, 386, 3339, 3481, 3859 | admissibility, Remark 14.1 | comparison_formulas | ✓n formulas; 14.3, 14.6 exponent steps ✓n given M2; 14.8 ✓ |

Not Part II but on the chain:
- Thm 17.1, the seed β_*(F) ≤ B♭ (l.11794; hx l.115, 3163), is all of Part II plus Part III.
- Thm 18.2 and Prop 19.1 (l.11906, 12045; hx l.4497) were audited sound in ../LiuHecke.

## 2. Falsification (part 2)
| script (code/) | statements | smallest cases | result (logs/) |
|---|---|---|---|
| gauss_q12.py | L12.9: both formulas, phase γ_i χ̄(h/p^{i−1})^i, CRT; L12.8 table at i = 1, 6, 7 | F = Q(ζ12), split primes of norm 13 (i ≤ 7: sums of 13^7 = 6.3·10^7 terms), 37 (i ≤ 4), CRT moduli up to 13^6·37; all valuations of h | 283/283 |
| corr_q12.py | L12.10: Fourier form, j = 0 (= φ(u)·1_{u=v}), closed form χ_{n1}(k)χ̄_{n2}(−k)L_C, complete-support form, c = 6 | u, v over primes of norm 13, 37, 61, 73, N(u)N(v) ≤ 5·10^6, plus u = v = π13^6 | 165/165 |
| exponents_m2.py | §12.5 (l.6280–6300), brace (l.6736–6750), diagonal excess (l.6909–6915), child width and allowance (l.6924–6960), restoration (l.7104–7126), exceptional excess, F_1, F_2 (l.7203–7262), centring sup (l.7316–7325), Thm 12.16 ranks and parameters (l.7609–7726) | symbolic; 60,000 random admissible rational instances; local tables to i < 200 | 30/30 |
| local_factors.py | Gauss-phase identities (l.8622–8631); C_{t,k}(J) series = closed P_j^* for j = 0…5 (l.8552–8636); Lemma 13.10 with C_1 = 56; G_p + χ̄_p(u)H_p (l.9287); σ_j table (l.9289–9298); principal box 907/500 and 437/500 (l.8800–8838) | sextic Gauss sums at q = 13, 37, 61, 73; q = 10^8…10^40 at 300 digits; exact monomial exponents | 14/14 |
| ctk_table.py | C_{t,k}(J) table (l.8578–8590): F(t,k,h) = Σ_{x mod π^{t+k}} g_{π^k}(−x)χ_π(x)^t e(hx/π^{t+k}) equals q^k C_{t,k}(J) ρ^{k−t}, J = v_π(h), ρ = χ_π(h/π^J) | q = 13: t+k ≤ 6, plus (t,k) = (6,1), (7,0) (moduli 13^7); q = 37: t+k ≤ 3; all J ≤ t+k+1, three unit parts | 453/453 |
| comparison_formulas.py | crossing r_*, its bounds, z_M, z_S (l.9226–9249); (13.16) = (14.3) with c = 2c_κ; h, C(B(ℓ)), 2D − 3P, 1/6 ≤ 1−δ ≤ R_κ ≤ 1 − 2δ/3, θ chain, endpoint E and E_small, full exponent (l.10016–10111); Prop 14.3 Hölder, balance, maximal-length gate, ω (l.10313–10441); Prop 14.6 margin 1279/750000 and E_small ≤ C_0 + A_smδ (l.10672–10697); B_∞ (Thm 13.15) | symbolic and rational grids | 25/25 |

Findings:
- **Lemma 12.10 depends on Lemma 3.4.** Its closed form holds exactly once the reciprocity factor χ_{n2}(n1)/χ_{n1}(n2) is inserted. For our generators, which are not allowed generators, this factor was a primitive sixth root of unity in 17 cases. The "allowed convention" of Lemma 3.4 is therefore genuinely used: Lemma 12.10 holds exactly when Lemma 3.4 holds.
- **Lemma 3.4 was not tested numerically.** Testing its isotropic complement would need wild Hilbert symbols above 2 and 3, which did not fit in the time. Its proof is standard: Hilbert reciprocity, local duality, and the dimensions 2(|S_f| + r_F) and |S_f| + r_F, which were re-derived.
- **Four FAILs in the first local-factor run were cancellation artefacts.** At 50 digits, H_p − 1 ≈ 10^{−72} at q = 10^40 lost all precision. At 300 digits, and in exact monomial arithmetic, all four checks pass.
- **The local correlation table is exact.** It needs one reading of the outer Poisson in I_η (l.8031–8035): modulus s·cn³, normalized by q_s, unit phase ρ^{k−t}. With that reading all 453 cases agree, including the edge cases 6 | t (q−1, τ_0 = −1) and 6 | t−1. Together with the series check this derives H_p and G_p from their definition.
- **The bounds are tight but correct.**
  - The σ_j bounds hold with equality at a = 1/2 for every j.
  - Lemma 13.10 holds with equality from the monomial ηTY_0.
  - The principal box worst monomial is (q−1)TY_0V, of exponent −907/500.
- **Cost.** One process at a time, single-threaded, at most 210 MB. gauss_q12 took 6 min, ctk_table 5 min; the others take seconds.

## 3. Referee notes (part 3)
Verdict per step: correct, gap, or error. No step was found to be an error.
- **Arithmetic.**
  - Lemma 3.4: correct.
    - E_S is Lagrangian by global reciprocity.
    - The pairing is alternating since −1 is a square and a cube.
    - Allowed generators are unique modulo U_S^6 by injectivity: the class group of R is trivial.
  - Lemma 12.9: correct. The digit sum gives P^{i−1}, then a field Gauss sum; CRT needs only multiplicativity.
  - Lemma 12.10: correct. Additive orthogonality; at j = 0, u = v. The local cases P−1, −1 and P−2 and the one-nonunit case P^{c−1}(P−1)1_{6|c} were all re-derived.
  - Lemma 12.8: correct. The table is Lemma 12.9 at i = 1, 6, 7. The factor χ_p(u)^{2i} = χ_{p^i}(u)χ_u(p^i) uses reciprocity.
- **Poisson.**
  - Prop 8.1: correct.
    - The unfolding runs over G/Γ.
    - R ⊂ F_S is self-dual of covolume 1.
    - Γ = U_S^6 fixes λ of order dividing 6.
    - The dual origin is the compact-quotient integral.
  - Lemma 9.1: correct (|gℓ_x|_S = x on F_S^1).
  - First bridge (l.6628): correct. Write x = abx_r + brx_a + arx_b; reciprocity cancels χ_a(b)χ̄_b(a).
  - Second bridge (l.6850): its normalization is correct.
  - Amplification identity (App. C, l.12815–12834): correct. Only i ∈ {1, 6, 7} survive, with P^{−i}|Δ|² = P^{−1}, (1−P^{−1})², P^{−1}; numerically confirmed through Lemma 12.9.
- **Reflection and prime sums.**
  - Lemma 12.5: correct (Möbius; the mass Σ_{d|R*} q_d^{−1/2} is subpower).
  - Lemma 12.6: correct. The displayed identity follows from the primitive Hecke FE, with C_k = |Δ_F|Nf_k/(2π)^{2r_F} and Γ(s)^{r_F} (F is totally complex and characters of finite order are trivial at ∞). The imprimitive factors are restored by the double sum over d and h.
  - Lemma 12.7: correct.
    - Borel–Carathéodory on radius r_2 < 2 − β_*; Hadamard three circles; then Cauchy.
    - Mellin inversion of −L′/L shifted to σ_0 = (1+κ)/2 + e.
    - The only hypotheses are nonprincipality, which holds off Θ, and β_* ≤ (1+κ)/2.
- **M2 (§12).**
  - All exponent accounting is correct (exponents_m2). Lemma 12.3, Prop 12.4, Lemma 12.12 (A_d ≤ C(1+C²)^d) and Lemma 12.14 are correct.
  - Ranks: the reflection and centring calls strictly lower R, and a width drop M′ ≤ M − σ/2 lowers b(M) by 2.
  - **G1.** Lemma 12.13 (l.7407–7514) and App. C.4 (l.12741–13100) state three things:
    - every nonexceptional returned norm is again an M2 family at the actual scales;
    - the enlargements (Jensen over the admissible pool, row-ball enlargement, Möbius extension to artificial pairs) are positive;
    - the weighted measures obey δ_edge ≤ C_*(η + ξ + ε_0).

    These are written as compact prose, not checked here. They are fixable in principle: the Q(√−3) model is Lean-checked.
  - Exceptional rows (l.7140–7192): correct. A row induces in the S-supported group iff its exterior valuations meet a prescribed class mod 6 (tame Kummer ramification), forced to v_p(h′) ≥ 4 at the new nonunit primes. The count Z^{(m′−2f)/6} is the sixth-power count.
- **G3, M1** (§§10–11, l.3661–5736): not refereed. It is modelled on OpenAI §17, which is Lean-checked over Q(√−3).
- **G4, theta and local correlations.** Prop 4.1, Lemma 4.2, Prop 5.2 and Lemma 13.3 were not refereed (l.2050–2989, 8238–8502). The local values C_{t,k}(J) (l.8578–8590) and their summation to the six J_j are now verified (ctk_table, local_factors); what stays open is that the theta coefficients are γ_2(c)ϑ(cn³)/(√q_c q_n) (l.8005, from Prop 4.1 and Lemma 4.2). ../LiuHecke found every field-dependent input in the literature for F ⊇ μ12.
- **Comparison (§13).**
  - Correct: local factors, principal box, Lemma 13.10, Lemma 13.1.
  - Lemma 13.7: its statement and log-control proof are standard.
  - Prop 13.4: the steps L(d) − L(0) ≤ −3d/8 and L(0) = l_x/2 + b/12 under l_y ≥ 11b/6 are correct. L(d) itself rests on the mean squares (13.6) and (13.7) (l.8245–8252), which are G4.
  - Lemma 13.11, Thm 13.15 and Lemma 14.8: correct (Mellin and contradiction with the supremum). σ_cont ≥ B > 437/500 is part of admissibility (Remark 14.1, l.10116), so H_η ≠ 0 there.
- **Refined numerator (§14).**
  - Prop 14.3: correct given M2 and Lemma 13.9. The block fourth moment uses M2 with n_1 = n_2 = m, z = z_p, 2m + 6κz_p ≤ 1. The Hölder exponent (3R+1)/4 − gz_p/2, the balance to (3R+1)/4 − θ/2, the maximal-length gate (R − (6κ−2g)L)/4 ≤ 0 and ω ≤ ε/60 were all recomputed.
  - Prop 14.6: its exponent steps are correct: the margin 1279/750000, E_small ≤ C_0 + A_smδ, and the floor and small-row lines of Remark 14.1.
  - Lemma 14.7: its statement matches Thm 12.16.
- **G5.** Lemma 13.9 (detector Sobolev) was not refereed. The derivation of (13.16) from the short and long counts was checked in ../Hexagon2 ("crossing in the old count R_s"), not here.
- **Constants.**
  - B_∞ = 0.8749570194200989461… is recomputed exactly as printed.
  - B♭ was reproduced earlier (../QRHImprove, ../LiuHecke).
  - κ_1 = 0.749774 lies in the admissible range [7497/10000, 3/4].
  - C_1 = 56 holds: the measured ratio |H_p − 1|q^{1+10e} is at most 1.09.
- **Strongest single input.** The z = 0 clause of M2 contains a Lindelöf-on-average fourth moment for the sextic family over F (Remark 12.18, l.7836–7848: Σ_{q_k ≤ U} |S_k(U^{1/2})|⁴ ≪ U^{1+ε}).
  - The paper says this does not follow from the published cubic large sieve; the cubic analogue is open.
  - Its Q(√−3) model, OpenAI's Lemma 18.1 for κ ≥ 3/4, is machine-checked in OpenAI's Lean proof (the width induction in Energy/, ../Kappa/REPORT.md). For F ⊇ μ12 nobody has checked it.

## 4. Literature (part 4)
- **Versions.**
  - arXiv 2610.26835: v1 19 Oct 2026, v2 29 Oct 2026 (latest).
  - Hexagon 2611.00082: v1 only, 4 Nov 2026.
  - arXiv 2610.12234: latest v2, 12 Oct 2026.
- **No erratum and no public doubt.**
  - We searched QuasiRiemannTracker issues and PRs (read-only gh api) for "Hecke" (17 hits), "number field" (2), "Liu" (11), "Dedekind" (2) and "26835" (0). None objects to Liu's number-field results.
  - PR 117 (closed) and PR 114 (open: AxiomMath 31/32, proposed general-field track) cite Liu's number-field bounds as unverified literature.
- **Lean.**
  - Nobody has Lean-checked any of Liu's number-field bounds. Liu's GitHub repository (last push 12 Oct) formalizes 2610.12234 only. For general fields, only Cureton–Thorner's 31/32 is Lean-checked.
  - The same numbers for Dirichlet L-functions are kernel-checked on the tracker: 0.874957069799 (Liu, 8 Oct) and B♭ = 0.874886647040 (Collet, 10 Oct); see ../Hexagon2/QRT-chronology.md. We did not check whether those Lean proofs follow Part II's inductions.
- **Other sources.**
  - arXiv API (Hecke zero-free, quasi-Riemann, au:Liu_Baiying): nothing newer than Chang–Wang 2611.00924 (already in our files).
  - Zenodo ("2610.26835", "2611.00082", Hecke zero-free): no record. HAL: nothing.

## 5. Dependency map (part 5)
Source: rg over DISCOVERIES, RECORDS, COROLLARIES, TREE, NOGO and discoveries/. Columns: what falls if Part II fails, and what survives.

| entry | file | Part II input | survives with |
|---|---|---|---|
| least prime represented by a form of discriminant D ≪ \|D\|^{1.7498+ε} | COROLLARIES l.4 | B♭ = 0.874887 (Thm 17.1) | \|D\|^{1.75+ε} from Hu 7/8 (human, Hexagon 2610.00216; ../Combos/NOTE.md Thm 2); \|D\|^{31/16+ε} from Cureton–Thorner 31/32 (Lean), by the same proof (not rerun) |
| Artin L holomorphic and zero-free on Re s > B_2 | COROLLARIES l.14; TREE R317 | B_2 (2611.00082) | Re s > 31/32 by the same Brauer argument (Lean input); 7/8 with Hu (his §9); abelian over Q: B_2 (Lean, Richey PR 129) |
| Chebotarev least prime (B): (0.73 log 3Δ_L)^{7.8978} | discoveries/chebotarev-numeric.md; TREE N11 | B_2 | part (A), (0.49 log 3Δ_L)^{32}, already certified (Lean input); exponent 8 with Hu (constant not recomputed); L/Q abelian: linnik-numeric 7.8978 (Lean) |
| verdicts "sound given 2610.26835v2": hexaudit2 (B_2), ChangWang (0.874462), V14 B_tw for Hecke | NOGO l.405, 432; ../Hexagon2, ../ChangWang | all of Part II | stay conditional; no number on our pages changes |

- Every other entry cites only Lean-checked half-planes (../HPAudit map) and is unaffected.
- Inside the chain:
  - B_2 needs Thm 12.16 at κ_1 = 0.749774 and the seed β_*(F) ≤ B♭, i.e. all of Part II.
  - With only Thm 13.15 (B_∞ = 0.874957), the hypothesis β_* ≤ (1+κ_1)/2 = 0.874887 fails. A rerun of Hexagon2 at κ = 2B_∞ − 1 = 0.749914 would be needed; this was not done.

## 6. Referee verdict
- Wherever we could check them, the parts of Liu's argument that B_2 uses hold up.
  - Every finite formula comes out exactly right in the smallest field, Q(ζ12).
  - Every exponent identity, every inequality and every constant we recomputed is right.
  - The steps specific to number fields are correct: Poisson summation over unit orbits, compact dual support, reflection, the pointwise prime bound, the local Euler factors and the contour region.
- The heart of the moment estimates is the claim that each double Poisson transform returns to the same family, with per-row volume bounds. It is written in compact prose and was not refereed line by line here. The same is true of the inverse moment, the theta mean squares and the detector Sobolev bound. The numerator average is correct once M2 is granted.
- Over Q(√−3) and for Dirichlet L-functions, the same scheme is machine-checked, and the same numbers are kernel-checked. Over larger fields, nobody has checked it. Its z = 0 clause alone already contains a strong sextic fourth moment that is not in the classical literature.
- Our three entries that rest on it (Artin half-plane, Chebotarev (B), binary forms 1.7498) should stay labelled "human, unrefereed". Each has a Lean-backed fallback, listed in §5.

RH relevance: N7/N11 (load-bearing input of Liu's half-planes B♭ and B_2 for every finite-order Hecke L-function over every number field, and of the Artin, Chebotarev (B) and binary-form entries).

Files: code/{gauss_q12.py, corr_q12.py, ctk_table.py, exponents_m2.py, local_factors.py, comparison_formulas.py}; logs/{gauss_q12, corr_q12, ctk_table, exponents_m2, local_factors, comparison_formulas}.log. Run each script with: uv run --no-project --with numpy (gauss_q12, corr_q12, ctk_table) | sympy (exponents_m2, comparison_formulas) | mpmath (local_factors) python3 code/<script>.py

