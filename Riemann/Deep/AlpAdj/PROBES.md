# AlpAdj PROBES.md

All runs single-threaded under `bash .scratch/memcap.sh 290 <s>` (PARI via `uv run --no-project --with cypari2`); code/probes.py; peaks: memcap 14 MB, /usr/bin/time max RSS 64 MB (logs/P1_rerun.log).

## P0 Literature (TeX of arXiv:2207.04925, code: .scratch/alpadj/tex.py, texl.py, thmno.py)
- Thorne, Theorem 5.1 (label thm_vanishing_over_CM_field, TeX l.1008-1010; intro theorem l.393-394): CM F, n ≥ 2, (π, χ) regular algebraic cuspidal polarized, r_{π,ι}|G_{F(ζ_{p^∞})} irreducible ⇒ H¹_f(F+, ad r_{π,ι}) = 0. l.400: no condition on the residual representation, 'which might even be trivial'. l.392: Newton–Thorne [New19a] needed r(G_{F(ζ_{p^∞})}) 'enormous' (p-adic image with enough regular semisimple elements), again not residual.
- Allen (arXiv:1411.7661) abstract: 'under an assumption on the residual image'. arXiv:2607.11769 (July 2026): non-polarized version, conditions only on the p-adic representation.

## P1 Weight one (test point a), logs/P1.log
- Dihedral weight-one newforms with a real quadratic field of induction (a_p = 0 at all p ≤ 400 inert in Q(√D), D > 0), N ≤ 150: 14 forms at 13 levels (PARI mfgaloistype dihedral, 4 except N = 145: 8), e.g. N = 39 (D = 13), N = 145 (RM by Q(√5) resp. Q(√29)). RM is detected empirically (ρ ⊗ η_D and ρ have equal traces up to 400), not proved.
- N = 39, l split in Q(√13): l = 17, 23, 29: a_l = 0, χ(l) = −1 (Satake roots ±i); l = 43: a_l = 2, χ(l) = 1 (roots 1, 1: irregular). Every Satake ratio is a root of unity (finite image) ⇒ FreeRatios fails; HT weights (0, 0) ⇒ h.Incr fails; π not regular algebraic ⇒ AutCsd0 fails. The Bellaïche–Dimitrov points are not points of the theorem.

## P2 Exceptional zero (test point b), logs/P2.log
- l = 3: E: y² − 11xy + 3y = x³ (3-torsion point (0,0)), N = 2118 = 2·3·353, a_3 = +1 (split multiplicative), v_3(j) = −3, j ∉ Z (non-CM); also [−11, 6] (N = 8958), [−11, 9] (N = 4722).
- l = 5: Tate normal form E(t) with 5-torsion: t = −25 (N = 4495 = 5·29·31), −20 (N = 6190), −30 (N = 36870), all a_5 = +1.
- These are residually reducible (ρ̄^ss = 1 ⊕ ω), non-CM (TypeI, N3 hold over any L), ordinary with weights (0, 1); ψ_1/ψ_2 = χ_cyc exactly (Tate curve, a_l = 1), so ψ_1 = ψ_2 on G_w ∩ G' and Free (OrdZero.lean:19-20) fails; OrdVan2Alt's ψ_i ≠ ψ_j·χ (OrdVan2.lean:19) fails too. Even without Free there is no counterexample: locally the fixed-weight ordinary line of H¹(G_l, ad⁰) is the image of the unit Kummer line = H¹_f (the class c_q of the Tate period dies in H¹(ad⁰)), and H¹_f vanishes globally (Thorne 2022).

## P3 Trace part (test point e), logs/P3.log
- v_l #Cl(L, l^k): Q(√−6), l = 5: 0, 2, 4, 6, 8, 10 (Z_5-rank 2 = d + 1, d = 1); Q(√−14), l = 3: 0, 2, ..., 10 (rank 2); Q(ζ_8), l = 17: 0, 3, 6 (rank 3 = d + 1, d = 2). l ∤ h(L) in all three.
- So with FREE weight the c = −1 (anticyclotomic) part of Hom(G_{L,S}, Q̄_l) has rank d, the dimension of the trace direction of a Hida family; the archive's OrdAt (diagonal of z vanishes on inertia, OrdZero.lean:27) and IsBK (crystalline = unramified for characters) leave Hom(Cl_L, Q̄_l)^{c=−1} = 0. Fixed weight kills the trace part.

## P4 Euler characteristic (test point f), logs/P4.log
- Greenberg–Wiles for ad~ over L+ with the fixed-weight ordinary condition, odd polarization: Σ_{v|l}(dim L_v − h⁰(G_v)) − Σ_{v|∞} h⁰(G_v, ad~) = d·n(n−1)/2 − d·n(n−1)/2 = 0 for n = 2, 3, 4, d = 1..4 (h⁰(ad~) = h⁰(ad~(1)) = 0 for absolutely irreducible ρ). Free weight: +d·n (Hida's weight space). Balanced: nothing forces H¹_ord ≠ 0.
