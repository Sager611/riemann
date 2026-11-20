R='/Users/adriansager/playground/riemann/'
T=open(R+'TREE.md',encoding='utf-8').read(); C=open(R+'COROLLARIES.md',encoding='utf-8').read()
o="hence lim D(λ)² log(1/λ) = 2 + γ − log 4π ⟺ RH + all zeros simple (COROLLARIES) ([DbarNB](Riemann/Deep/DbarNB/REPORT.md), DISCOVERIES)."
n=("hence lim D(λ)² log(1/λ) = 2 + γ − log 4π ⟺ RH + all zeros simple (COROLLARIES) ([DbarNB](Riemann/Deep/DbarNB/REPORT.md), DISCOVERIES); "
 "the same for every primitive Dirichlet χ under GRH for L(s,χ): lim D_χ(λ)² log(1/λ) = Σ m_ρ²/|ρ|², = 2Re L′/L(1,χ) + log(q/π) + ψ((1+a)/2) ⟺ GRH for χ + simple zeros (nbdir, 20 Nov 2026, routine, three referees CORRECT; lower half de Roton 2006; COROLLARIES, [NBDir](Riemann/Deep/NBDir/REPORT.md)).")
assert T.count(o)==1
line=("- The continuous Nyman–Beurling rate for Dirichlet L-functions (N5; nbdir, 20 Nov 2026; routine): for every primitive Dirichlet character χ mod q ≥ 3, real or complex, of parity a, GRH for L(s,χ) implies "
 "lim_{λ→0} D_χ(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|², where D_χ(λ) = dist_{L²(0,∞)}(1_{(0,1)}, span{A_χ(θ/x) : λ ≤ θ ≤ 1}) and A_χ(y) = Σ_{n≤y} χ(n); hence the limit equals C_χ = 2Re L′/L(1,χ) + log(q/π) + ψ((1+a)/2) "
 "(C_{χ₋₄} = 0.1555679799…) if and only if GRH holds for L(s,χ) and all its zeros are simple (fixed χ; ineffective, not uniform in q). "
 "Proof: the ∂̄ upper bound of our ζ theorem (discoveries/nyman-beurling-continuous-rate.md, audited) with L(s,χ) for ζ (the trivial zero s = 0 of even χ, a possible central zero, one-sided zero counts, two-sided good heights and both half-plane contour families for complex χ); "
 "the lower half is de Roton's (Bull. SMF 134 (2006) 417–445, Thm 2, Selberg class), re-proved with dilated cosine/sine kernels. Before: de Roton's liminf only. "
 "Reviewers: three independent block referees CORRECT, ROUTINE (§8 after one repair); orchestrator reran nbdir_D.py from a copy: identical; numerics at log(1/λ) ≤ 16 consistent (bounds below C_χ, increasing). [NBDir](Riemann/Deep/NBDir/REPORT.md)")
assert 'nbdir, 20 Nov' not in C
if not C.endswith('\n'): C+='\n'
open(R+'TREE.md','w',encoding='utf-8').write(T.replace(o,n,1)); open(R+'COROLLARIES.md','w',encoding='utf-8').write(C+line+'\n'); print('ok')
