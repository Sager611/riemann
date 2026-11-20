T = open('TREE.md').read(); C = open('COROLLARIES.md').read()
pg = 'discoveries/weil-simple-ground-state.md'; PG = open(pg).read()
L = '[WSGSAudit](Riemann/Deep/WSGSAudit/REPORT.md)'
t_old = "also κ(n) brackets giving 8 new terms of OEIS A401172)"
t_new = "also κ(n) brackets giving 8 new terms of OEIS A401172; audited line by line, wsgsaudit, 20 Nov 2026: AUDIT PASSED, fidelity to CCM (3.10), lower bounds on all of L², deflation, chain logic, N₀ at 52 windows and the eight terms, no false step, " + L + ")"
c_old = "eight certified OEIS A401172 terms. [page](discoveries/weil-simple-ground-state.md)"
c_new = "eight certified OEIS A401172 terms; audited line by line (wsgsaudit, 20 Nov 2026: AUDIT PASSED, " + L + "). [page](discoveries/weil-simple-ground-state.md)"
tok = '**RH relevance:**'
assert T.count(t_old) == 1 and C.count(c_old) == 1 and PG.count(tok) == 1 and '**Audit.**' not in PG
p_new = ("**Audit.** Line by line by wsgsaudit (Orca, 20 Nov 2026): AUDIT PASSED, no false step. Fidelity: weil_matrix.py's entries are CCM's QW_λ (arXiv 2511.22755 eq. (3.10)): pole term ±2⟨f, cosh or sinh(x/2)⟩², the archimedean term (identical to −W_R of CCM's own matrix formula), prime terms Λ(n)n^{−1/2} for n < c; 25 entries by quadrature to 2.5·10⁻¹⁴. "
 "Lower bounds on all of L² (lattice identity, minorant, Woodbury with Chekmachev's Lemma 5.3, Cholesky residuals, rscan.c's rounding bound, the 401-node scan) and the rank-one deflation: correct (its reviewer on the chain: CORRECT). "
 "Logic: U < min(L_odd, L₂) forces a simple even ground state; sector min-max values are non-increasing in c; 37 contiguous chain steps 2 → 20; the all-N claim is exactly the hypothesis of CCM's Thm 5.10, N₀ ∈ [30, 160] at 52 windows. "
 "OEIS A401172: all eight terms recomputed; one gap filled (C_c^∞ of the open interval is a form core). Falsification: negative controls (ε above a truncated eigenvalue, a window just outside the first step, K2 = 201 at c = 7, g ⊥ ξ_N) fail as a sound certificate must. "
 "Reruns reproduce the logged margins (c = 2, 29/4 → 31/4, 16 → 133/8, 9, 12, 19, OEIS n = 16; the orchestrator reran 29/4 → 31/4 from a copy: identical); windows 18, 20 and the steps 37/2 → 20 (312–676 MB) were not rerun under the 300 MB cap (same code path; reproduced on 6 Nov). "
 "Nits: U_up strings are 12-digit roundings (every margin ≥ 2); the integer rows c = 2–6, 8, 10 used K2 = 201 (valid, weaker minorant). [WSGSAudit](../Riemann/Deep/WSGSAudit/REPORT.md)\n\n")
open('TREE.md', 'w').write(T.replace(t_old, t_new, 1)); open('COROLLARIES.md', 'w').write(C.replace(c_old, c_new, 1))
open(pg, 'w').write(PG.replace(tok, p_new + tok, 1)); print('ok')
