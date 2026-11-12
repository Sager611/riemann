# explore45 STATE (start 16:53 UTC 12 Nov 2026; write-up by 19:23; end 19:53)
- 16:55 setup. Tool: Baernstein star function / symmetrization. Node: N8 (Levinson concavity step).
- Plan: (1) LP barrier: linewise 2nd moments + any inequality valid for ALL subharmonic u -> optimum = Levinson; extremal u = m*(sigma) (t-independent, convex) is subharmonic.
  (2) numbers: actual Jensen gap of V*Psi at modest T, zeta vs DH. (3) quantized inverse problem: polynomials with same line moments, max zeros in |w|<=1 vs Levinson bound.
- 17:50 finite-T runs done (T=1e4, H=4000, dt=0.01, theta=1/2, P=x, Q=1-1.032x, R=1.246; zeta peak 207 MB, DH 241 MB).
  zeta: V has 0 zeros right of 1/2 (winding), Levinson finite-T bound nu<=0.249 (LP over lines 0.2426 at x=-1);
  loss = concavity 0.144 + Littlewood weights (V zeros in (a,1/2)) 0.106. J_f(-R)=0.179, var 0.233.
  DH: V has 275 right zeros (nu=0.057), bound 0.534, concavity 0.197, weights 0.261, m(0) 0.020.
- next: quantized inverse problem (periodic polynomial competitor with same line second moments, max #roots in |w|<=1).
- 17:46 IP results (D3 in DERIVATION): competitors with same line 2nd moments reach nu = 0.200 (M=12/60), LP 0.2407. Next: Pz=100, verify competitor, lit search.
- 19:10 REPORT.md written (BARRIER). Pz=100 IP: 21 +0.0011 (unconverged), 22 +0.009, 23 +0.023. Competitor verified (verify.py). Big npz deleted. Done.
