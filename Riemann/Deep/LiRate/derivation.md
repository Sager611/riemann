# LiRate: heuristic derivations (rigorous parts: proofs.md)

Notation as in proofs.md: θ_γ = 2 arctan(1/(2γ)), λ_n = Σ_γ 2(1 − cos nθ_γ) (RH), R(n) = λ_n − (n/2)(log n − 1 − log 2π + γ_E).

## H1. Mean square law (RH + pair correlation F(α) = min(|α|,1), heuristic)
(1/N)Σ_{n≤N} 4 cos nθ cos nθ' ≈ 2 sinc(N(θ − θ')), θ_γ − θ_{γ'} ≈ (γ' − γ)/γ². At height t (density ρ = L_t/2π, L_t = log(t/2π)) the
kernel has α-support α_0(t) = N/(t² L_t); per zero ∫ min(|α|,1) r̂ = g(α_0) = α_0/2 (α_0 ≤ 1), 1 − 1/(2α_0) (α_0 ≥ 1).
M_c(N) := (1/N)Σ_{n≤N}(R − 7/4)² ≈ 2∫ ρ(t) g(α_0(t)) dt.  With L frozen: ∫_0^∞ g(1/x²)dx = 5/6 + 1/2 = 4/3, so
M_c(N) ≈ (4/(3π)) T_1 L_1 (T_1² L_1 = N) → (2√2/(3π)) √(N log N) = 0.300 √(N log N). Half the 4/3 is the zero diagonal (F = 1),
half the prime part (|α|); the prime part equals the sum of the resolved spectral lines of H2.
Finite-N (exact density from γ_1 on, pred_ms.py): 9.29, 46.97, 195.6 at N = 10³, 10⁴, 10⁵; /√(N log N) = 0.112, 0.155, 0.182,
0.201 (10⁶), 0.225 (10⁸), 0.249 (10¹²): the asymptotic constant is approached at the speed 1/log N.
Data (centred): 7.08, 30.2, 144.3, i.e. data/pred = 0.76, 0.64, 0.74 (logs/analyze.txt): order and growth right, constant
~30% high at N ≤ 10⁵ (low-height pair statistics; unresolved).

## H2. Explicit oscillatory term (prime side)
Bombieri–Lagarias: λ_n = (archimedean) − lim_X[Σ_{m≤X} (Λ(m)/m) L^{(1)}_{n−1}(log m) − n log X]. Hilb/Fejér for fixed m:
(Λ(m)/m) L^{(1)}_{n−1}(log m) ≈ (Λ(m)/√m)(n/log m)^{1/2} J_1(2√(n log m)) ≈ π^{−1/2} b_m n^{1/4} cos(2√(n log m) − 3π/4),
b_m = Λ(m) m^{−1/2} (log m)^{−3/4}. Zero side gives the same term: S(t) ≈ −π^{−1}Σ Λ(m) sin(t log m)/(√m log m) in −∫ S g_n',
∫_0^∞ sin(t log m) sin(n/t) dt/t² = (π/2)(log m/n)^{1/2} J_1(2√(n log m)) (GR 3.871).
Prediction: F(u) := (R(u²) − 7/4)/√u has lines at ω_m = 2√(log m), complex amplitude (π^{−1/2} b_m/2) e^{iπ/4}.
Data (u ∈ [10, 316], Hann): |C| = 0.1820, 0.1667, 0.0765, 0.1421, 0.1261 at m = 2, 3, 4, 5, 7 vs 0.1820, 0.1667, 0.0765, 0.1421,
0.1259; phases 0.783–0.789 for m ≤ 23 (π/4 = 0.785); all 19 periodogram peaks > 0.03 sit at ω_m (m ≤ 53); median |C| 0.0006.
The sum over m does not converge absolutely (Σ b_m² ~ 2√(log M)): lines are resolved only for p ≲ √N (frequency spacing
√(log p)/p vs 1/√N); the rest blends into the zero diagonal of H1. Not a pointwise formula.

## H3. Size heuristics
RMS: (R − 7/4) ≈ √(M_c(N)) ≈ 0.55 (N log N)^{1/4} asymptotically (0.37 at N = 10⁵); max_{n≤10⁵}|R|/(n log n)^{1/4} = 1.56.
Conjecture (cf. Montgomery for ψ): R(n) = O(n^{1/4+ε}); proved range of the pointwise exponent: [1/4, 1/2] (proofs.md, Barrier).
