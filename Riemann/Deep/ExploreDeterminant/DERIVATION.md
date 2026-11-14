# explore49 DERIVATION (determinant method -> N9 second spacing problem)

## 1. Where lattice points near curves enter 13/84 (Bourgain 1408.5794 = JAMS 2017; brief's "1407.2567" is a wrong id: Witten index paper)
S = Σ_{m~M} e(T F(m/M)), M ≤ √T. Short intervals I of length N; R = ⌈(2M³/(cNT))^{1/2}⌉ (arc J(I) = f''(I)/2 has length > 1/R²);
minor arcs I ∈ 𝓘(Q,ℓ): reduced a/q ∈ J(I), Q ≤ q < 2Q, Q ≥ R; H ≍ NQ/R². Poisson on I: Σ_{h≤H} e(x(I)·(h,h²,h^{3/2},h^{1/2})).
x(I) ∈ [-1/2,1/2]²×[-X3,X3]×[-X4,X4], X3 = X4 = (R/Q)²H^{1/2}; Y = (6H, 6H², 6H^{3/2}, 6H^{1/2}).
Double large sieve (B-I Lemma 2.4): |Σ_I Σ_h|² ≪ A·B1·Π(X_jY_j+1).
- A (FIRST spacing, the h-vectors) = A_6 ≪ δ²H^{10+ε}, δ = Q²/(H²R²): Bourgain's decoupling input, essentially optimal (major arcs).
- B1 (SECOND spacing) = #{(I,I'): |x_j(I) − x_j(I')| < 1/(2Y_j), j ≤ 4}: pairs of rationals a/q, a'/q' with
  x2 = −\bar{4a}/q (mod 1), x1 = −\bar{2a}b/q (mod 1), b = nearest integer to q f'(m0), f''(m0)/2 = a/q,
  x3 ∝ q^{−3/2}(f'''(m0))^{−1/2}, x4 ∝ {q f'(m0)}·x3   [Gauss sum + stationary phase; standard B-I shape].
  So B1 counts integer points (q,a,q',a') close to a 2-dim set in R⁴, defined by modular-inverse and smooth conditions.
Bourgain uses Huxley–Watt 1988 [H-W] for B1: B1 ≪ Δ1Δ2 (M/N)²(Q/R)⁴ = M²R²Q²/N⁶ at N = M T^{−2/7}
 (Δ1 = 1/(X2Y2) = 1/(3H²), Δ2 = 1/(X3Y3) = Q²/(6R²H²)); the newer resonance-curve bound [H1] only extends the M-range.

## 2. Exponent bookkeeping at M = √T (worst case Q = R, so H = N/R, |𝓘| ≍ M/N)
|S|⁶ ≪ M⁶/N³ (short-interval loss, (3.5)) + (M⁵R⁴N/Q⁷)·B1^{1/2}  [A·Π ≍ H^{12}].
B1 ≥ diagonal |𝓘| ≍ M/N  → term D = M^{11/2}N^{1/2}/R³ ∝ N² (increasing in N).
[H-W] off-diagonal M²R⁴/N⁶ → term M⁶/(RN²) ∝ N^{−3/2} (decreasing). D = [H-W] term exactly at N = T^{3/14} (R = T^{1/7}):
 |S|⁶ ≪ T^{17/7} → |S| ≪ M^{1/2}T^{13/84}.  (Checked: Bourgain (3.12)-(3.13).)
Ideal B1 = |𝓘| + |𝓘|²Π_j 1/(X_jY_j) (all four coordinates equidistributed): off-diagonal term = |𝓘|²Q⁴/(324H⁶R⁴) ≍ M²/(N²H⁶),
 its contribution to |S|⁶ ≍ M⁶/N³ = the short-interval term. Then |S|⁶ ≪ T³/N³ + T²N² → N = T^{1/5}, |S| ≪ T^{1/4+3/20}.
 ⇒ the WHOLE gap 13/84 − 3/20 = 1/420 sits in the off-diagonal count of B1 (second spacing), for 6th-power Hölder.


## 3. For ζ the (x2, x3) conditions are integer points on conics (exact algebra)
f(u) = −(t/2π)log u, a/q = f''(m0)/2 = t/(4π m0²) ⇒ m0 = (tq/(4πa))^{1/2}; |γ| = |f'''(m0)|/2 = t/(2π m0³);
q f'(m0) = −(t a q/π)^{1/2}; D := q^{3/2}|γ|^{1/2} = (2π)^{−1/2}(4π)^{3/4} t^{−1/4}(aq)^{3/4}.
⇒ x3 = (2/3)/D and x4 = η/D depend on (a, q) only through P = aq; η = {(taq/π)^{1/2}} (nearest-integer residue), b = round((taq/π)^{1/2}).
x3-closeness |x3 − x3'| < 1/(12H^{3/2}) ⟺ |P − P'| ≲ P/(9 x3 H^{3/2}) ≍ Q^{7/2}/(t^{1/4}H^{3/2}) ≍ R⁴/N² (Q = R).
x2: with A = 4a, γ = [[A,c],[q,d]] ∈ SL2(Z) (d = \bar{A} mod q), likewise γ'; δ := γ'γ^{−1} = [[α,β],[g,e]], (A',q')ᵀ = δ(A,q)ᵀ,
 g = q'd − d'q, and ‖\bar{A}/q − \bar{A'}/q'‖ = ‖g/(qq')‖ (d' mod q' shifts g by qq'); x2-closeness ⟺ |g| ≤ qq'/(12H²) for one choice.
Identity (αe − βg = 1): Aq − A'q' = −B_δ(A,q),  B_δ(x,y) = αg x² + 2βg xy + βe y²,  disc B_δ = −4βg.
⇒ an off-diagonal (x2,x3)-pair = (δ, (A,q)) with |g| ≲ Q²/H² and |B_δ(A,q)| ≲ 4R⁴/N²: integer points on the conics B_δ = n, |n| small.
Per conic (fixed δ, n, B_δ indefinite, disc not a square): solutions = O(τ(n)·log Q) points (unit orbits) in the box — divisor-type
bound, far below Bombieri–Pila's B^{1/2+ε} for degree 2. So for these two conditions a determinant count is dominated by arithmetic.


## 4. Exchange rate: second-spacing saving -> exponent (exact, from Bourgain (3.12)-(3.16) bookkeeping, Q = R binding)
Exponents of T, M = T^α, N = T^n, R = T^{(3α−n−1)/2}, H = N/R:
 trivial E_t = 6α − 3n; diagonal E_d = α + 2n + 3/2; off-diagonal with B1_off ≪ [H-W]·H^{−σ}: E_o = 9α/2 − 3n/2 + 1/2 − (σ/4)(3n − 3α + 1).
 Optimum E_d = E_o: n = (14α − 4 − σ + 3σα)/(14 + 3σ); then |S| ≪ M^{1/2} T^{θ(σ)+ε},
 **θ(σ) = (26 + 5σ)/(12(14 + 3σ))**, independent of α, valid for (6+σ)/(14+3σ) ≤ α ≤ 1/2, 0 ≤ σ ≤ 2 (E_t binds beyond σ = 2).
 σ = 0 ([H-W]): 13/84 (α ≥ 3/7: Bourgain's range ✓). σ = 1 (x1 OR x4 fully used): 31/204 = 0.15196. σ = 2 (all four, ideal): 3/20.
 θ'(σ) = −24·... < 0: every σ > 0 improves 13/84.
V-route (Huxley's resonance curves, x2 tightened by V = N/Q ≍ H, VB_V ≪ V·diag + Δ1Δ2Δ4^{s}(M/N)²(Q/R)⁴), α = 1/2:
 E_dV = 15/8 + 11n/4, E_o(s) as above with σ = s ⇒ n = (7+s)/(34+6s), θ_V(s) = 1/16 + 11(7+s)/(24(34+6s)).
 s = 2/3 (Huxley 2004/H1): 53/342 = 0.15497 (= Bourgain (3.18) ✓, loses to 13/84); s = 1: 37/240 = 0.154167; s = 2: 7/46 = 0.15217.
 Threshold: θ_V(s) < 13/84 ⟺ **s > 3/4**. General tightening V = H^v with perfect counting: ceiling (6+v)/(40+6v) (v=0: 3/20, v=1: 7/46).


## 5. Exact counts (spacing2.c; logs/run0-run6; t = πK; Q = R; q odd; 4 random shuffles of (x1, x2, η) = random model)
Normalization of x(I) is mine (Gauss sum + stationary phase), so constants are mine; exponent conclusions are not.
- Conic identity Aq − A'q' = −B_δ(A,q), αe − βg = 1: verified on every (x2,x3)-pair (≤ 10⁶ per run, 0 failures).
- x2x3/[H-W]-shape Δ1Δ2(M/N)²(Q/R)⁴ = 0.0159 ± 0.0002 for t = 10^{18..22}, ν ∈ {2/7, 0.29825, 0.3}, V ∈ {1, N/Q}, α ∈ {0.45, 0.5}
  (0.0218 at α = 0.42): [H-W]'s bound has exactly the shape of the truth for (x2,x3) — nothing left to gain there.
- all4/([H-W]-shape·H^{−2}) = 0.00104–0.00109 (same runs): the truth is σ = 2 (ideal); Huxley's Δ4^{2/3} shape overestimates by a
  growing factor (all4/(shape·Δ4^{2/3}) = 0.00066 → 0.00034 for t = 10^{18} → 10^{22}).
- true/random: x2x3 1.17→1.24 (saturates, t = 10^{12}→10^{22}), all4 ≈ 1.38: a constant-factor arithmetic excess, no growth.
- off(all4)/diag at ν = 0.3 (N = T^{1/5}): 0.147, 0.1505, 0.1504 (t = 10^{18,20,22}) — balance point confirmed; at ν = 2/7: ~10^{−4}.
- Per curve: generic pairs (disc B_δ = −4βg > 0, non-square) ≤ 2 per conic (δ, n) and ≤ 12 per matrix δ (α = 0.5, 0.45, 0.42).
  Degenerate: square disc (rational line pairs through small-denominator points, e.g. a/q ≈ 1): 0.4–3% of pairs;
  g = 0 translations (q' = q, a' = a + jq; exist only for α < 1/2, arc range > 1): 0 (α = 0.5), 2.1% (0.45), 32% (0.42).
- Smallest case t = 10^{10} (H = 2.04, 77 arcs, 12 off-diagonal (x2,x3)-pairs, 12 distinct conics): 41/51 ~ 44/43 with
  δ = [[169,−540],[41,−131]], Aq − A'q' = 796 = −B_δ(164, 51), disc = 88560.
## 6. Why a determinant count cannot supply σ > 0 (freedoms: curve, coordinates, box, degree, real/p-adic)
(a) First spacing: A_6 ≪ δ²H^{10+ε} is attained by major arcs (Bourgain): no count can beat it.
(b) (x2,x3): [H-W] = truth × 63 (shape constant), so no method gains here; the gap 13/84 → 3/20 is entirely the x1/x4 saving
    (measured 0.066/H², = σ = 2), worth exponent θ(σ) = (26+5σ)/(12(14+3σ)).
(c) x1 = 2b·x2 (mod 1) exactly (\bar{2a} ≡ 2·\bar{4a} mod q): x1 tests the integer part b of √(taq/π) modulo qq'/(2g), x4 its
    fractional part η: both are Diophantine conditions on √(t a'q'/π) along the (x2,x3)-pairs.
(d) Curve-by-curve (conic B_δ = n, strip of δ, hyperbola a'q' = P' for fixed (a,q)): generic curves carry ≤ 2 counted points.
    A determinant bound (Bombieri–Pila B^{1/d+ε}, Heath-Brown p-adic, Salberger, any degree d) is a sum over curves of integers
    ≥ 1 whenever the curve can carry a point, so it is ≥ #generic (x2,x3)-pairs ≍ [H-W]: σ = 0. (D = (d+1)(d+2)/2 > 2 points are
    needed for one vanishing determinant.)
(e) x4 as integer points near y = √(KP) in the variable P (any shear y − kP): g''(P) = √K/(4P^{3/2}) ∈ [0.0156H, H] (Q = R; √K = √2HR³).
    det[1,P_i,m_i] = Π(P_j − P_i)·g''(ξ)/2 + O(4ε) ≥ g'' − 4ε for distinct integers: once 0.0156H > 1 + 4ε (H ≳ 64, asymptotic)
    no degree-1 (divided-difference) determinant can vanish. Degree d > 1 is NOT excluded by curvature (reviewer): y² = (t/π)P is a
    conic; the x4 points sit at distance ≍ 1/H from it (thick regime); not analysed further.
    Removing g'' needs a rational approximation of g''/2 and progressions mod its denominator = a second Bombieri–Iwaniec step,
    not a determinant step. In (P, n = P − P') coordinates ∂²/∂n² ≍ g'' as well.
(f) Only for fixed n is the x4 curve flat (curvature ≈ nH/Q²); there a rounding count over all integers P works, but the
    x2-admissible P are an arithmetic set of density 10^{−5}–10^{−7} per n (runs): the curve count loses the x2 saving (≫ [H-W]).
(g) p-adic: x4 is archimedean; x1, x2 are congruences to the varying modulus q, not to a fixed prime: Heath-Brown's grouping by
    residues mod p only re-sorts the same ≤ 2-point conics.
(h) Degenerate families (translation lines, rational line pairs) are lines: there a 1-D rounding count along the line can save x4
    (curvature ≈ 0.15j at α = 0.42), but they are ≤ 1/3 of B1 (0 at α = 1/2, share ∝ T^{4/5−2α}); the generic part alone keeps σ = 0.
V-route (Huxley): tightening x2 by V = N/Q makes resonance-curve counting possible (s = 2/3) but inflates the diagonal by V;
even perfect counting gives 7/46; beating 13/84 there needs s > 3/4 — a resonance-curve question, not decided here.


## 7. Review (reviewer subagent, 19:33Z)
Algebra correct: θ(σ) = (26+5σ)/(12(14+3σ)) for 0 ≤ σ ≤ 2; V-route θ_V(s) = (32+5s)/(12(17+3s)) (53/342 at s = 2/3, 13/84 at s = 3/4);
conic identity and disc = −4βg; x1 = 2b·x2 (mod 1). Verdict on "BARRIER": too strong — the fibre argument obstructs only separate
curve-by-curve bounds; regrouping across fibres, other coordinates, joint congruence/archimedean arguments stay open; higher-degree
exclusion by curvature withdrawn. Report relabelled NO IDEA with the fibrewise obstruction; numerics labelled as observations.
