# grhsub working notes (agent grhsub, 20 Nov 2026)
Inputs: θ = θ₁₄₁ = 0.8733458561067185768 (L(s,χ) ≠ 0, Re s > θ, every Dirichlet χ); 31/32 (ζ_K ≠ 0, Re s > 31/32, every number field K).
Zeros symmetric under ρ ↦ 1 − ρ̄, so every nontrivial zero has 1 − θ ≤ β ≤ θ and |ρ| ≥ 1 − θ.

## Item 1: 2610.25793 (main-revised.tex; line numbers of this file, brief's are +2)
- GRH hits (rg GRH|Riemann|zero-free|quasi|Chebotarev|Serre): l.91 Thm maingeneralGRH; l.182 heading; l.197–200 Lemma generalGRH;
  l.289–295 Prop chebotarev (Serre (20_R), K = Q); l.297–313 proof. No other use: Lemmas generalevend/quadraticrep (l.185–288) use
  Pollack's Thm 1.1 + Pólya–Vinogradov + divisor bounds only (unconditional); 'Arguing as in Lemma generalevend' reuses that count.
- GRH enters only in (eq:usedcheb) π_C(x) ≥ Li(x)/p − A_d x^{1/2} log|pmx|, E = Q(m^{1/p}, ζ_p), G = Gal(E/Q) of order p(p−1)
  (T^p − m irreducible ⇒ [K:Q] = p, and p ∤ ... gives full G), C = {σ_{j,1}: j ≠ 0}, |C|/|G| = 1/p; need (eq:crucialinequality)
  π_C(|m|^{1/6}) ≥ 2 + ω(p|m|).
- Replacement (31/32 for ζ_E). Deuring + Lagarias–Odlyzko explicit formula (LO 1977, Thm 7.1 / §9, unconditional form):
  for x ≥ 2, ψ_C(x) = (|C|/|G|)x − (|C|/|G|) Σ_{ψ∈Ĥ} ψ̄(g) Σ_{ρ_ψ, |γ|<x} x^ρ/ρ + O((|C|/|G|)(log x)(log Δ_E + n_E log x) + log x log Δ_E),
  H = ⟨g⟩, ρ_ψ zeros of the Hecke L(s, ψ, E/E^H) (entire for ψ ≠ 1); Π_ψ L(s,ψ) = ζ_E, no poles off s = 1, so ζ_E ≠ 0 on Re s > 31/32
  ⇒ every ρ_ψ has β ≤ 31/32 (no exceptional zero). Zero sum ≤ x^{31/32} Σ_{ρ: ζ_E(ρ)=0, |γ|<x} 1/|ρ| ≪ x^{31/32} log x (log Δ_E + n_E log x)
  (LO Lemma 5.4: #{ρ: |γ − t| ≤ 1} ≪ log Δ_E + n_E log(|t| + 2); |ρ| ≥ 1/32 for |γ| ≤ 1). Prime powers and ramified primes: ≪ n_E x^{1/2} + log x log Δ_E.
  So |ψ_C(x) − (|C|/|G|)x| ≤ A (x^{31/32} log x (log Δ_E + n_E log x) + n_E x^{1/2}), A absolute; uniform in E.
  This is the shape of Serre's (20_R) proof with x^{1/2} → x^{31/32}; ChebExp (DERIVATION.md) proves the least-prime corollary numerically.
- Discriminant (paper, l.299): ramified primes divide pm; Serre Prop. 6 (Hensel): log Δ_E ≤ (n_E − 1) log(p|m|) + n_E log n_E ≪_p log|m|.
- Conclusion: with x = |m|^{1/6} (or any x ≥ (log|m|)^{33}): ψ_C(x) ≥ x/p − O_p(x^{31/32} log x log|m|) ≥ x/(2p);
  π_C(x) ≥ (ψ_C(x) − O(√x log x))/log x ≫_p x/log x > 2 + ω(p|m|) (ω ≪ log|m|). The rest of the proof is unchanged ⇒ Lemma generalGRH and
  Thm maingeneralGRH hold unconditionally (modulo the Lean-checked 31/32 input): for every fixed d ≥ 2 finitely many generalized
  Euclidean pure orders Z[m^{1/d}].
- Note: the paper's own (log|m|)^3 under GRH: x^{1/2} log|pmx| ≪ (log|m|)^{5/2} ✓. With 31/32, x = (log|m|)^{33}: x^{31/32} log|m| = x^{31/32+1/33}
  = x^{0.99905}·, times log x: still o(x) ✓ (exponent margin 1/1056 − o(1)).

## Item 3: 2610.23660 (On_a_series_related_to_the_Mobius_function.tex)
- GRH hits: l.379–398 (θ(x,χ) ≪ x^{1/2}(log kx)^2, Davenport ch. 20; partial summation; (eq:T) l.396: Σ_{p>X} χ(p)ζ(1,p)/p ≪ (log kX)^2/X^{1/2}, X ≥ k);
  l.400–510 Thm grh (#{l: C(k,l) ≥ 0} ≪ φ(k)/log k; proof X = k(log k)^3); l.294 is the classical zero-free region (Thm thm, unconditional).
  ζ(1,x) = Π_{q≤x}(1 − 1/q)^{−1} (l.146), telescoping 1 + Σ_{p≤Y} ζ(1,p)/p = ζ(1,Y) (l.362). w_p := ζ(1,p)/p.
- With θ: θ(x,χ) ≪ x^θ (log kx)^2 (truncated explicit formula, T = x, |ρ| ≥ 1 − θ, Σ_{|γ|≤x} 1/|ρ| ≪ log x log kx), hence
  (eq:T_θ) T_X(χ) := Σ_{p>X} χ(p) w_p ≪ (log kX)^2 X^{θ−1} (X ≥ k), constants depending on θ only.
- FAILS for the 1 − o(1) claim: the (eq:GRH) step needs (φ(k)/log k)^2 max|T_X(χ)|^2 ≪ φ(k)/log k, i.e. X^{2(1−θ)} ≫ φ(k)(log kX)^4/log k,
  X ≥ k^{1/(2(1−θ))+o(1)} = k^{3.9479...}; but the (eq:BT) step (Selberg sieve on pairs p₂ = p₁ + mk, m < X/k) needs X/k ≤ (log k)^{O(1)}:
  for X = k^A the uncentered bound Σ_l F(l)^2 ≥ (Σ_l F(l))^2/φ(k) ≍ (log k)^2/φ(k) gives no o(φ(k)) (reviewer's rephrasing). Zero
  density does not rescue 1 − o(1) (that needs c ≤ 2 on [1/2, θ]) but shortens the tail to X ≥ k^{1.331}: refinement below.
- SURVIVES (Markov): for every ε > 0, k ≥ k₀(ε): #{1 ≤ l ≤ k: (l,k) = 1, C(k,l) < 0} ≥ (2(1 − θ) − ε) φ(k), 2(1 − θ₁₄₁) = 0.2533082877865628464.
  Proof. X = k^A, A = (1+η)/(2(1−θ)). Lemma lem:lem at X: C(k,l) = 1_{l prime} w_l + F(l) − ζ(1,X)/φ(k) + T_X(k,l), F(l) := Σ_{k<p≤X, p≡l} w_p ≥ 0,
  T_X(k,l) = φ(k)^{−1} Σ_{χ≠χ₀} χ̄(l) T_X(χ). E₁ := {l prime, w_l ≥ ε log k/φ(k)}: |E₁| ≤ π(A₀φ(k)/ε) = o(φ(k)) (paper l.437–441).
  E₂ := {l: |T_X(k,l)| ≥ ε log k/φ(k)}: |E₂| ≤ (φ(k)/(ε log k))^2 max_χ |T_X(χ)|^2 ≪ φ(k)^2 (A log k)^4/(ε^2 (log k)^2 k^{1+η}) = o(φ(k)).
  For l ∉ E₁ ∪ E₂ with C(k,l) ≥ 0: F(l) ≥ (ζ(1,X) − 2ε log k)/φ(k). Σ_l F(l) = ζ(1,X) − ζ(1,k) (each p > k lies in one reduced class).
  Markov: #{...} ≤ φ(k)(ζ(1,X) − ζ(1,k))/(ζ(1,X) − 2ε log k). Mertens ζ(1,x) = e^γ log x (1 + O(1/log x)) ⇒ ratio → (A − 1)/(A − 2εe^{−γ}).
  Hence #{C(k,l) < 0} ≥ φ(k)(1 − 2εe^{−γ})/A − o(φ(k)) ≥ (2(1 − θ) − ε')φ(k). □  (31/32 would give 1/16.)

## Item 2: 2609.30573 (main.tex)
- GRH hits: l.262 Thm Quanticep (GRH), l.268 (GRH = for the Dedekind zeta functions to which Bach–Sorenson is applied via Thm ECDTA),
  l.1376 Thm ECDTA (GRH for ζ_{K̃}, K̃ = K(√m); least p ∤ m unramified in K with Frob ∈ C, p ≪ (log d_K̃ + [K̃:Q])^2).
  No other GRH use (l.245, 336–340 use qualitative Chebotarev only). Proof l.1388–1450: n ≥ 3, M = Π_{ℓ≤n} ℓ, K = Q(E[M]),
  C = {A ∈ G_E(M): det(I − A) ∈ (Z/M)^×} ≠ ∅, m = rad(2N_E M); [K̃:Q] ≤ 2M^4 (E:DegBd), log d_K̃ ≪_E M^4 log M (Serre Prop 6) ⇒ p ≪_E M^8(log M)^2.
- Replacement: chebotarev-numeric (A) (Lean-checked 31/32 for ζ_{K̃} only): for the Galois K̃/Q and the class C̃ of a lift à ∈ Gal(K̃/Q) of
  some A ∈ C there is p unramified in K̃ with Frob_p = C̃ and p ≤ (0.484 log 3d_K̃)^32. m squarefree even ⇒ m ≡ 2 (4), disc Q(√m) = 4m, every
  prime dividing m ramifies in Q(√m) ⊆ K̃ ⇒ p ∤ m = rad(2N_E M); Frob_p|_K = class of A ⊆ C.
- So log p ≤ 32 log log d_K̃ + O(1) ≤ 32(4 log M + log log M) + O_E(1) = (128 + o(1)) n (PNT: log M ~ n). c_E(p) > n > κ log p for κ < 1/128.
  Infinitely many p: same last paragraph. ⇒ Quanticep holds unconditionally for 0 < κ < 1/128 (CM and non-CM).
- Remark: Liu's unrefereed B₂ (chebotarev-numeric (B), exponent 7.8978) would give κ < 1/31.59; GRH κ < 1/8.
- Before: Jones (MR2727656 Appendix) unconditional unboundedness only; unconditional LMO/Kadiri–Ng–Wong least prime ≪ d^{16}
  gives only c_E(p) > (1/4 − ε) log log p i.o. Possible sieve results for CM curves (Iwaniec–Jiménez Urroz type): paper l.268 says
  sieve methods may give more but are not pursued. Non-CM unconditional quantitative: none known to us.

## Item 4: 2611.10489 (Manuscript.tex)
- GRH hits: l.1111 Thm application 2 (modular forms; GRH for Dedekind zetas; if I ≠ [0,π] also Hypothesis A + GRH for Sym^m f ⊗ χ);
  l.2616 LO77 Thm 1.1 under GRH: #{p ∈ Π*(x,[0,π];f): 𝔬 | F(a_f(p))} = h(𝔬)π(x) + O(x^{1/2}(N𝔬)^4 log(xN N𝔬)) (deg ≤ (N𝔬)^4, Ser81 §1.4 disc);
  l.2634 Ser81 §7.2 under GRH: #{p ≤ x: F(a_f(p)) = 0} ≪ x^{7/8}(log x)^{1/2}; l.2643 I ≠ [0,π]: LL26 Thm 1.3 Chebotarev–Sato–Tate (GRH + Hyp A);
  l.2653 Tho21 Thm 1.3 (GRH for symmetric powers) for |Π_good| (l.1901: Π_good = {p ∈ Π: p ≫ 1}). l.249/302 are TeX comments.
  Thm application 1 (l.1019, integers n ≤ x) and Thm main_geometric (l.926) are unconditional; the latter needs only SOME c ∈ (0,1).
- Replacement, I = [0,π] (μ_ST = 1, no Sato–Tate; |Π_good| = π(x) − #{F(a_f(p)) = 0} − O(1)):
  (a') half-plane LO (Lemma H, item 1): error (|C|/|G|)x^{31/32} log(Δ x^n) + log Δ ≪ x^{31/32}(N𝔬)^4 log(xN N𝔬).
  (b') Serre §7.2 with Lemma H: for a root a ∈ O_{K_f} of F, ℓ ≍ (x^{1/32}/log x)^{1/4} split in K_f with large image, C = {tr ≡ a mod 𝔩},
       |C|/|G| ≪ 1/ℓ, n ≤ ℓ^4: #{p ≤ x: a_f(p) = a} ≪ x/(ℓ log x) + x^{31/32}ℓ^3 log(ℓNx) ≪ x^{127/128}(log x)^{O(1)} (GRH: 7/8).
  ⇒ h(x,𝔬) = h_{F,f}(𝔬) + O(x^{−1/128+ε}(N𝔬)^4); h_{F,f}(𝔬) ≥ (N𝔬)^{−4} ⇒ ratio 1 + O(x^{−1/128+ε}(N𝔬)^8) = 1 + o(1) for N𝔬 ≤ x^c, c < 1/1024;
  Σ_{q_𝔩 ≤ y}|h(x,𝔩^n) − h(𝔩^n)| ≪ x^{−1/128+ε} y^{4n+1} = o(1) for y^n ≤ x^c, c < 1/640. Rest of the proof (l.2681–2705) unchanged.
  ⇒ Thm application 2 for I = [0,π] holds unconditionally (given 31/32): E_{Π(x,[0,π];F,f)}[d_{K_f}(F(a_f(p)))^m] ≍ (log x)^{K_m}.
- FAILS for I ≠ [0,π]: joint Chebotarev–Sato–Tate needs zero-free regions for L(s, Sym^m f ⊗ χ) (and Hyp A); our inputs cover
  Dedekind/Hecke/Dirichlet only; unconditional Sato–Tate (Thorner) saves only powers of log x, but (ii)/(sieve dimension) needs x^{−δ}.

## Lemma H (citable): Frei–Sofos–Widmer, arXiv 2610.18491, Thm thrm:generaldelta (TeX l.194–207)
If ζ_L has no zeros with Re s > 1 − δ, δ ∈ (0, 1/2], then for every Galois L/K, class C, x > 2:
|π_C(x, L/K) − (#C/#G) Li(x)| ≤ c₁ (#C/#G) x^{1−δ} log(Δ_L^{1/δ} x^{n_L}), c₁ absolute effective (also MR4610084 Thm 4, no proof).
δ = 1/32 (Cureton–Thorner): error ≤ 32c₁(#C/#G) x^{31/32} log(Δ_L x^{n_L}). Their Cor. cor:quantitative_splitting_qrh: a prime of
each class with N𝔭 ≤ c(δ)(log Δ_L)^{1/δ}(log log Δ_L)^{2/δ}. Used in items 1, 2 (alternative to chebexp (A)), 4 (a'), (b').
Item 3 uses instead the Dirichlet explicit formula (Davenport ch. 19–20) with |x^ρ| ≤ x^θ₁₄₁: ψ(x,χ) ≪ x^θ (log kx)^2.

## Item 3 refinement: zero density raises 2(1 − θ) to 2(2 − θ)/3 (draft; reviewer pending)
(ZD) Montgomery (Topics LNM 227, Thm 12.1; Iwaniec–Kowalski Thm 10.4 type): N_k(σ,T) := Σ_{χ mod k} N(σ,T,χ) ≪ (kT)^{c(σ)(1−σ)}(log kT)^B,
c(σ) = 3/(2 − σ), 1/2 ≤ σ ≤ 1, T ≥ 2. Note c(σ)(1 − σ) ≤ 1 for σ ≥ 1/2.
Claim (MV): for A > A₀ := 3/(2(2 − θ)) = 1.33137..., y ≥ k^A: Σ_{χ≠χ₀} |P_y(χ)|² ≪ y^{−η₁}, P_y(χ) := Σ_{p>y} χ(p)/p, η₁ = η₁(A,θ) > 0.
Proof. Y := y^{10}. Tail p > Y: pointwise |P_Y(χ)| ≪ (log kY) Y^{θ−1} (Q2), so Σ_χ ≪ φ(k) Y^{−2(1−θ)} (log)^2 ≤ y^{1/A − 2.5}(log)^2.
p ∈ (y, Y]: P = ∫_y^Y dθ(t,χ*)/(t log t) (χ* primitive inducing χ; p > k so χ(p) = χ*(p)); θ − ψ ≪ √t contributes ≪ y^{−1/2}
(Σ_χ: φ(k)/y ≤ y^{1/A − 1}). Explicit formula (Davenport ch. 19–20), T₀ = Y²: ψ(t,χ*) = −Σ_{|γ|≤T₀} t^ρ/ρ + O(log(kY) log t) for
t ≤ Y (b(χ) ≪ log k/(1 − θ) as |ρ| ≥ 1 − θ); the O-term contributes ≪ log(kY)/y. Zero term: |∫_y^Y t^{ρ−2}(log t + 1)(log t)^{−2} dt|
≤ 2 y^{β−1}/((1 − θ) log y), so |Z_χ| ≪ (log y)^{−1} Σ_{|γ|≤T₀} y^{β−1}/|ρ|; Cauchy–Schwarz and Σ_{|γ|≤T₀} 1/|ρ| ≪ log²(kY):
Σ_χ |Z_χ|² ≪ log²(kY) · W, W := Σ_χ Σ_{|γ|≤T₀} y^{2(β−1)}/|ρ| ≤ Σ_{U=2^j ≤ T₀} (2/U) Σ_χ Σ_{|γ|≤2U} y^{2(β−1)}  (|ρ| ≥ 1 − θ for U = 1).
Zeros with β < 1/2: ≤ φ(k) U log(kU) y^{−1} per U. β ≥ 1/2 (and β ≤ θ by Q2): y^{−1}N_k(1/2,2U) + 2 log y ∫_{1/2}^θ y^{2(σ−1)} N_k(σ,2U) dσ,
and y^{2(σ−1)} (2kU)^{c(σ)(1−σ)}/U ≤ 2 [(2k)^{c(σ)} y^{−2}]^{1−σ} ≤ 2·[2^3 y^{c(σ)/A − 2}]^{1−σ}; c(σ) ≤ c(θ) < 2A on [1/2, θ], so with
1 − σ ≥ 1 − θ (Q2!) this is ≪ y^{−(2 − c(θ)/A)(1 − θ)}. Summing ≪ log T₀ dyadic U: W ≪ (log y)^{B+3} y^{−η₀}, η₀ = min((2 − c(θ)/A)(1 − θ), 1 − 1/A).
So Σ_χ |P_y|² ≪ y^{−η₀/2} (log kY ≪ log y). □
Then T_X(χ) = ζ(1,X) P_X(χ) + Σ_{q>X} (ζ(1,q)/q) P'_q(χ) (paper l.388), P'_q = Σ_{p≥q}; Cauchy–Schwarz with weights q^{∓η}:
Σ_χ |T_X(χ)|² ≪ (log X)² X^{−η₀/2} + X^{−η} Σ_{q>X} (log q) q^{η−1−η₀/2} = o(1) for η < η₀/4. So |E₂| ≤ (φ(k)/(ε log k))² φ(k)^{−1} Σ_χ |T_X|²
= o(φ(k)) already for X = k^A, A > A₀ (instead of A > 1/(2(1 − θ)) = 3.948). Markov as before: #{C(k,l) < 0} ≥ (1/A − ε)φ(k).
RESULT: ≥ (2(2 − θ)/3 − ε) φ(k) = (0.7511027625955209488 − ε) φ(k). Remarks: a uniform exponent c (Huxley 12/5) gives 2/c = 5/6; the
density hypothesis (c = 2 on [1/2, θ]) would give 1 − ε, the paper's GRH conclusion. Q2 is used only to cap β ≤ θ (power saving at σ → 1).
