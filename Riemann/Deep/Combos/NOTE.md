# Combining the 2026 half-planes with family zero density: proofs

Date 2026-10-20. Brief: `.scratch/orca/combos-brief.md`. The inputs (L) Liu, (O) OpenAI 7/8, the
zero count (Z), the functional equation (FE), the test function and the explicit formula are as in
`Riemann/Deep/ArtinHalfPlane/NOTE.md` §0–§1. One more input:

- **(Hu)** Hu, Hexagon 2610.00216 (13 Oct 2026), Thm 1.1: for every number field F and every unitary
  Hecke character ψ of F, of finite or infinite order, L(s,ψ) ≠ 0 for Re s > 7/8.

(L), (O) and (Hu) are claimed, partly Lean-checked, not refereed. Everything below that is not marked
"(L)" is unconditional. B♭ = 0.874886647039964…, δ♭ = 1 − B♭, 2B♭ = 1.749773294…

## Verdicts

| # | Candidate | Previous record | Best from the 2026 inputs | Verdict |
|---|---|---|---|---|
| 1 | least prime p(Q) represented by a primitive positive definite form of discriminant D | \|D\|^455 (Zaman, PhD thesis 2017, Cor 1.4.1); \|D\|^694 (Thorner–Zaman 2017); \|D\|^A with A unspecified (Hexagon 2610.00178 §54, Oct 2026) | under (L): p(Q) ≪ \|D\|^{2B♭+ε} (§1, Thm 2); under (Hu): \|D\|^{7/4+ε} | **NEW (main; conditional on (L))** |
| 1' | zero density for ring class L-functions | weaker in \|D\| for every σ: Pascadi–Thorner 2508.14888 Thm 1.3; Kowalski–Michel 2002; Weiss/Thorner–Zaman (log-free); see [NOVELTY.md](NOVELTY.md) | Σ_χ N(σ;t) ≪ \|D\|^{2σ(1−σ)+ε}T³ per unit window, unconditional (§1, Thm 1) | **NEW (tool)** |
| 2 | least prime ideal in a ray class | \|d_K\|^694 N𝔪^521 (TZ17 Thm 1.1, K imaginary quadratic) | (\|d_K\|N𝔪)^{3+ε} for K imaginary quadratic (§2, Thm 8); general K: substitution only | new, minor |
| 3 | Goldbach exceptional set | X^{7/12+ε} (Pintz–Althöfer 2610.15559) | N^{1/2}(log N)^100: Li–Liu arXiv 2610.26756 (19 Oct) | known |
| 4 | least QNR, least prime QR | (log q)^8 (Hexagon 178, Lemma 1.6) | (log q)^{1/δ♭} = (log q)^{7.993} | substitution; barrier |
| 5 | short intervals, almost all intervals, BV beyond 1/2 | x^{17/30}, x^{2/15} (Guth–Maynard); 1/2 | unchanged | no |
| 6 | explicit class number bound | 0.655√\|D\|/(πe log\|D\|) (Hexagon 178, Lemma 19.2) | explicit loglog constant: routine | no |
| 7 | Piltz exponents α_k | 1 − 1.421(k−1.18)^{−2/3} (Bellotti 2023), k ≥ 30 | α_k ≤ B♭ for all k (§7, Thm 9) | new, minor (trivial corollary) |

## 0. Notation and imported facts

- D < 0, D = f²d_K, K = Q(√D), O = Z + fO_K. G := Cl(O) ≅ I_K(f)/P_{K,Z}(f) [Cox, *Primes of the
  form x²+ny²*, Prop 7.22]; h = h(D) = |G|.
- For χ ∈ Ĝ and an O_K-ideal 𝔞: χ(𝔞) := χ([𝔞]) if (𝔞,f) = 1, else 0. χ* is the primitive Hecke
  character inducing χ, with conductor 𝔣_χ | fO_K and A(χ*) := |d_K|N𝔣_χ ≤ |D|. Distinct χ have
  distinct χ*.
- L(s,χ) := Σ_{(𝔞,f)=1} χ(𝔞)N𝔞^{−s} = L(s,χ*)E_χ(s) with E_χ(s) = Π_{𝔭|f}(1 − χ*(𝔭)N𝔭^{−s}).
  E_χ has ≤ 2ω(f) factors, no zeros in Re s > 0, and |E_χ| ≤ 4^{ω(f)} ≤ d(f)² on Re s = 1/2.
  So L(s,χ) and L(s,χ*) have the same zeros, with multiplicity, in 0 < Re s < 1.
- m(𝔞) := the largest m ∈ Z_{≥1} with 𝔞 ⊂ mO_K; 𝔞 = m(𝔞)𝔞₀ with 𝔞₀ primitive; m(𝔞)² ≤ N𝔞. If
  (𝔞,f) = 1 then (m(𝔞),f) = 1, (𝔞₀,f) = 1, and χ(𝔞) = χ(𝔞₀), because mO_K ∈ P_{K,Z}(f).
- χ̄(𝔞) = χ(𝔞̄) for (𝔞,f) = 1, because 𝔞𝔞̄ = (N𝔞) ∈ P_{K,Z}(f).
- N_D(σ;t) := Σ_{χ∈Ĝ} #{ρ = β+iγ : L(ρ,χ*) = 0, β ≥ σ, |γ − t| ≤ 1}, with multiplicity. T := |t| + 3.
- d_k(𝔞) = number of ordered factorisations 𝔞 = 𝔞₁⋯𝔞_k; d = d₂. C(ε), c(ε) are effective constants
  depending only on ε and the fixed smooth functions; c₀, C, C₀ are absolute and effective.

Imported facts.
- **(F1)** [Cox Prop 7.20, 7.22] 𝔞 ↦ 𝔞 ∩ O is a norm-preserving bijection from O_K-ideals prime to f
  onto O-ideals prime to f, compatible with I_K(f)/P_{K,Z}(f) ≅ Cl(O). For (m,f) = 1:
  𝔞 ⊂ mO_K ⟺ 𝔞 ∩ O ⊂ mO (because mO_K ∩ O = mO).
- **(F2)** [Cox Thm 7.7 and proof] For a reduced primitive form Q = (a,b,c) of discriminant D
  (|b| ≤ a ≤ c), 𝔞_Q := Za + Z(−b+√D)/2 is a proper O-ideal of norm a, and every class of Cl(O) is
  [𝔞_Q] for some reduced Q. For α = xa + y(−b+√D)/2 one has N(α) = a·Q(x,−y). Hence the integral
  invertible O-ideals in [𝔞_Q]^{−1} are the 𝔟 = α𝔞_Q^{−1}, α ∈ 𝔞_Q∖{0}, each coming from |O^*|
  values of α, with N𝔟 = Q(x,−y); and 𝔟 ⊂ mO ⟺ α ∈ m𝔞_Q ⟺ m | gcd(x,y).
- **(F3)** h(D) ≤ 3|D|^{1/2}log²(3|D|). [h(D) = h(d_K)f[O_K^*:O^*]^{−1}Π_{p|f}(1 − (d_K/p)/p) ≤
  h(d_K)σ(f), σ(f) ≤ f(1 + log f), h(d_K) ≤ (3/π)|d_K|^{1/2}(log|d_K| + 2).]
- **(F4)** [Iwaniec–Kowalski Prop 5.7] For a primitive Hecke character ψ of K:
  #{ρ : |γ − t| ≤ 1} ≤ c₀ log(A(ψ)(|t|+3)²).
- **(F5)** [IK Thm 5.3, Prop 5.4] For primitive ψ ≠ 1 of K, A = A(ψ), s = 1/2 + iu:
  L(s,ψ) = Σ_𝔞 ψ(𝔞)N𝔞^{−s}V_s(N𝔞/√A) + ε(ψ,s)Σ_𝔞 ψ̄(𝔞)N𝔞^{s−1}V_{1−s}(N𝔞/√A), |ε(ψ,s)| = 1,
  V_s(y) = (2πi)^{−1}∫_{(2)} y^{−w}e^{w²}(Γ_C(s+w)/Γ_C(s))dw/w, and for every B:
  |V_s(y)| ≤ C_B min(1, (y/(|u|+2))^{−B}).
- **(F6)** Explicit formula (ArtinHalfPlane NOTE, Lemma 1, first display): w ∈ C_c^∞((1/2,1)), w ≥ 0,
  w ≢ 0, W(s) = ∫w(t)t^{s−1}dt, |W(σ+it)| ≤ C_w(1+|t|)^{−6} on −1 ≤ σ ≤ 2. For primitive ψ of K and
  x ≥ 1: S_ψ(x) := Σ_𝔞 Λ_K(𝔞)ψ(𝔞)w(N𝔞/x) = 1_{ψ=1}xW(1) − Σ_ρ x^ρW(ρ) + R with |R| ≤ C(log A(ψ)+2).
  The sum runs over nontrivial zeros and converges absolutely.

## 1. Primes represented by binary quadratic forms

### Lemma 1 (primitive lattice points)
For Q reduced of discriminant D and N ≥ 1:
P_Q(N) := #{(x,y) ∈ Z² : gcd(x,y) = 1, Q(x,y) ≤ N} ≤ 2 + 8N|D|^{−1/2} + 4N^{1/2}|D|^{−1/4}.

*Proof.* Q(x,y) = a(x + by/(2a))² + |D|y²/(4a). So Q ≤ N forces |y| ≤ V := (4aN/|D|)^{1/2} and puts x
in an interval of length 2(N/a)^{1/2}. Row y = 0 contains only the primitive vectors (±1,0). Each of the
≤ 2V rows with y ≠ 0 holds ≤ 2(N/a)^{1/2} + 1 points. So P_Q(N) ≤ 2 + 2V(2(N/a)^{1/2} + 1) =
2 + 8N|D|^{−1/2} + 2V. Reduction gives |D| = 4ac − b² ≥ 3a², so V ≤ 2·3^{−1/4}N^{1/2}|D|^{−1/4}. □

### Lemma 2 (weighted large sieve for ring class characters)
For 0 < ε ≤ 1/2, N ≥ 1 and complex c_𝔞 supported on (𝔞,f) = 1, N𝔞 ≤ N:

  Σ_{χ∈Ĝ} |Σ_𝔞 c_𝔞χ(𝔞)|² ≤ (30/ε) log²(3|D|) (N + |D|^{1/2}) Σ_𝔞 m(𝔞)^{1+2ε}|c_𝔞|².

*Proof.* Write 𝔞 = m𝔟 with 𝔟 primitive. Since χ(m𝔟) = χ(𝔟), Σ_𝔞 c_𝔞χ(𝔞) = Σ_𝔟 χ(𝔟)B_𝔟 with
B_𝔟 = Σ_m c_{m𝔟}. Orthogonality and Cauchy–Schwarz inside each class give
Σ_χ|Σ_𝔟 χ(𝔟)B_𝔟|² = hΣ_{C∈G}|Σ_{𝔟∈C}B_𝔟|² ≤ h·max_C π_C·Σ_𝔟|B_𝔟|², where π_C counts primitive 𝔟
prime to f with [𝔟] = C and N𝔟 ≤ N. Choose a reduced Q with [𝔞_Q] = C^{−1}; by (F1), (F2),
π_C ≤ P_Q(N). Lemma 1 and (F3) give
h·π_C ≤ 3log²(3|D|)(2|D|^{1/2} + 8N + 4N^{1/2}|D|^{1/4}) ≤ 30log²(3|D|)(N + |D|^{1/2}).
Finally |B_𝔟|² ≤ ζ(1+2ε)Σ_m m^{1+2ε}|c_{m𝔟}|² and ζ(1+2ε) ≤ 1/ε. □

*Remark.* Lemma 2 is a weighted inequality, and the weight m(𝔞)^{1+2ε} cannot be dropped. Take c_𝔞 = 1
on the R ideals (n), n ≤ N^{1/2}, (n,f) = 1: every χ sees the sum R, so the left side is hR² while
Σ|c_𝔞|² = R. An unweighted constant is therefore ≥ hR ≍ hN^{1/2}φ(f)/f, which exceeds N + |D|^{1/2} by
a power of |D| when |D|^ε < N < |D|^{1−ε}. The decomposition 𝔞 = m𝔟 moves this loss into the weight.
For coefficients bounded by divisor functions (every use below) Lemma 3 shows the weight costs only
N^{2ε}; then h ≈ |D|^{1/2} characters behave like the characters to a modulus ≈ |D|^{1/2}. Published
large sieves for this family have constant ≍ N + |D| (Huxley/Schaal/Schumer specialised, as quoted by
Ditchen 1312.1502; Pascadi–Thorner 2508.14888, Thm 1.1).

### Lemma 3 (coefficient sums)
For k ≥ 1, ε > 0, M ≥ 2: Σ_{N𝔞≤M} m(𝔞)^{1+2ε}d_k(𝔞)²/N𝔞 ≤ C(k,ε)M^{2ε}, and the same sum without
the factor 1/N𝔞 is ≤ C(k,ε)M^{1+2ε}.

*Proof.* m(𝔞)^{2ε} ≤ M^ε. Write 𝔞 = m𝔟. Then d_k(m𝔟) ≤ d_k(mO_K)d_k(𝔟) ≤ d_k(m)²d_k(𝔟), and
d_k² ≤ d_{k²}. The number of ideals of norm n is ≤ d(n), so Σ_{N𝔞≤Y}1/N𝔞 ≤ (1 + log Y)², and
Σ_{N𝔟≤Y}d_{k²}(𝔟)/N𝔟 ≤ (1 + log Y)^{2k²}. Hence the first sum is
≤ M^ε Σ_{m≤√M} d_k(m)⁴m^{−1}·(1 + log M)^{2k²} ≤ M^ε(1 + log M)^{k⁴+2k²}. The second sum is ≤ M
times the first. □

### Lemma 4 (mollified second moment)
Let 1 ≤ X ≤ |D| and M_X(s,χ) := Σ_{N𝔡≤X,(𝔡,f)=1} μ_K(𝔡)χ(𝔡)N𝔡^{−s}. For u ∈ R and ε > 0:

  Σ_{χ∈Ĝ, χ≠1} |L(1/2+iu,χ)M_X(1/2+iu,χ)|² ≤ C(ε)|D|^{1/2+ε}X^{1+ε}(|u|+3)^{1+ε}.

*Proof.* Put s = 1/2 + iu.
1. |L(s,χ)| ≤ d(f)²|L(s,χ*)|, and d(f) ≤ C(ε)|D|^ε.
2. Group the χ by their conductor 𝔣 (≤ d(fO_K) ≤ d(f)² groups); inside a group A := A(χ*) is fixed.
   By (F5), L(s,χ*) = S₁(χ) + ε(χ*,s)S₂(χ) with |ε| = 1.
3. Write each 𝔞 in S₁, S₂ as 𝔞 = 𝔢𝔞′, with 𝔢 composed of primes dividing f and (𝔞′,f) = 1. Then
   χ*(𝔞) = χ*(𝔢)χ(𝔞′) and χ̄*(𝔞) = χ̄*(𝔢)χ(𝔞̄′). Reindex 𝔞′ ↦ 𝔞̄′ in S₂. This gives
   S_j(χ) = Σ_𝔢 κ_j(𝔢)N𝔢^{−1/2}S_{j,𝔢}(χ) with |κ_j(𝔢)| ≤ 1, where
   S_{1,𝔢}(χ) = Σ_{(𝔞′,f)=1} χ(𝔞′)N𝔞′^{−s}V_s(N𝔢N𝔞′/√A) and
   S_{2,𝔢}(χ) = Σ_{(𝔞′,f)=1} χ(𝔞′)N𝔞′^{s−1}V_{1−s}(N𝔢N𝔞′/√A).
   Cauchy–Schwarz with Σ_𝔢 N𝔢^{−1/2} = Π_{𝔭|f}(1 − N𝔭^{−1/2})^{−1} ≤ d(f)⁴ gives
   |S_jM_X|² ≤ d(f)⁴Σ_𝔢 N𝔢^{−1/2}|S_{j,𝔢}M_X|².
4. S_{j,𝔢}(χ)M_X(s,χ) = Σ_{(𝔟,f)=1} c_𝔟χ(𝔟), with c_𝔟 independent of χ and
   |c_𝔟| ≤ C₀d(𝔟)N𝔟^{−1/2}. Put P := X|D|^{1/2+ε}(|u|+3)^{1+ε}. If N𝔟 > P then every term of c_𝔟 has
   N𝔞′ > P/X ≥ (|D|(|u|+3))^ε√A(|u|+3). With B = 10/ε in (F5) the tail Σ_{N𝔟>P}|c_𝔟| is
   ≤ C(ε)X(|D|(|u|+3))^{−9} ≤ C(ε)|D|^{−8}. Extend the χ-sum to all of Ĝ and apply Lemma 2 and
   Lemma 3 to the part N𝔟 ≤ P:
   Σ_{χ∈Ĝ}|S_{j,𝔢}M_X|² ≤ C(ε)log²(3|D|)(P + |D|^{1/2})P^{2ε} + 2h·C(ε)|D|^{−16} ≤ C(ε)|D|^ε P^{1+3ε}.
5. Sum over 𝔢, j and the groups, and include step 1: the total is
   ≤ Cd(f)^{14}|D|^εP^{1+3ε} ≤ C(ε)|D|^{15ε}P^{1+3ε} ≤ C(ε)|D|^{1/2+19ε}(X(|u|+3))^{1+5ε} (ε ≤ 1/3).
   Rename ε. □

### Lemma 5 (zero detector)
Let 0.6 ≤ σ ≤ 1 and 2 ≤ X ≤ Y. For (𝔞,f) = 1 put a_𝔞 := Σ_{𝔡|𝔞, N𝔡≤X} μ_K(𝔡); so a_{O_K} = 1,
a_𝔞 = 0 for 1 < N𝔞 ≤ X, and |a_𝔞| ≤ d(𝔞). Let χ ≠ 1 and let ρ = β + iγ satisfy L(ρ,χ*) = 0 and
β ≥ σ. Then
- (I) |Σ_{N𝔞>X} a_𝔞χ(𝔞)N𝔞^{−ρ}e^{−N𝔞/Y}| ≥ 1/4, or
- (II) ∫_R |L(1/2+i(γ+u),χ)M_X(1/2+i(γ+u),χ)| e^{−|u|}du ≥ c₁Y^{σ−1/2}, with c₁ > 0 absolute.

*Proof.* For Re s > 1, L(s,χ)M_X(s,χ) = Σ_𝔞 a_𝔞χ(𝔞)N𝔞^{−s}. With e^{−y} = (2πi)^{−1}∫_{(2)}Γ(w)y^{−w}dw,
Σ_𝔞 a_𝔞χ(𝔞)N𝔞^{−ρ}e^{−N𝔞/Y} = (2πi)^{−1}∫_{(2)} L(ρ+w,χ)M_X(ρ+w,χ)Γ(w)Y^w dw.
Move the line to Re w = 1/2 − β ∈ [−1/2, −1/10]. L(·,χ) is entire with polynomial growth in vertical
strips, and Γ decays exponentially. The only pole crossed is w = 0, with residue
L(ρ,χ)M_X(ρ,χ) = E_χ(ρ)L(ρ,χ*)M_X(ρ,χ) = 0. Hence
e^{−1/Y} + Σ_{N𝔞>X}(…) = (2π)^{−1}∫ L(1/2+i(γ+u),χ)M_X(1/2+i(γ+u),χ)Γ(1/2−β+iu)Y^{1/2−β+iu}du.
By Stirling, |Γ(v+iu)| ≤ C₃e^{−|u|} for −1/2 ≤ v ≤ −1/10. If (I) fails, the left side has modulus
≥ e^{−1/2} − 1/4 > 0.35, which gives (II) with c₁ = 0.7π/C₃. □

### Theorem 1 (zero density for ring class L-functions; unconditional)
For every ε > 0, every D < 0 and every t ∈ R:
- (a) N_D(σ;t) ≤ C(ε)|D|^{2σ(1−σ)+ε}T³ for 0.6 ≤ σ ≤ 1;
- (b) N_D(σ;t) ≤ C|D|^{1/2}log³(3|D|)log T for every σ.

Height dependence: (a) is per unit window. Summing over the windows |t − k| ≤ 1, |k| ≤ T, gives
Σ_χ N(σ,T,χ) ≤ C(ε)|D|^{2σ(1−σ)+ε}T⁴ for T ≥ 3, where N(σ,T,χ) counts the zeros of L(s,χ*) with
β ≥ σ and |γ| ≤ T. Only the discriminant aspect is claimed. The T-aspect is crude (Pascadi–Thorner's
family bound has T^{8(1−σ)/(3−2σ)}); Theorem 2 needs only bounded heights, since |W(ρ)| ≪ (1+|γ|)^{−6}.

*Proof.* (b) follows from (F3) and (F4). For (a) we may assume |D| ≥ D₀(ε); otherwise use (b).
Put ε₀ := ε/100, X := |D|^{ε₀}, Y := |D|^{σ−1/2} ≥ |D|^{1/10}, Z := 2Y log²Y. The trivial character
has ≤ c₀log(|D|T²) zeros in the window. For every other χ with a zero ρ in the window with β ≥ σ,
choose one such zero ρ_χ = β_χ + iγ_χ, and let 𝒳 be the set of these χ. Then
N_D(σ;t) ≤ c₀log(|D|T²)(|𝒳| + 1). Lemma 5 splits 𝒳 = 𝒳_I ∪ 𝒳_II.

*Class II.* Square (II), apply Cauchy–Schwarz (∫e^{−|u|} = 2), substitute v = γ_χ + u, and use
e^{−|v−γ_χ|} ≤ e^{1−|v−t|}. With Lemma 4:
c₁²Y^{2σ−1}|𝒳_II| ≤ 2e∫e^{−|v−t|}Σ_{χ≠1}|L(1/2+iv,χ)M_X(1/2+iv,χ)|²dv ≤ C(ε)|D|^{1/2+2ε₀}T^{1+ε₀}.
Since (σ − 1/2)(1 − 2σ) + 1/2 = 2σ(1−σ), this gives |𝒳_II| ≤ C(ε)|D|^{2σ(1−σ)+2ε₀}T².

*Class I.* The number of ideals of norm n is ≤ d(n) and d(𝔞) ≤ d(N𝔞)². So the tail
Σ_{N𝔞>Z} d(𝔞)N𝔞^{−σ}e^{−N𝔞/Y} ≤ Σ_{n>Z} d(n)³n^{−1/2}e^{−n/Y} ≤ 1/16 once |D| ≥ D₀(ε).
Fix φ ∈ C_c^∞((1/2,2)), 0 ≤ φ ≤ 1, with Σ_{j∈Z}φ(v/2^j) = 1 for v > 0. Let
𝒥 := {j : 2^{j+1} > X, 2^{j−1} < Z} and J := |𝒥| ≤ 2log Z. Every 𝔞 with X < N𝔞 ≤ Z is seen only by
j ∈ 𝒥. So for χ ∈ 𝒳_I, (I) and the tail bound give |Σ_{j∈𝒥}F_j(χ,ρ_χ)| ≥ 1/4 − 2/16 = 1/8, where
F_j(χ,ρ) := Σ_𝔞 a_𝔞e^{−N𝔞/Y}φ(N𝔞/2^j)χ(𝔞)N𝔞^{−ρ}. Choose j ∈ 𝒥 with |F_j(χ,ρ_χ)| ≥ 1/(8J), and
let 𝒳_{I,j} be the set of χ with this j.

*Removing β.* For v > 0, φ(v/2^j)v^{−β} = 2^{j(σ−β)}v^{−σ}Φ_β(v/2^j) with Φ_β(v) := φ(v)v^{σ−β}.
Mellin inversion gives Φ_β(v) = (2π)^{−1}∫Φ̂_β(τ)v^{−iτ}dτ with Φ̂_β(τ) = ∫_0^∞Φ_β(v)v^{iτ−1}dv. Two
integrations by parts in log v give |Φ̂_β(τ)| ≤ C_φ(1+|τ|)^{−2}, uniformly for 0 ≤ β − σ ≤ 1/2.
Since |2^{j(σ−β)}| ≤ 1 and Φ_β(N𝔞/2^j) = 0 unless 2^{j−1} < N𝔞 < 2^{j+1},
|F_j(χ,ρ)| ≤ (2π)^{−1}∫|Φ̂_β(τ)||P_j(χ,γ+τ)|dτ, with
P_j(χ,t′) := Σ_{2^{j−1}<N𝔞<2^{j+1}} a_𝔞e^{−N𝔞/Y}χ(𝔞)N𝔞^{−σ−it′}. This no longer depends on β.

*k-th power.* If 2^{j−1} < |D|^{1/2}, let k be the least integer with 2^{(j−1)k} ≥ |D|^{1/2};
otherwise k := 1. Since 2^{j−1} > X/4 ≥ |D|^{ε₀/2}, k ≤ ⌈1/ε₀⌉. Put M := 2^{(j+1)k}. Then
|D|^{1/2} ≤ M ≤ 4^k|D|^{1/2}Z: if k ≥ 2, minimality gives 2^{(j−1)(k−1)} < |D|^{1/2}; if k = 1,
M = 2^{j+1} ≤ 4Z. Also P_j^k(χ,t′) = Σ_{M/4^k<N𝔞<M} b_𝔞χ(𝔞)N𝔞^{−σ−it′} with |b_𝔞| ≤ d_{2k}(𝔞).
Hölder and the uniform bound on Φ̂_β give
(8J)^{−2k} ≤ |F_j|^{2k} ≤ (2π)^{−2k}(πC_φ)^{2k−1}C_φ∫(1+|τ|)^{−2}|P_j(χ,γ_χ+τ)|^{2k}dτ.
Sum over χ ∈ 𝒳_{I,j}. Since γ_χ + τ ∈ I_τ := [t+τ−1, t+τ+1],
|𝒳_{I,j}| ≤ (4JC_φ)^{2k}∫(1+|τ|)^{−2}Σ_χ sup_{t′∈I_τ}|P_j^k(χ,t′)|²dτ.

*Sobolev and the large sieve.* On an interval I of length 2, sup_I|G|² ≤ ∫_I((3/2)|G|² + |G′|²). For
each fixed t′, Lemma 2 applies to the coefficients b_𝔞N𝔞^{−σ−it′} and (−i log N𝔞)b_𝔞N𝔞^{−σ−it′}.
With M ≥ |D|^{1/2} and Lemma 3:
Σ_χ sup_{I_τ}|P_j^k|² ≤ 2(3/2 + log²M)(30/ε₀)log²(3|D|)·2M·(M/4^k)^{−2σ}Σ_{N𝔞<M}m(𝔞)^{1+2ε₀}d_{2k}(𝔞)²
≤ C(ε)M^{2−2σ+4ε₀}.
Hence |𝒳_{I,j}| ≤ C(ε)J^{2k}M^{2−2σ+4ε₀} ≤ C(ε)(|D|^{1/2}Z)^{2−2σ+5ε₀}. Since Z ≤ |D|^{σ−1/2+ε₀},
|𝒳_I| ≤ J·max_j|𝒳_{I,j}| ≤ C(ε)|D|^{2σ(1−σ)+10ε₀}.

Adding the classes, N_D(σ;t) ≤ c₀log(|D|T²)(C(ε)|D|^{2σ(1−σ)+10ε₀}T² + 1), which is (a). □

### Theorem 2 (least prime represented by a binary quadratic form)
Assume (L). For every ε > 0 there is an effective C(ε) such that every primitive positive definite
binary quadratic form Q of discriminant D < 0 represents a prime p ∤ D with

  p ≤ C(ε)|D|^{2B♭+ε},  2B♭ = 1.749773294…

Moreover, for x ≥ C(ε)|D|^{2B♭+ε} the number of primes p ∤ D in (x/2, x] represented by Q is
≥ c_w x/(h(D)log x). Assuming only (Hu), the same holds with 2B♭ replaced by 7/4.

*Proof.* Equivalent forms represent the same primes, so let Q be reduced, and put C := [𝔞_Q]^{−1} ∈ G.
Call a prime ideal 𝔭 ⊂ O_K *good* if N𝔭 = p is prime, p ∤ D and [𝔭] = C. For a good 𝔭, (F1) and (F2)
give 𝔭 ∩ O = α𝔞_Q^{−1} with Q(x,−y) = N(α)/a = p, so Q represents p. Put
S(x) := Σ_{(𝔞,f)=1, [𝔞]=C} Λ_K(𝔞)w(N𝔞/x) = h^{−1}Σ_χ χ̄(C)S_χ(x), with
S_χ(x) := Σ_{(𝔞,f)=1} Λ_K(𝔞)χ(𝔞)w(N𝔞/x). Only powers of the 𝔭 | f separate S_χ from S_{χ*}, so
|S_χ − S_{χ*}| ≤ 2‖w‖ω(f)log x. By (F6) with A(χ*) ≤ |D|,

  hS(x) = xW(1) − Σ_χ χ̄(C)Σ_{ρ∈Z(χ*)} x^ρW(ρ) + O(h(log|D| + ω(f)log x)).   (1)

*The zero sum.* |W(ρ)| ≤ 2⁶C_w(1+|k|)^{−6} when |γ − k| ≤ 1, k ∈ Z. By (L) every zero has β ≤ B♭.
With x^β ≤ x^{1/2} + 1_{β>1/2}log x∫_{1/2}^{β}x^σdσ,

  Σ_χΣ_ρ x^β|W(ρ)| ≤ CΣ_k (1+|k|)^{−6}[x^{1/2}N_D(0;k) + log x∫_{1/2}^{B♭}x^σN_D(σ;k)dσ].   (2)

Write x = |D|^{2B♭+ε}y with y ≥ 1, and let |D| ≥ D₀(ε). Put T_k = |k| + 3.
- By Theorem 1(b), x^{1/2}N_D(0;k) ≤ x|D|^{−0.37}log T_k, and x^σN_D(σ;k) ≤ x|D|^{−0.19}y^{−0.4}log T_k
  for 1/2 ≤ σ ≤ 0.6, because x^{−0.4} ≤ |D|^{−0.69}y^{−0.4}. Keep the factor y^{−0.4}: the integral in
  (2) carries log x = log(|D|^{2B♭+ε}y), and y^{−0.4}log x ≤ C log(3|D|) for every y ≥ 1.
- For 0.6 ≤ σ ≤ B♭, Theorem 1(a) with ε′ = εδ♭/4 gives
  x^σN_D(σ;k) ≤ C(ε)x(|D|^{2σ}/x)^{1−σ}|D|^{ε′}T_k³ ≤ C(ε)x(|D|^{−ε}y^{−1})^{δ♭}|D|^{ε′}T_k³,
  since |D|^{2σ} ≤ |D|^{2B♭} and 1 − σ ≥ δ♭. As y^{−δ♭}log(|D|³y) ≤ C(ε)|D|^{ε′}, the integral term
  in (2) is ≤ C(ε)x|D|^{−εδ♭/2}T_k³.

Summing (2) over k gives Σ_χΣ_ρ x^β|W(ρ)| ≤ C(ε)x|D|^{−c(ε)}. By (F3) the error term in (1) is
≤ x|D|^{−1}. So S(x) = (xW(1)/h)(1 + O_ε(|D|^{−c(ε)})).

*Bad terms.* Ideals above primes dividing d_K contribute ≤ 4‖w‖log|D|·log x. Prime powers 𝔭^k with
k ≥ 2 and primes with N𝔭 = p² contribute ≤ 6‖w‖√x. Both are ≤ (xW(1)/h)|D|^{−0.2}, because
h ≤ 3|D|^{1/2}log²(3|D|) and x ≥ |D|^{1.7}. So for |D| ≥ D₀(ε),
Σ_{good 𝔭, x/2<N𝔭<x} log N𝔭·w(N𝔭/x) ≥ xW(1)/(2h). Each p lies under at most two good 𝔭, so there
are ≥ xW(1)/(4h‖w‖log x) good primes p.

*Small discriminants.* Let |D| < D₀(ε), and use Theorem 1(b) for every σ in (2). The zero sum in (1)
is ≤ Cx^{B♭}(1 + log x)D₀^{1/2}log³(3D₀), the error term in (1) is ≤ Ch·log(3D₀)log x, and the bad
terms are ≤ C(log D₀·log x + √x). Since 1 ≤ h ≤ 3D₀^{1/2}log²(3D₀), the first two are ≤ xW(1)/8 and
the third is ≤ xW(1)/(8h) once x ≥ C′(ε), an effective constant depending only on D₀(ε). Then the
good terms are ≥ xW(1)/(2h) as before. So both claims hold for |D| < D₀(ε) once C(ε) ≥ C′(ε).
Under (Hu) replace B♭ by 7/8 throughout. □

### Corollary 3 (ring class PNT in a power range)
Assume (L). For every ε > 0 there is c(ε) > 0 such that, uniformly for |D| ≤ x^{1/(2B♭)−ε}
(1/(2B♭) = 0.571503…) and C ∈ G,

  π_C(x) := #{𝔭 : N𝔭 ≤ x prime, (𝔭,f) = 1, [𝔭] = C} = (Li(x)/h(D))(1 + O_ε(x^{−c(ε)})).

So #{p ≤ x : Q represents p} = (Li(x)/h(D))(1 + O(x^{−c})) if C ≠ C^{−1}, and half of that if
C = C^{−1}, up to the ≤ log|D| primes dividing D. Before, an asymptotic was known only for
|D| = x^{o(1)}: Thorner–Zaman 2019, Thm 1.4, has relative error exp(−c log x/log|D|) plus an
exceptional-zero term, which Hexagon 2610.00178 §54 removes using (O). Under (L), the published family
density of Pascadi–Thorner (2508.14888, Thm 1.3) would give only |D| ≤ x^{(3−2B♭)/3−ε} = x^{0.4167−ε},
for fundamental D (not stated there).

*Proof.* Put η := εδ♭/100, y := x^{1−η}, T₀ := x^{2η}. Take the smooth weights g_± of the proof of
Theorem D′ in the ArtinHalfPlane NOTE. They satisfy 0 ≤ g₋ ≤ 1_{(y,x]} ≤ g₊ ≤ 1, have supports in
[y/2, 2x], and their Mellin transforms W_±(s) = ∫g_±(xt)t^{s−1}dt obey xW_±(1) = x(1 + O(x^{−η})),
|W_±(s)| ≤ 16 for Re s ≥ 1/8, and |W_±(σ+it)| ≤ C_m(x^η/|t|)^m for |t| ≥ 1. Run (1)–(2) with W_±;
the remainder in (F6) becomes O(x^η(log|D| + 2)).
- Zeros with |γ| > T₀ contribute ≤ C_mx^{B♭+ηm}Σ_{|k|>T₀}|k|^{−m}N_D(0;k). With m = ⌈2/η⌉ this is
  ≤ x^{1/2}.
- Zeros with |γ| ≤ T₀ contribute ≤ 16(2T₀ + 3)max_{|k|≤T₀+1}[x^{1/2}N_D(0;k) + log x∫x^σN_D(σ;k)dσ].
  Here |D|^{2σ} ≤ |D|^{2B♭} ≤ x^{1−2B♭ε}, so Theorem 1 gives x^σN_D(σ;k) ≤ C x^{1−2B♭εδ♭+ε′}T₀³ for
  σ ∈ [0.6, B♭]. Also x^{0.6}|D|^{1/2+ε′} ≤ x^{0.89} and x^{1/2}|D|^{1/2+ε′} ≤ x^{0.79}.
So Σ_{[𝔞]=C} Λ_K(𝔞)g_±(N𝔞) = (x/h)(1 + O(x^{−c})), which sandwiches ψ_C(x) − ψ_C(y), where
ψ_C(u) := Σ_{N𝔞≤u, (𝔞,f)=1, [𝔞]=C} Λ_K(𝔞). By the row count of Lemma 1 without primitivity, the
ideals of norm ≤ u in C number ≤ (2V+1)(2(u/a)^{1/2}+1) ≤ 8u|D|^{−1/2} + 6u^{1/2} + 1. So
ψ_C(y) ≤ C log x(x^{1−η}|D|^{−1/2} + x^{1/2}) ≤ (x/h)x^{−η/2}, using h ≤ 3|D|^{1/2}log²(3|D|) and
|D|^{1/2} ≤ x^{0.286}. Partial summation, and removal of prime powers (≤ 6√x), give π_C(x). □

**Corollary 4.** Assume (L).
- The least prime of the form x² + ny² is ≪_ε n^{2B♭+ε}. Here D = −4n; trivially such a prime is ≥ n.
- Let L/Q be generalized dihedral with imaginary quadratic subfield K, and let L lie in the ring class
  field of conductor f. Then for every σ ∈ Gal(L/K), the least p with Frob_p = {σ, σ^{−1}} is
  ≪_ε (f²|d_K|)^{2B♭+ε}. Lift σ to the ring class field and apply Theorem 2.

**Remarks.**
1. The exponent is 2B + ε (the ε is the |D|^ε of Theorem 1(a)), so every improvement of the
   half-plane lowers it. B = 1/2 would give 1 + ε, the GRH exponent up to ε. For Linnik's problem, by contrast, the density peak at σ = 5/7 makes B irrelevant.
2. Simple variant: take X = |D|^{1/2}, no k-th powers, and convexity instead of Lemma 4. This gives
   N_D ≪ |D|^{2(1−σ)+ε} unconditionally and p(Q) ≪ |D|^{2+ε} under (L).
3. With (L) alone and the trivial zero count, x^{δ♭} ≫ h log|D| gives only |D|^{1/(2δ♭)} = |D|^{3.996}.
   The same comes from the least Chebotarev prime (FSW Cor 1.2, Hu (9.10)) on the ring class field.
   The best published family densities do better but stay above 2B♭: (L) with Pascadi–Thorner
   (2508.14888, Thm 1.3) gives |D|^{3/(3−2B♭)+ε} = |D|^{2.40} for fundamental D (up to |D|^{4.0} for
   D = f²d_K with f large), and (L) with Kowalski–Michel (Pacific J. Math. 207 (2002), Thm 2) gives
   |D|^{11}. Neither is stated in those papers.
4. Three structural facts enter: h ≈ |D|^{1/2} (the family is the square root of its conductor);
   χ(m𝔞) = χ(𝔞) (Lemma 2); and χ̄(𝔞) = χ(𝔞̄), which makes the dual AFE sum a ring class sum (Lemma 4).

**Records for candidate 1.**

| Bound on p(Q) | Source |
|---|---|
| \|D\|^455 | Zaman, PhD thesis (Toronto 2017), Cor 1.4.1: best unconditional; Gaudet (arXiv 2408.01610) calls it the record |
| \|D\|^694 | Thorner–Zaman, ANT 11 (2017) (arXiv 1604.01750), Thm 1.2: best in a journal, unconditional |
| \|D\|^{20/3+ε} on average over D; \|D\|^{3+ε} for most forms of most D | Ditchen, arXiv 1312.1502 |
| \|D\|^L, L unspecified | Gaudet, arXiv 2408.01610 (prime D, sieve; L = 1.3·10^10 in his 2023 thesis); Hexagon 2610.00178 §54 (Oct 2026, uses (O)) |
| \|D\|^{3.996} | substitution of (L), or of FSW/Hu with (L); not stated anywhere found |
| \|D\|^{2.40} (fundamental D) | (L) with Pascadi–Thorner's family density (Remark 3); not stated anywhere |
| **\|D\|^{1.7498}** | **this note, Thm 2, under (L)** |
| \|D\|log⁴\|D\| | GRH (Lagarias–Odlyzko/Bach–Sorenson on the ring class field; Sardari arXiv 1802.06193) |

## 2. Ray classes over imaginary quadratic fields (candidate 2)

Setting: K imaginary quadratic, 𝔪 ⊂ O_K, G_𝔪 := Cl_𝔪(K) = I(𝔪)/P_{𝔪,1}, H := |G_𝔪|,
A := |d_K|N𝔪. Then H = h_KΦ(𝔪)/[O_K^*:O_{𝔪,1}^*] ≤ h_KN𝔪 ≤ 2|d_K|^{1/2}log(3|d_K|)N𝔪. For
χ ∈ Ĝ_𝔪: A(χ*) ≤ A, and L(s,χ), E_χ and N_𝔪(σ;t) are defined as in §0, with 𝔪 in place of f.

### Lemma 6 (ray class lattice count and large sieve)
For C ∈ G_𝔪 and N ≥ 1, #{𝔞 ∈ C integral : N𝔞 ≤ N} ≤ 8N/(N𝔪|d_K|^{1/2}) + 4.31(N/N𝔪)^{1/2} + 1.
Hence Σ_{χ∈Ĝ_𝔪}|Σ_{N𝔞≤N} c_𝔞χ(𝔞)|² ≤ 21log(3A)(N + A)Σ|c_𝔞|².

*Proof.* Fix an integral 𝔟 ∈ C^{−1} prime to 𝔪. Every integral 𝔞 ∈ C satisfies 𝔞𝔟 = (α) with
α ∈ 𝔟 and α ≡ 1 mod 𝔪, i.e. α in a coset α₀ + 𝔟𝔪 (CRT), with N(α) = N𝔞N𝔟. Embed K ⊂ C. The lattice
Λ = 𝔟𝔪 has covolume N(𝔟𝔪)|d_K|^{1/2}/2 and minimum λ₁ ≥ N(𝔟𝔪)^{1/2}. Take a reduced basis v₁, v₂
(|v₁| = λ₁). The coset points in the disc of radius R = (N·N𝔟)^{1/2} lie on
≤ 2Rλ₁/covol + 1 lines parallel to v₁, with ≤ 2R/λ₁ + 1 points on each. Hermite's bound
λ₁² ≤ (2/√3)covol gives the first claim. For the second: H times the first bound is
≤ 2log(3|d_K|)(8N + 4.31(|d_K|N N𝔪)^{1/2} + |d_K|^{1/2}N𝔪) ≤ 21log(3A)(N + A); then use
orthogonality and Cauchy–Schwarz as in Lemma 2. □

### Theorem 7 (zero density for the ray class family; unconditional)
N_𝔪(σ;t) ≤ C(ε)A^{3(1−σ)+ε}T³ for 0 ≤ σ ≤ 1.

*Proof.* For σ ≤ 0.6 use the trivial bound H·c₀log(AT²) ≤ CA^{1+ε}log T ≤ A^{3(1−σ)+ε}log T. For
σ ≥ 0.6, use Lemma 5 verbatim (f replaced by 𝔪) with X := A, Y := A^{3/2}, Z := 2Y log²Y.
- Class II: by (F5) and absolute values, |L(1/2+iv,χ*)| ≤ 2Σ_n d(n)n^{−1/2}|V(n/√A(χ*))| ≤
  C(ε)(A^{1/2}(|v|+3))^{1/2+ε} (convexity), and |E_χ| ≤ C(ε)A^ε. Lemma 6 gives
  Σ_χ|M_X(1/2+iv,χ)|² ≤ 21log(3A)(X + A)(1 + log X)². So |𝒳_II| ≤ C(ε)Y^{1−2σ}A^{3/2+ε}T².
- Class I: all blocks have 2^{j+1} > X = A, so take k = 1 and M = 2^{j+1} ∈ (A, 4Z]. Lemma 6 has
  constant ≤ 63log(3A)M, the coefficient sum is ≤ C(ε)M^{1−2σ+ε}, and Sobolev is as in Theorem 1.
  So |𝒳_I| ≤ C(ε)Z^{2−2σ+2ε}.
With Y = A^{3/2} both bounds are ≤ C(ε)A^{3(1−σ)+3ε}T². □

### Theorem 8 (least prime ideal in a ray class)
Assume (L). For K imaginary quadratic, 𝔪 ⊂ O_K and C ∈ Cl_𝔪(K), there is a prime ideal 𝔭 ∈ C of
degree one with N𝔭 ≤ C(ε)(|d_K|N𝔪)^{3+ε}.

*Proof.* As for Theorem 2, with x = A^{3+ε}y. For 0.6 ≤ σ ≤ B♭, Theorem 7 gives
x^σN_𝔪(σ;k) ≤ x(A^{−ε}/y)^{δ♭}A^{ε′}T_k³. For σ ≤ 0.6, x^{0.6}A^{1+ε′} ≤ xA^{−0.19}y^{−0.4}; keep y^{−0.4} to absorb the log x in (2), as in
Theorem 2. The error
H log(Ax) is ≤ xA^{−1}. Prime powers, degree-two primes and ramified primes contribute
≤ (x/H)A^{−0.4}. □

Comparison:
- Thorner–Zaman 2017, Thm 1.1 (with F = K, L = the ray class field, A = Gal(L/K)) gives
  ≪ |d_K|^694 N𝔪^521.
- Substituting (L) into the least Chebotarev prime (ArtinHalfPlane Thm C) gives
  ≪ (C h_𝔪 log(3A))^{1/δ♭} ≈ (|d_K|^{1/2}N𝔪)^8.
- GRH gives ≪ (H log A)².
- Example: for K = Q(i) and 𝔪 = (q), the least Gaussian prime ≡ α (mod q) has norm ≪ q^{6+ε}.
- For general K only the substitution is available: the lattice count needs units and higher-rank
  geometry.

The ring class exponent 2B is not reached here. The large sieve loses |d_K|^{1/2} (thin coset
lattices), and the decomposition 𝔞 = m𝔟 fails because χ((m)) ≠ 1.

## 3. Goldbach exceptional set (candidate 3): known

- Pintz–Althöfer (arXiv 2610.15559, Thm 1): H_{θ₀} for some θ₀ < 1 gives E(X) ≪ X^{7/12+ε}. Here
  7/12 = 2·19/24 − 1, and 19/24 is the Maynard–Pandey–Radziwiłł minor-arc exponent. PA themselves
  note that improving 7/8 does not change 7/12 (§3), and that Guth–Maynard does not help their
  Thm 3 (§5).
- Li–Liu (arXiv 2610.26756, 19 Oct 2026, abstract): if every nontrivial zero of every primitive
  Dirichlet L-function has real part ≤ Θ < 11/12, then E₊(N), E₋(N) ≪ N^{1/2}(log N)^100, via a
  fourth-moment minor-arc bound. With (O) this holds unconditionally. N^{1/2} is the GRH
  (Hardy–Littlewood) exponent; no 2026 input goes below it.
- Checked: treating the arcs P < q ≤ X/P through L-functions, using (O) and the
  Montgomery/Huxley/CGL densities, gives S(a/q) ≪ X^{0.807} at q ≈ X^{5/9}, which is worse than 19/24.

## 4. Least quadratic non-residue, least prime residue (candidate 4): substitution

- Hexagon 2610.00178, Lemma 1.6, from (O): every nonprincipal χ mod q has a prime
  p ≤ C₀(log 3q)^8 with χ(p) ≠ 1. The same argument, after subtracting prime squares (≤ 3√x), gives
  a prime with χ(p) = 1, ≤ C(log 3q)^8. (L) contains the Dirichlet case, so 8 can be replaced by
  1/δ♭ = 7.9928: a substitution.
- Barrier: if χ(p) = 1 for all p ≤ x, then xW(1) = −Σ_ρ x^ρW(ρ) + O(log q). The right side is only
  ≤ x^{1−δ}·#{zeros near height 0} ≍ x^{1−δ}log q. Since W(7/8) ≥ W(1) for w ≥ 0, no weight helps.
  Going below 1/δ needs zero density for one L-function on the log q scale (the same dead end as the
  least Chebotarev prime, ArtinHalfPlane §6). Vinogradov's smooth-number trick needs bounds for
  L(s,χ) itself near Re s = 1 − δ; the half-plane gives only exp(C(log q)^{(1−σ)/δ}): no gain.

## 5. Short intervals, almost all intervals, BV beyond 1/2 (candidate 5): no

- A half-plane at B only removes density estimates for σ > B. The maximum of the best density
  exponent, 30/13, sits at σ = 7/10 < 7/8: Ingham's 3/(2−σ) and Guth–Maynard's 15/(3+5σ) meet there.
  So h = x^{17/30+ε} (all intervals) and h = x^{2/15+ε} (almost all intervals) are unchanged.
- BV beyond 1/2: zero-based methods reach at most Q ≤ x^{1/2}, the GRH range. Dispersion and
  Kloosterman methods do not use zero-free regions.
- PNT in progressions for q ≤ x^{3/7−ε}, and its short-interval forms, are CGL
  Cor 1.4 with (O) (ArtinHalfPlane D′). The x²(log x)^{−A} term in Barban–Davenport–Halberstam and
  Montgomery–Hooley becomes x^{2−c}: routine.

## 6. Class numbers (candidate 6): no

- The project's Lean bound Re L(1,χ) ≥ e^{−101}/log q (FORMALIZATIONS.md) is far from the printed
  state of the art:
  - Hexagon 2610.00178, Lemma 19.2 (Tatuzawa + (O)): h(D) > 0.655|D|^{1/2}/(πe log|D|) for
    |D| > e^{11.2}, explicit.
  - Hexagon 2610.00178, Corollary 1.5: h(D) ≥ c|D|^{1/2}/log log|D|, with c not explicit.
  - OpenAI QRH (5 Oct) announces h(D) ≫ √|D|/log log|D| and completes Euler's 65 idoneal numbers.
- An explicit loglog constant is bookkeeping: explicit zero counts, explicit Mertens, and the smooth
  explicit formula at x = (C log q)^8. It changes no order of magnitude. Using it for the class
  number problem beyond Watkins' h ≤ 100 would need class numbers of all |D| up to ≈ (150H)²: out
  of scope (CPU).

## 7. Other constants (candidate 7)

### Theorem 9 (Piltz divisor problem)
Assume (L), for ζ only. For every k ≥ 2, Δ_k(x) := Σ_{n≤x}d_k(n) − xP_{k−1}(log x) ≪_{k,ε} x^{B♭+ε}.
So α_k ≤ B♭ = 0.8749 and β_k ≤ B♭. (O) gives 7/8.

*Proof.* Lindelöf holds right of the half-plane: ζ(σ+it) ≪_{η,ε} |t|^ε for σ ≥ B♭ + η, |t| ≥ 2. This
is Titchmarsh's §14.2 argument; the Lean statement `zeta_lindelof` in
Riemann/Deep/QRHCor/Lindelof.lean covers any zero-free half-plane. Apply truncated Perron
(Titchmarsh Lemma 3.12) with c = 1 + 1/log x and T = x:
D_k(x) = (2πi)^{−1}∫_{c−iT}^{c+iT}ζ(s)^kx^s ds/s + O_k(x^{1+ε}/T). Move the line to Re s = B♭ + η.
The residue at s = 1 is xP_{k−1}(log x); the horizontal sides are ≪ x^{1+ε}T^{−1+kε}; the vertical
side is ≪ x^{B♭+η}T^{kε}log T. Take ε = η/k. □

Comparison with the Analytic Number Theory Exponent Database (teorth/expdb, `divisor_sum.tex`,
updated July 2026). For large k its best bounds are:
- Bellotti 2023: α_k ≤ 1 − 1.421(k−1.18)^{−2/3} for k ≥ 30, and 1 − 1.889k^{−2/3} for sufficiently
  large k;
- Ivić's table: (k−1)/(k+4) for 26 ≤ k ≤ 50, then (31k−98)/(32k), then (7k−34)/(7k);
- Heath-Brown: 1 − 0.849k^{−2/3}.

B♭ beats every one of them for k ≥ 59. For 40 ≤ k ≤ 58 it beats every bound stated with an explicit
range: Bellotti's k ≥ 30 bound is 0.8761 at k = 40. Under the Lindelöf hypothesis, α_k ≤ 1/2. This
is a two-line corollary of Lindelöf right of B♭, which the project already has in Lean at 7/8. It is
recorded only because it improves a tabulated classical exponent. Not found in Hexagon 2610.00178,
the OpenAI preprints, or arXiv 2026 (abs:"divisor problem", abs:"Piltz").

### Leads, not pursued
- (Hu) covers Hecke characters of infinite order. The record for Gaussian primes in narrow sectors is
  angle X^{−0.381} (Harman–Lewis 2001, sieve). The ring class method applied to the angular family
  {λ^m : |m| ≤ M} (size M, conductor ≍ M², angular large sieve N + M) gives only angle
  x^{−1/(4B)} = x^{−0.286}, which is worse. A Huxley-type density for {λ^m} plus (Hu) might give
  x^{−5/12}; this needs a literature check (Coleman 1990, Ricci).
- Indefinite forms (real quadratic orders): the decomposition 𝔞 = m𝔟 still works for narrow ring
  class characters (m is totally positive). The lattice count needs a unit-domain argument.

## 8. Novelty audit (2026-10-20)

- **Records found.**
  - Thorner–Zaman, arXiv 1604.01750v2, Thm 1.2: |D|^694.
  - Ditchen, arXiv 1312.1502: averages.
  - Zaman, arXiv 1710.08914: cites "O(D^700)".
  - Thorner–Zaman, arXiv 1803.02823, Thms 1.4 and 1.6: relative error exp(−c log x/log|D|), so an
    asymptotic only for |D| = x^{o(1)}.
  - Gaudet, arXiv 2408.01610: sieve, unspecified L.
  - Sardari, arXiv 1802.06193 and 1803.03218: GRH, h²log², and √D log^{2+ε} for almost all forms.
  - Hexagon 2610.00178 §54: |D|^A with A unspecified.
- **Sources read.**
  - Hexagon 2610.00178: all 55 section titles; §§1, 19, 54 read.
  - FSW arXiv 2610.18491v2, Hu Hexagon 2610.00216, Liu 2610.26835 (tex): none mentions quadratic
    forms, ring classes or divisor problems.
  - OpenAI QRH, 30 Sep and 5 Oct: class numbers and idoneal numbers only.
  - OpenAI preprint list: 749 titles, keyword-filtered; none on least primes, forms or Piltz.
  - Pintz–Althöfer 2610.15559.
- **arXiv API searches.**
  - abs:"least prime" AND abs:"binary quadratic form": 3 hits, all ≤ 2016.
  - abs:"quadratic form" AND abs:"least prime": 3.
  - abs:"prime represented" AND abs:"quadratic form": 9, all ≤ 2021.
  - abs:"binary quadratic" AND abs:"Linnik": 3.
  - abs:"ring class" AND abs:"least prime": 0.
  - abs:"class group" AND abs:"zero density": 9, none with small exponents.
  - abs:"class group L-functions": 5.
  - abs:"quasi-Riemann": 12 since May 2026. This found Li–Liu; none is on forms, ray classes or
    divisor problems.
  - abs:"divisor problem", abs:"Piltz": nothing using QRH.
- **Other searches.**
  - Hexagon site search for quadratic form, least prime, binary quadratic, class group, Linnik, ring
    class, represented: only Naslund's 3.99.
  - GXL evidence API, two claims: only the records above.
  - Keyless web search.
- **Second audit (2026-10-20): [NOVELTY.md](NOVELTY.md).** The closest densities (Pascadi–Thorner
  2508.14888, Kowalski–Michel 2002, Weiss/Thorner–Zaman) are weaker in |D| for every σ. The
  unconditional record for p(Q) is |D|^455 (Zaman's thesis), not |D|^694.

## 9. Dead ends (NOGO candidates)

- Ingham/Huxley for the ring class family: the fourth moment via the large sieve is |D|, the size of
  the conductor, not |D|^{1/2}. Halász–Montgomery off-diagonals are N^{1/2}(|D|T²)^{1/4}. Neither
  improves 2σ.
- Ray classes: the coset large sieve costs |d_K|^{1/2}, so the exponent is 3, not 2B.
- Goldbach below N^{1/2}: GRH barrier. 7/12 is superseded by Li–Liu.
- Least QNR below (log q)^{1/δ}: needs density on the log q scale.
- Short intervals: the density maximum sits at σ = 7/10, below every 2026 half-plane.
