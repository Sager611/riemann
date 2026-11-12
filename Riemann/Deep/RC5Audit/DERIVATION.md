# RC5Audit derivations (2026-11-12)

Referee of [RCApps5 NOTE](../RCApps5/NOTE.md) §2–5 (line numbers = that file). [R3] = [RCApps3 NOTE](../RCApps3/NOTE.md), [R] = [RCApps NOTE](../RCApps/NOTE.md). Scripts and outputs: [checks/](checks/).

## D1. Lemma P (l.50–77): re-derived, correct

- (a) p | f, p | nn′ ⇒ p odd (n odd) ⇒ p² | Δ ⇒ (Δ/p) = 0, so only c with (c, f) = 1 survive. T² − Δy² = 4M ⇒ T² ≡ 4M (p²); p ∥ M ⇒ p | T ⇒ p² | 4M: impossible. v_p(nn′/c²) = v_p(n) + v_p(n′) ∈ {0,1,2}, so p | n ⟺ p | n′ and e := (n,f) = (n′,f) (squarefree). With n = em, n′ = em′: (mm′, f) = 1, (n,n′) = e(m,m′), (c,e) = 1 ⇒ c | (m,m′), (Δ/c) = (e²/c)(Δ_e/c) = (Δ_e/c).
- (b) e odd ⇒ e² ≡ 1 (8) ⇒ Δ_e ≡ Δ (8), conductor f/e. p | e: p² | Δy², p² | 4e²M′ ⇒ p | T ⇒ e | T. T/e ≡ T ≡ yΔ ≡ yΔ_e (2) ⇒ γ/e ∈ O_{Δ_e}; e·O_{Δ_e} ⊂ O_Δ. Bijection on elements of norm e²M′ ↔ M′: r_Δ(e²M′) = r_{Δ_e}(M′).
- (c) γ ∈ O_{Δ_e} = Z + f_eO_K ⇒ γ ≡ a (f_eO_K), Nγ ≡ a² (f_e); (Nγ, f_e) = 1 ⇒ (a, f_e) = 1 ⇒ γO_K ∈ P_{K,Z}(f_e), prime to f_e. Generators inside O_{Δ_e}: one O_{Δ_e}^×-orbit (u = a′/a mod f_eO_K ⇒ u ∈ O_{Δ_e}). So #{γ ∈ O_{Δ_e}: Nγ = M} = w_e[1](T_M). [R3] Lemma 1 for Δ_e uses only [(p)] = 1 (also |Δ_e| ≤ 4: G trivial). [1](SS*) = Σ_C S_C·conj(S_C) via 𝔞 ↦ 𝔞̄. Hence B_e = (w_e/2)Σ_C|S_C|² ≥ 0.
- Own check (`lemmaP_mine.py`, independent code): LHS by direct (T,y) counts, RHS from reduced primitive forms of each Δ/e² (r_F(m)/w_e ideals per class): 231 orders Δ ∈ [−2600, −7] with an odd prime ≤ 13 dividing f, random complex b on the 28 odd squarefree n ≤ 3000 built from 3–13; all non-maximal, all with a nonzero e > 1 block, 14 with Δ_e ∈ {−3, −4}. Max relative deviation 2.7·10⁻¹⁵; no negative block.

## D2. Corollary P (l.79–89): correct

- Identity: O_Δ = Z[(Δ+√Δ)/2], so (T + y√Δ)/2 ∈ O_Δ ⟺ T ≡ yΔ (2); summing r_Δ over Δ with Φ gives ½𝒩_c.
- Inequality: for D ∈ 𝔇(X) (fundamental, |D| > X/2 ≥ 16) only e = 1, w = 2: Σ_χ|Σbλ_χ|² = h(D)Q′_D(b). Φ(|D|/X) = 1 there, every other Φ(|Δ|/X)Q′_Δ ≥ 0, |Δ| ≤ 4 excluded by X ≥ 32. Cosmetic: [C] (F3) gives h ≤ 3|D|^{1/2}log²(3|D|) (one log more than written); absorbed in X^δ.

## D3. Theorem 1 (l.99–122): correct for the detector sequences

- 𝒬′ = ½Σ_c Q_c (n = cn₂, n′ = cn₂′). supp W_{j,t′} ⊂ [1/(4ℓ²), 4] = [1/16, 4], so n₂ ∈ [L/64, L], L = 4N/c. Q_c is [R3]'s Q_{c,1} (β = b(c·), X′ = X).
- (i) Lemma 6 at w = 1: (c, M) = 1 (n = cn₂ squarefree); M ≤ L² ≤ X²; c√M ≤ cL = 4N, 64c√MX^δ ≤ 256X^{1−3δ} ≤ X iff X^{3δ} ≥ 256. [R3]'s own check at general w (64c√MX^δ ≤ 256NX^δ/w ≤ X/w²) is the same inequality at w = 1.
- Lemma 6 re-derived: ∫f = (y²X′/(2√M))K(y²X′/(4M)); local sum 4A_M(y)μ(c) for (y,c) = 1 (Jacobsthal Σ_T((T²−4M)/p) = −1, p ∤ 4M), 0 for p | (y,c) (Σ_j((A₁ + 2T₁j)/p) = 0, p ∤ 2T₁); per y: 2·4A_Mμ(c)/(4y²c)·∫f = μ(c)X′A_M K/(c√M). Errors: ‖f^{(i)}‖₁ ≪ √M(4√M/(y²X′))^i, |f̂(k/q)| ≪ √M(16c√M/(|k|X′))^i ≤ √M(X^{−δ}/(4|k|))^i.
- (ii) M = n₂n₂′ odd, cube-free. The 2-part A^{(2)}_M(2^a) depends on T₀ mod 2^{2a+2} only (y₁² ≡ 1 (8) for odd y₁), so A_M = A^{(2)}·A^{odd} and §3.5's factorisation holds.
- (iii) v ≥ X/(4L²) = c²X/(64N²); S(u) = 0 for u ≥ √8 ⇒ 2^ate < √8·8N/(c√X) = 22.6N/(c√X) ≤ E₀, uniformly in M.
- (iv) E₀² ≤ L ⟺ N ≤ 4cX/529. κ | e, e ≤ E₀ ⇒ κE₀ ≤ L ⇒ E₀ + L/κ ≤ 2L/κ.
- §3.5: Š(z) = (1−2^{−z})ζ(z)ǧ(z), ǧ poles at 0, −2, … (g even); z = 0 cancelled; only z = 1 crossed on the way to Re z = −δ; residue gives 2c₁√X′𝔐.
- (S) re-checked case by case: h(p) = (n₂/p)(n₂′/p) − 1[p|n₂] − 1[p|n₂′] + (p+1)1[p|n₂]1[p|n₂′] gives (M/p), −1, p−1 for p ∤ M, p ∥ M, p² ∥ M; h(p²) matches. Notation clash in [R3] (not RCApps5): F_{κ,q,s} is used as F_{T₂T₄T₅,T₁,T₅}, i.e. its character is (n/T₁), while the triple uses q = T₁T₅; F(κm) = sign·(m/T₁T₅) as stated. Cosmetic.
- §3.6: coefficients γ_m = β(κm)(κm)^{(z−1)/2}ψ(κm)·sign·(m/T₅) are fixed as T₁ varies; Cauchy–Schwarz over T₁ then [HB95, Thm 1] (odd squarefree q, m). Bound uses no structure of b.
- §3.7 / Lemma 9: a_V(κm) = Σ_{k₁|κ}Σ_{k₂|m, k₁k₂≤V}μ(k₁)μ(k₂) ((κ,m) = 1); 𝒟(s′) = L(s′,ψ̃)/L(2s′,ψ̃²)·Π_{p|2Pκ}(1+ψ̃(p)p^{−s′})^{−1}·(double sum ≤ d(κ)V^{1/2}log V). A pole only for q = 1, ψ = ψ₀ (ψ mod 8 real, (·/q) non-principal for squarefree q > 1); residue ∝ Ŵ(1−s) = ∫W(v)v^{−σ−it′}dv = 0 by (Z). Twisted mean value Σ_{q≤E₀}q^{−1}|L(½+δ+iτ, ψ(·/q))|² ≪ E₀^δ(1+|τ|)^{1+δ}.
- §3.8: y = 0 ⇒ n₂ = n₂′, ≤ 2Σ_ΔΦ ≤ 4X per n₂.
- No step of §3.4–§3.8 uses w > 1 or coprimality to a conductor. Σ_c|Q_c| ≪ X^{12δ}(X + √XV)‖b‖²; × CX^{1/2}log X gives X^{3/2+13δ}N^{1−2σ}.
- For arbitrary b the statement is false (trivial and genus characters, §N2 below); RCApps5 restricts to b_{j,t′}, ∂_{t′}b_{j,t′} (l.100), and the discovery file says so.

## D4. Theorem 2 (l.126–141): correct

- V = X^{1/2}/2 < |D|^{1/2} (|D| > X/2), as [R3] Lemma 3(c) needs; 4Z = 8X^{1−5δ}log²Y ≤ X^{1−4δ}.
- Class II: 3/2 + 2δ + (1−5δ)(½+δ−σ) + δ = 2 − σ + 4δ + 5δ(σ−½−δ) ≤ 2 − σ + 6δ.
- Imported at Re s = ½+δ: Σ_{χ≠1}|L(½+δ+iv,χ)|² ≪ |D|^{1/2+δ}(|v|+3)^{1+δ}, cited as [D2] (2.1), which is stated on Re s = ½ (Density2 NOTE l.150). Same AFE proof at ½+δ. Routine.
- Class I: |τ| ≤ X^δ/2 and |t| ≤ X^δ/4 ⇒ |t′| ≤ X^δ (Theorem 1 range); Σ_j 2^{j(1−2σ)} ≤ (V/4)^{1−2σ}/(1−2^{1−2σ}) (σ ≥ ½+3δ, absorbed in C(δ)); X^{3/2}V^{1−2σ} = 2^{2σ−1}X^{2−σ}. Large |τ|: per-D sieve × X^{−3}.
- Detector ([R3] Lemma 4) re-checked: (Z) by u = ℓ²v; ∂_{t′}W = 2i log ℓ·ℓ^{2−2σ−2it′}e^{−ℓ²v2^j/Y}φ₂(ℓ²v), so (Z) holds for W^♯.

## D5. Theorem 3 (l.145–171): correct

- 1′: δ = (A−1)/8000 ≤ A/4400 and ≤ (A−1)/80.
- 2′: on [σ₀, 9/10], φ = 1 − (A−1)(1−σ) ≤ 1 − (A−1)/10 = 1 − 8ε′; on (9/10, 1), F₂′ = −10/(3−2σ)² < 0, F₂(9/10) = 5/6 ≤ 1 − 8ε′ ⟺ A ≤ 8/3. Theorem 2's window condition: T_X + 3 = 2X^{0.26δ} + 3 ≤ X^δ/4. Markov: X^{φ+2ε′+61δ+2(κ+ε₁)} ≤ X^{1−6ε′+62δ} = X^{1−5.38ε′} (the extra 2(κ+ε₁) is slack).
- 3: 4Aδ + κ + ε₁ ≤ 5.06δ ≤ 6δ; (A−1)/2 − 6δ = 3994δ > δ/2.
- 4 and [R] §2.4(a), (d), §2.5: ε′ = 100δ > δ/2; 0.005A − κ − ε₁ > δ/2; x ≥ |D|^{1+2δ} since A ≥ 1 + 8000δ; c = δ/2 ≤ min{ε′, 0.505δ}.
- Effective but huge X₀: Theorem 1 needs X^{3δ} ≥ 256 (and more), so X₀ ≥ 2^{8/(3δ)} = 2^{64000/(3(A−1))} (A = 1.1: 2^{213333}). Not a flaw.

## D6. §6 (D > 0, sectors): reduction checked, [R4] inputs not re-refereed

- Lemma P⁺: γ ↦ γ/e keeps γ ≫ 0 and u_γ; each block is [R4]'s nonnegative Q for the non-square order Δ_e. Theorem 1⁺: X ≥ 128c√MX^δ as c√M ≤ 4N.
- Lemma P^∠: ℓ even (γ ↦ −γ invariant); for Δ_e ∈ {−3,−4}, Σ_u u^ℓ = w_e1[w_e | ℓ]. Theorem 1^∠: 64c(1+2ℓ²)√MX^δ ≤ X for |ℓ| ≤ X^δ, N ≤ X^{1−4δ}.

## N1. Detector sequences, brute force (`fam.c` mode 0, `run_det.sh`, `an_det.py`; `det_*.out`, `det_table.txt`)

Family 𝔇(X) exactly (no Φ, no non-maximal orders): classes = reduced forms (Σh checked against PARI qfbclassno at X = 2000: 5444), v_C(n) = r_F(n)/2, b on admissible n (prime factors ≥ 19), Y = X, φ₂ flat on [½, 2]. R1 := Σ_DΣ_χ|Σbλ_χ|²/(Σ_D h(D)·Σ|b|²).

| sequence | R1 range (X = 2000–32000, N = 2^j ≤ ~X) | slope in N on [√X, X^0.95] (X = 2000 → 32000) | slope in X at N/X = 1/64, 1/16, 1/4, 1 |
|---|---|---|---|
| σ = 5/6, V = √X/2, killed | 0.81–1.91 (> 1.1 only at N ≲ 2√X, where χ = 1 carries 18–56 %) | −0.21 → −0.09 | −0.23, +0.04, −0.02, +0.05 |
| σ = 0.55, V = √X/2, killed | 0.78–1.76 | −0.14 → −0.05 | −0.24, +0.05, −0.01, +0.05 |
| σ = 0.55, V = 2, killed | 0.76–1.61 | −0.11 → −0.01 | −0.19, 0.00, 0.00, +0.01 |
| σ = 5/6, t′ = 5, killed | 1.00–1.06 | −0.01 → 0.00 | 0.00, −0.02, 0.00, 0.00 |
| ∂_{t′}b, σ = 5/6, killed | 0.79–1.78 | −0.15 → −0.06 | −0.24, +0.05, −0.01, +0.05 |
| σ = 5/6, pole (no ℓ-term) | 1.35–53.8 | +0.67 → +0.72 | +0.16, +0.08, +0.28, +0.35 |
| σ = 0.9, V = 2, pole | 1.31–88.2 | +0.68 → +0.83 | +0.25, +0.42, +0.49, +0.51 |

Trivial-character share of the killed sums 4–15 % (N ≥ 4√X), of the pole sums up to 99.4 %. X = 64000 (`det64_64000.out`, N = 512 … 32768): every killed sequence R1 = 0.95–1.10; pole 2.05, 2.99, 11.6, 45.7. t′ scan (X = 16000, N = 2048, 8192; t′ ∈ {0.5, …, 40} incl. 14.1347, 21.022, 25.011; σ = 0.51): max R1 = 1.045 (`det_tp_16000.out`).

## N1′. The dual form 𝒬′ of Corollary P (`fam.c` mode 2, `run_dual.sh`; `dual_*.out`, `dual_table.txt`)

All negative discriminants with |Δ| ∈ [X/4, 2X], Φ flat on [½, 1] with smooth steps on [¼, ½] and [1, 2]. Per Δ: conductor f, squarefree e | f with prime factors ≥ 19 (other e give b_e = 0), primitive reduced forms of Δ_e = Δ/e², S_C = Σ_{m ≤ 4N/e, (m,f)=1} b(em)r_F(m)/w_e, block (w_e/2)Σ_C|S_C|². 𝒬′(b)/(X·Σ|b|²):

| sequence | X = 2000 | 4000 | 8000 | 16000 | 32000 |
|---|---|---|---|---|---|
| killed, σ = 5/6 (N = 2.1√X … X) | 0.47–0.58 | 0.49–0.59 | 0.51–0.62 | 0.53–0.62 | 0.54–0.61 |
| killed, σ = 0.55 | 0.45–0.56 | 0.48–0.60 | 0.50–0.62 | 0.52–0.63 | 0.53–0.62 |
| killed, ∂_{t′}b | 0.46–0.57 | 0.48–0.60 | 0.51–0.62 | 0.52–0.63 | 0.54–0.62 |
| pole, σ = 5/6, at N ≈ X | 16.3 | 20.6 | 27.0 | 34.7 | 41.6 |

At N = 1.4√X (X = 2000) and 2.0√X (X = 4000) the killed ratio is 0.67–1.21. Agrees with RCApps5 famq.c (0.50–0.53 at X = 20000, 40000) under a different Φ and support.

## N2. Large-sieve constants (`fam.c` mode 1, `eigs.py`, `run_eig.sh`; `eig_*.out`, `eig_fit.txt`)

Top eigenvalue λ₁ of the Gram operator Σ_DΣ_χ λ_χλ_χ^* on admissible n ∈ [N, 2N), scipy eigsh, in units of the family size Σ_D h(D) ≈ X^{3/2}; "notriv" drops χ = 1, "nongenus" drops all χ² = 1 (2^{t−1} genus characters per D from assigned characters; orthogonality checked for every D). Validation: Σ_C v_C(n) = Π_{p|n}(1 + (D/p)) on 2000 random (D, n) per run: 0 failures.

| X | N | full | notriv | nongenus | max diagonal |
|---|---|---|---|---|---|
| 2000 | 1024 | 10.62 | 3.30 | 1.38 | 1.49 |
| 4000 | 2048 | 15.05 | 4.50 | 1.37 | 1.56 |
| 8000 | 4096 | 21.14 | 5.92 | 1.32 | 1.37 |
| 16000 | 512 / 2048 / 8192 | 2.76 / 8.43 / 29.65 | 1.41 / 2.92 / 8.11 | 1.00 / 1.06 / 1.31 | 1.13 / 1.16 / 1.33 |

- full: slope in N 0.65 → 0.85 (X = 2000 → 16000), in X at fixed N/X = ¼, ½: 0.49, 0.49 ⇒ λ₁ ≍ X·N: the trivial character (pole). notriv: slopes 0.37 → 0.63 in N and 0.19–0.43 in X (N/X = 1/16 … 1/2): genus characters (b = (·/q) on D with q* | D).
- nongenus: slope in N 0.14 → 0.10 (tracks the max diagonal, a divisor effect), in X at fixed N/X: +0.04, +0.01, +0.02, −0.03; λ₁ = 0.88–0.99 × max diagonal.
- Projecting out ψ(n)(n/q), q ≤ N/40 odd squarefree, ψ mod 8 (the directions of [R3]'s main term 𝔐): full drops to 1.05–2.04, close to the nongenus level.
- Max memory 259 MB (X = 16000, N = 8192), below the 300 MB cap.
