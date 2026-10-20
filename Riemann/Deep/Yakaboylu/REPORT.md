NO DISCOVERY: the RH claim fails at an intertwining/positivity step that just assumes RH, and the author's own latest version (v17) drops the claim (now only "W > 0 ⟺ RH"). Lean proofs below.

Paper: E. Yakaboylu, arXiv 2408.15135. 17 versions. Read: v17 (latest, 2026-08-18, "Nontrivial Riemann Zeros as Spectrum"), v11 (last version under the briefed title "On the Existence of the Hilbert–Pólya Hamiltonian"), v16 (last version that claims RH). Sources: `.scratch/yakaboylu/{hilbert_polya_17.tex, v11/, v16/}`.
Claim history (from arXiv abstracts): v1–v6 say "a step toward RH"; v7–v16 claim RH for simple zeros; v17 only claims the equivalence.

Setup (all versions): a non-normal operator R = −D − iμ(T) on L²[0,∞). Ψ_ρ are eigenvectors of R with eigenvalue i(1/2−ρ), Φ_ρ are eigenvectors of R† with eigenvalue −i(1/2−ρ̄), and ⟨Φ_ρ'|Ψ_ρ⟩ = δ (biorthogonal). The plan: find W ≥ 0 with W R = R† W, then h = W^{1/2} R W^{−1/2} is self-adjoint, so the zeros are real.

Where it fails:
1. v11, text before eq. (23), eq. (23) W := Σ_ρ |Φ_ρ⟩⟨Φ_ρ|, and eq. (24): the claim R† = W R W⁻¹. In fact W R Ψ_ρ = λ_ρ Φ_ρ but R† W Ψ_ρ = conj(λ_ρ) Φ_ρ. So the identity holds iff every eigenvalue is real, i.e. iff RH. The only reason given is that the eigenvalue sets of R and R† coincide (eq. (21)). That does not give the identity, because the matching between eigenvalues is a permutation ρ ↦ 1−ρ̄, not the identity map. So the step is circular. The symmetry of h fails the same way: h has orthonormal eigenvectors with eigenvalues i(1/2−ρ).
   Lean: `v11_intertwine_iff` (for any invertible P and any spectrum d: W R = Rᴴ W ↔ ∀k, d_k real) and `v11_lemma_false` (spectrum {i, −i}, which is closed under conjugation, breaks it).
2. v16, Lemma "Positive Semidefiniteness of W": it writes W = ω(T)^{−2} K(R) and then W = (K^{1/2})† ω^{−2} K^{1/2} ≥ 0, using a "conjugation-symmetric branch" of the square root. No such branch exists: at −1 it would have to be real and square to −1.
   Lean: `no_conj_symm_sqrt`, `v16_lemma_false` (V > 0 and V K = Kᴴ V do not force V K ⪰ 0), and `v16_sqrt_rewrite_false` (no matrix M has Mᴴ M = diag(1, −1)).
3. Side gap, used again in v17: v11 says the biorthogonal set "forms a complete basis" with Σ|Ψ⟩⟨Φ| = I on the closed span, and v17's corollary on non-orthogonal projections says the partial sums converge strongly to the identity on the closed span. That is false in general.
   Lean: `resolution_false` uses ℓ² with Ψ_n = e_n − e_{n+1} and Φ_n = e_0+…+e_n. These are biorthogonal and e_0 lies in the closed span of the Ψ_n, but the partial sums are e_0 − e_N, which is not a Cauchy sequence.
   This does not affect RH in v17, since v17 claims no RH. It does undercut v17's W⁻¹, the self-adjointness argument for h, and the claimed spectrum.

v17 status: Theorem 4.4 (W > 0 ⟺ Re ρ = 1/2) is essentially correct on finite combinations; its own proof shows that W has a negative direction whenever an off-line pair of zeros exists. That is exactly the gap behind item 1. No unconditional result touches RH.

Prior critiques: I searched the web, the arXiv API ("Yakaboylu") and Watkins' list of proposed RH proofs, and found no third-party critique. But the author's v17 retreat concedes that positivity, the step the proof depends on, is equivalent to RH. Steps of this kind (the conclusion built into the intertwiner) are routine, so I see nothing significant here.

Lean: `Riemann/Deep/Yakaboylu/Paper.lean`. `lake build Riemann.Deep.Yakaboylu.Paper` passes, with no sorry. #print axioms gives only propext, Classical.choice and Quot.sound for all six theorems. The file is not imported from `Riemann.lean`; that is the orchestrator's call.
