/-!
# Round 26 — flaw triage of recent claimed RH proofs (suspected, unverified)

| claim | mechanism | suspected failing step | class |
|-------|-----------|------------------------|-------|
| Christou Micheas, arXiv 2607.04338 | Hadamard factorization, zero pairing | pairing κ ↔ 1-κ is symmetry, not fixedness | symmetry ≠ fixedness |
| Covei, arXiv 2606.24924 | Riccati–Gamma concavity, symmetric cancellation | localisation of each zero assumed | symmetry ≠ fixedness |
| Gaber, arXiv 2607.26114 | ℓ-totient criteria | criterion stated, its hypothesis not proved | criterion ≠ proof |
| El Bachir, viXra 2505.0054 | real-series identity extended to complex s | convergence/continuation ignored | σ-blindness |
| Ip, Zenodo 15334802 | arithmetic self-adjoint operator | self-adjointness/completeness asserted | false self-adjointness |
| Gunes, JRSMAMS 2025 | eta main/remainder split, limiting ratio | uniformity of the limit missing | finite → limit |
| Geisner, Zenodo 17086205 (Lean) | Lean 4 development | five admit statements | unverified transitions |
| Tamburini, arXiv 2503.09644 | Majorana Hamiltonian in Rindler space | surjectivity onto all zeros not shown | spectral completeness gap |
| Desogus, arXiv 2609.20367 | Weil positivity on all windows | Gate II uniform positivity | under deep audit (Riemann/Deep/Desogus) |
-/
