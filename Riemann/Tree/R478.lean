/-!
# Round 478 — node N7 (zero-free regions and zero-density constants), round 18, falsification of our own fresh correction of Bellotti's explicit log-free zero density (bellotti2, 14 Nov): the (u+w)⁴ → (v+w)⁴ slip in (3.10), the code's d₃ = TeX d₃/log²T₀, the (4.9) repair with F_A ≤ 1.1362, the one-sided near count (c₃ doubles, c₄'s factor 2.4 covers both sides) and the certified |G| bound all HOLD in the TeX and numerically; the corrected Table 2 stands; attractive coupling of exclusion/zero-range processes, Perron–Frobenius for sofic shifts and Gelfand–Tsetlin interlacing transfer no tool

* Node: N7 (18th round; R469 before). Round type: falsification of our own fresh result (bellotti2 corrected a refereed paper and
  withdrew our older 'factor 2 in c₄' objection, lesson 267), aimed at the parts no reviewer had checked; sources: Bellotti's TeX
  (arXiv 2405.12545v2, final.tex) and Riemann/Deep/BellottiFix/NOTE.md.
* [falsify] (3.10) / (estimatesx): TeX l.532–538 bound |M(1)| ≤ d₂ log⁴(VW) with VW = T^{v+w} but print (u+w)⁴; l.503–505 fix
  0 < u < v, so the printed factor is the smaller (unsafe) one: our NOTE §3(a) HOLDS.
* [falsify] d₃: TeX (specialcase) l.464–487 gives d₃ = A/log T; the code's Psisquare (constants_orig.py l.95–97) divides by log³ t:
  code d₃ = TeX d₃/log²T₀ (1/825.39 at T₀ = 3·10¹²), the unsafe direction: NOTE §3(b) HOLDS.
* [falsify] (4.9), TeX l.830: X^y/y is minimal at y = 1/log X; the repair max_{[2ν, 2δ]} X^y/y ≤ F_A X^{2δ}/(2δ), F_A = e^{a−1}/a
  (a = 2ν log X < 1) holds on the grid a ∈ {0.3, 0.6, 0.9}, δ/ν ∈ {1.5, 3, 10}; the per-side count (H_{⌊1/δ⌋} + F_A − 1)/δ and
  c₂ → (5F_A − 4)c₂ re-derived; mono.log rows 26–38 F_A = 1.0011–1.1362: HOLDS.
* [falsify] (4.11) / (4.13): TeX l.843–849 take one representative per absolute-distance band, no doubling anywhere in (4.9)–(4.11)
  or boundc3, so c₃ doubles; (4.13) l.858–871 carries 2.4 ≥ 2 + 4(1 − α₀)e^{1/6−π/2}/I = 2.0587 (α₀ = 0.9927, I = 0.1222258)
  and 2.1199 (α₀ = 0.985, I = 0.1229359): HOLDS.
* [falsify] |G|: exact float64 sums A(W) = Σ_{j,k≤W} |θ_j||θ_k|[j,k]^{−1/2}: A(6836) = 638.7127 < 9.7572 W/log²W = 855.4802
  (ratio 0.7466), A(3000) = 332.3953 < 456.6415; max of S(y)/√y over squarefree-gap ends y ≤ 10⁶ is 2.428788 < 2.4297: HOLDS.
* [explore] exclusion / zero-range processes (attractiveness, basic coupling) → N8: ζ's zeros carry no order-preserving Markov
  dynamics: no transfer.
* [explore] sofic shifts (Perron–Frobenius for sofic presentations) → N1: finite-matrix dynamical zeta poles (golden-mean shift:
  |t| = 0.618, 1.618) put no constraint on Weil positivity: no transfer.
* [explore] Gelfand–Tsetlin patterns (Cauchy interlacing) → N13: needs a Hermitian principal-minor model of ξ's zeros, which is
  the Hilbert–Pólya operator itself: no transfer.
-/

namespace Fayed.R478
end Fayed.R478
