/-!
# Round 131 — explicit Siegel-zero constants, Linnik–Goldbach and quasi-RH refinements (8 scouts, concurrency 4)

* 2610.15784 Li–Liu, 6 powers of 2 suffice in Linnik–Goldbach (5 under GRH). The scout's flaw ("z = 1 lies in
  D_{5/6}, so |Q|² < 2/15 fails") is REJECTED: β_λ = (r + √(r²+4r))/2 with r = 25/72 is 0.788, not 1, so
  z = 1 ∉ D_{5/6}. The unconditional input (Pintz–Ruzsa exceptional moduli) is used in its stated form, and the
  constant identity 4·(121/256) − 2·(121/256 + 15/32) = 1/128 checks. No flaw at triage depth.
* 2610.15483 Banks, an explicit version of OpenAI's no-Siegel-zero proof:
  - (1−β) log q ≥ 3.0·10⁻³ for q ≥ 10^10, and with Lu–Zaman–Zhao for every q ≥ 3;
  - recomputed: δ₀ = 81/(5888 e^{3/2}) = 0.0030696; no flaw.
* 2404.06482 Thorner (joint Sato–Tate, exceptional Rankin–Selberg zeros): no ideal-class or uniform-in-D
  statements, so our Sato–Tate-in-every-class novelty is unaffected.
* 2610.18491 Frei et al.: Chebotarev under QRH(δ), small splitting primes, ℓ-torsion; one displayed integral
  has the wrong range ([x, 2x] for [x/2, x]). Trivial.
* 2610.15291 Xu–Sang, a 0.8748867 half-plane, the same constant as Liu's B♭ = 0.874886647. It extends the
  imported fourth-moment lemma (OpenAI's Lemma 18.1) from κ ≥ 3/4 to κ ≥ 13/18 by "re-running" an induction
  that is never written out. Possible shared gap behind every improvement on 7/8. → queued for an Orca audit.
* 2610.21511 Kimura: not reported (scout cut off).
-/

namespace Fayed.R131
end Fayed.R131
