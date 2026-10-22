/-!
# Round 126 — quick audits of 8 September math.NT papers on zeros and explicit bounds (scouts at concurrency 4)

* 2609.00537 (explicit approximate functional equation for ζ): the statement of Theorem AFE1 sums over
  1 < n ≤ ct, but the proof and the corollary use 1 ≤ n.
  - As printed it is false: at σ = 1/2, t = 14.13473, c = 1 the scaled error is 3.765 > m(1) ≈ 1.228
    (mpmath); with n = 1 included it is 1.053.
  - Misprint.
* 2609.35950 (primes in arithmetic progressions and Siegel zeroes): sloppy and likely AI-written.
  - The displayed goal (l. 362–370) has the same sum on both sides.
  - D* is used but never defined.
  - The error term in Main Theorem (ii) lacks a factor x.
  - Not a mathematical result as written.
* 2609.02713 Grobner: the folklore formula π(x) = R(x) − Σ_ρ R(x^ρ) DIVERGES; its partial sums are not O(T^θ)
  for θ < Θ = sup Re ρ. This is plausible: R(x^ρ) only decays like |γ|^{Θ−1}, via the tail Σ_{n>N} μ(n)/n with
  N ≈ |γ| log x. It is consistent with von Mangoldt's (convergent) formula for J and ψ. No flaw found; not ours.
* No flaw found:
  - 2609.23222, explicit PNT bound for x ≥ 21 (small x tested);
  - 2609.24167 Wang, > 67.25% simple + on the line, > 83.625% distinct;
  - 2609.33043 Knausgård (rests on an external computer-assisted proposition);
  - 2609.07918, short-interval proportions;
  - 2609.14471 Gu–Nguyễn (GRH Sárközy, N^{89/90}).
-/

namespace Fayed.R126
end Fayed.R126
