/-!
# Round 124 — quick audits of 8 more new arXiv math.NT papers (scouts at concurrency 4; flaws re-checked in the TeX)

* 2610.23280 Yang, largest prime factors of consecutive integers: STATEMENT ≠ PROOF.
  - Theorem 1.3 claims positive density for P⁺(n−1) < P⁺(n) < P⁺(n+1) > max_i P⁺(n+a_i).
  - The proof ((3.6)–(3.7)) gives P⁺(n−1) < P⁺(n+1) < P⁺(n) > max_i …, a different pattern; the second
    pattern also does not match.
  - The stated patterns are unproved as written; a "similar" argument might fix it. Minor or moderate.
* 2610.19574 Li, integers not a prime plus an odd Zumkeller number (≪ x/(log log x)^{0.2378}): the chain
  log z_r ≪ log x ⇒ exp(D·2^r) ≪ log x ⇒ 2^r ≪ log log x runs the wrong way. It only has an upper bound
  log z_r ≪ exp(D·2^r), and it then needs a LOWER bound on 2^r. Repair: take r maximal with
  exp(D·2^r) ≤ log x/(4L). Minor.
* 2610.27591 small values of ζ near the critical line (RH; ChatGPT/Claude-assisted): the scout's "false
  equation" is REJECTED. x₀ = exp(α^{ν−1}), so log log x₀ = (1−ν) log(1/α) exactly.
* 2610.24900 Grobner: totally real abelian fields with h = 1 are Euclidean. It is conditional on "Scholium I", a
  uniform zero-free strip for finite-order Hecke L over cyclotomic fields (quasi-RH type; cf. Liu 2610.26835,
  under Orca audit). Unconditional only if that holds.
* 2610.19534 (partitions into distinct primes; GPT-assisted), 2610.28001 Diao (S⁺(X) average class numbers),
  2610.26705 (ζ at zeros of real-rooted trigonometric polynomials), 2610.22747 (Hua, prime squares):
  nothing found; small cases and constants rechecked.
-/

namespace Fayed.R124
end Fayed.R124
