/-!
# Round 276 — Fasi-Higham off the diagonal: the first failing k is about 2m + 13 (data for m <= 49)

* Data (scouts; mpmath at 500 and 900 digits, several confirmed with exact rationals; n <= 400): the smallest k with a
  non-positive Taylor coefficient of 1/q_{k,m} is k0(m) = 66, 44, 42, 43, 44, 48, 50, 54, 57, 61, 65, 69, 73, 77, 81, 83,
  87, 91, 95, 99, 103, 111 for m = 5, 7, ..., 45, 49, and k0(m) - 1 is positive in range. For 23 <= m <= 33,
  k0 = 2m + 15; for 35 <= m <= 49, k0 = 2m + 13; the first failing n is about 3m + 10. So every k <= 2m + 12 stays
  positive in these ranges: the diagonal theorem (discoveries/fasi-higham-diagonal.md) sits well inside the positive
  region. Candidate theorem for a later attempt: positivity for all odd m and all k <= m (or k <= 2m).
* The failures themselves are known in kind (Scardino, 2610.39051, refutes the off-diagonal conjecture).
-/

namespace Fayed.R276
end Fayed.R276
