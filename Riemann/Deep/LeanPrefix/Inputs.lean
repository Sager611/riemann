import Riemann.Deep.LeanPrefix.Lemma7A
import Riemann.Deep.LeanPrefix.Defs
import Mathlib.GroupTheory.SpecificGroups.Alternating
import Mathlib.GroupTheory.GroupAction.Primitive
import Mathlib.Data.ZMod.Basic

/-!
# The two literature inputs (named Props, passed as hypotheses of the main theorem)

* LZZW2026_Thm3_1: Lei, Zhao, Zhu, Wu, arXiv 2610.30892 (unrefereed), Theorem 3.1, first sentence.
* Jones2014_Thm1_2_k0: G. A. Jones, Bull. Aust. Math. Soc. 89 (2014), Theorem 1.2, case k = 0.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv
open scoped LinearAlgebra.Projectivization

/-- **Literature input** (unrefereed preprint by others). Yexiao Lei, Jinxing Zhao, Ziyu Zhu,
Yang Wu, *On the generation of the symmetric group by index relations of three prefix reversals*,
arXiv:2610.30892v1 (Oct 2026), Theorem 3.1 (TeX lines 221-223), first sentence, verbatim:
"If $m+k<n-1$, then $\left \langle r_{n},r_{m},r_{k} \right \rangle <S_{n}$."
Definition (TeX lines 87-95): $r_i(j)=i+1-j$ for $1\le j\le i$, $r_i(j)=j$ for $i<j\le n$, which
is rev n i on Fin n (0-based: x ↦ i-1-x for x < i). Standing assumption (TeX lines 109-112):
$G_{n}:=\langle r_n,r_m,r_k\rangle$, $2\le k<m<n$. "< S_n" is ≠ ⊤ (G is a subgroup of S_n).
The second sentence of the theorem is not used. -/
def LZZW2026_Thm3_1 : Prop :=
  ∀ n m k : ℕ, 2 ≤ k → k < m → m < n → m + k < n - 1 → G n m k ≠ ⊤

/-- **Literature input** (published). Gareth A. Jones, *Primitive permutation groups containing a
cycle*, Bull. Aust. Math. Soc. 89 (2014), 159-165, arXiv:1209.5169, Theorem 1.2 (TeX lines 54-72),
used for k = 0 only. Verbatim: "Let $G$ be a primitive permutation group of finite degree $n$, not
containing the alternating group $A_n$. Suppose that $G$ contains a cycle fixing $k$ points, where
$0\le k\le n-2$. Then one of the following holds: 1. $k=0$ and either (a) $C_p\le G\le AGL_1(p)$
with $n=p$ prime, or (b) $PGL_d(q)\le G\le P\Gamma L_d(q)$ with $n=(q^d-1)/(q-1)$ and $d\ge 2$
for some prime power $q$, or (c) $G=L_2(11)$, $M_{11}$ or $M_{23}$ with $n=11, 11$ or $23$
respectively. 2. $k=1$ and … 3. $k=2$ and …" (cases k = 1, 2 not used).

Lean form, k = 0 (a cycle fixing 0 points: an n-cycle), with weaker conclusions:
(a) n prime and, under a bijection with ZMod n, every element of G is affine x ↦ a x + b, a ≠ 0
(G ≤ AGL_1(p); the lower bound C_p ≤ G is dropped);
(b) a field F with q elements, d ≥ 2, n = (q^d - 1)/(q - 1), and a bijection with the points of
PG(d-1, F) = ℙ F (Fin d → F) under which every element of G is induced by a semilinear bijection
of F^d (G ≤ PΓL_d(q); the lower bound PGL_d(q) ≤ G is dropped);
(c) n = 11 or n = 23 (the identification of G is dropped). -/
def Jones2014_Thm1_2_k0 : Prop :=
  ∀ (n : ℕ) (H : Subgroup (Perm (Fin n))),
    MulAction.IsPreprimitive H (Fin n) →
    ¬ alternatingGroup (Fin n) ≤ H →
    (∃ c ∈ H, Perm.IsCycle c ∧ c.support = Finset.univ) →
    (n.Prime ∧ ∃ e : Fin n ≃ ZMod n, ∀ h ∈ H, ∃ a b : ZMod n, a ≠ 0 ∧ ∀ x, e (h x) = a * e x + b) ∨
    (∃ (q d : ℕ) (F : Type) (_ : Field F) (_ : Fintype F), Fintype.card F = q ∧ 2 ≤ d ∧
      n = (q ^ d - 1) / (q - 1) ∧
      ∃ e : Fin n ≃ ℙ F (Fin d → F), ∀ h ∈ H, IsCollineationVia e h) ∨
    (n = 11 ∨ n = 23)

end Riemann.Deep.LeanPrefix

