
BT = chr(96)
def block(name, q, e, Pdesc):
    xs=["x%d"%i for i in range(q)]
    L=[]
    L.append("/-- %s, %sP(x) = x^%d%s (%s, not of the form %sc - x%s): no circular order works. -/" % ("F₈" if q==8 else "F₉", BT, e, BT, Pdesc, BT, BT))
    L.append("theorem not_conj_A229827_i_%s : ¬ Conj_A229827_i := by" % name)
    L.append("  intro h")
    L.append("  have hcard : Fintype.card %s = %d := %s.card_eq" % (name, q, name))
    L.append("  have hdeg : (X ^ %d : %s[X]).degree < ((Fintype.card %s - 1 : ℕ) : WithBot ℕ) := by" % (e, name, name))
    L.append("    rw [degree_X_pow, hcard]; decide")
    L.append("  have key : ∀ c : %s, ∃ x : %s, x ^ %d ≠ c - x := by decide +kernel" % (name, name, e))
    L.append("  have hnot : ∀ c : %s, (X ^ %d : %s[X]) ≠ C c - X := by" % (name, e, name))
    L.append("    intro c hc")
    L.append("    obtain ⟨x, hx⟩ := key c")
    L.append("    apply hx")
    L.append("    simpa using congrArg (eval x) hc")
    L.append("  obtain ⟨a, ha, hall, hcov⟩ := h %s (X ^ %d) hdeg hnot" % (name, e))
    L.append("  simp only [eval_X_pow] at hcov")
    L.append("  have huniv : a.toFinset = Finset.univ :=")
    L.append("    Finset.eq_univ_iff_forall.2 (fun x => List.mem_toFinset.2 (hall x))")
    L.append("  have hlen : a.length = %d := by" % q)
    L.append("    rw [← List.toFinset_card_of_nodup ha, huniv, Finset.card_univ, hcard]")
    L.append("  obtain ⟨k, hk⟩ := exists_rotate_head (hall 0)")
    L.append("  have hb : (a.rotate k).Nodup := List.nodup_rotate.2 ha")
    L.append("  have hlenb : (a.rotate k).length = %d := by rw [List.length_rotate, hlen]" % q)
    L.append("  have hcovb : ∀ y : %s, y ∈ List.zipWith (fun x z => x ^ %d + z) (a.rotate k) ((a.rotate k).rotate 1) := by" % (name, e))
    L.append("    intro y; rw [sums_rotate]; exact List.mem_rotate.2 (hcov y)")
    L.append("  have hS := nodup_of_cover _ (by simp only [List.length_zipWith, List.length_rotate, hlen, min_self, hcard]) hcovb")
    L.append("  have hfun : (fun x z : %s => x ^ %d + z) = %s.sT := by" % (name, e, name))
    L.append("    funext x z; exact (%s.sT_eq x z).symm" % name)
    L.append("  rw [hfun] at hS")
    L.append("  obtain ⟨%s, hx⟩ := eq_of_length_%s _ hlenb" % (", ".join(xs), "eight" if q==8 else "nine"))
    L.append("  rw [hx] at hb hS hk")
    L.append("  simp only [List.head?_cons, Option.some.injEq] at hk")
    L.append("  subst hk")
    xz = ["0"] + xs[1:]
    sT = lambda i,j: "%s.sT %s %s" % (name, xz[i], xz[j])
    for kk in range(1, q):
        L.append("  have g%d : %s ∉ [%s] := (List.nodup_cons.1 (nodup_rev_take hb %d)).1" % (kk, xz[kk], ", ".join(xz[j] for j in range(kk-1,-1,-1)), kk+1))
    for kk in range(2, q):
        L.append("  have s%d : %s ∉ [%s] := (List.nodup_cons.1 (nodup_rev_take hS %d)).1" % (kk, sT(kk-1,kk), ", ".join(sT(j,j+1) for j in range(kk-2,-1,-1)), kk))
    L.append("  have sc : %s ∉ [%s] := (List.nodup_cons.1 (nodup_rev_take hS %d)).1" % (sT(q-1,0), ", ".join(sT(j,j+1) for j in range(q-2,-1,-1)), q))
    args=[]
    for kk in range(1, q):
        args.append("%s (%s.mem_LL %s) g%d" % (xs[kk], name, xs[kk], kk))
        if kk>=2: args.append("s%d" % kk)
    args.append("sc")
    L.append("  exact %s.no_order_zero %s" % (name, " ".join(args)))
    return "\n".join(L)

hdr = '''import Riemann.Deep.SunPrim.Lean.F8Top
import Riemann.Deep.SunPrim.Lean.F9Top
import Riemann.Deep.SunPrim.Lean.CircAux
import Mathlib.Algebra.Polynomial.Degree.Defs
import Mathlib.Algebra.Polynomial.Eval.Defs

/-!
# Sun, OEIS A229827, Conjecture (i), general part (finite fields)

Stated for every finite field ¤F¤ and every polynomial ¤P ∈ F[X]¤ (Mathlib ¤Polynomial¤) of degree
¤< q - 1¤ not of the form ¤c - x¤; refuted in ¤F₈¤ with ¤P = x²¤ and in ¤F₉¤ with ¤P = x³¤. ¤F₈¤, ¤F₉¤ are
explicit tables on ¤Fin 8¤, ¤Fin 9¤ whose field axioms are checked by ¤decide¤ (¤F8.lean¤, ¤F9.lean¤;
every field with ¤q¤ elements is isomorphic to these). ¤GaloisField 2 3¤ is not used: its import does
not fit the 290 MB cap and its elements are not computable. A circular order is rotated to start
at 0 (¤CircAux.lean¤); the search below 0 is ¤F8Top.lean¤, ¤F9Top.lean¤.
-/

namespace SunLean

open Polynomial

/-- Z.-W. Sun, OEIS A229827, comment "Conjecture: (i) …" (text format l. 9-11): for any finite field
¤F(q)¤ with ¤q¤ elements and a polynomial ¤P(x)¤ over ¤F(q)¤ of degree smaller than ¤q - 1¤, if ¤P(x)¤
is not of the form ¤c - x¤, then there is a circular permutation ¤a_1, …, a_q¤ of all the elements of
¤F(q)¤ with ¤{P(a_1) + a_2, …, P(a_(q-1)) + a_q, P(a_q) + a_1}¤ equal to ¤F(q)¤. -/
def Conj_A229827_i : Prop :=
  ∀ (F : Type) [Field F] [Fintype F] (P : F[X]),
    P.degree < ((Fintype.card F - 1 : ℕ) : WithBot ℕ) → (∀ c : F, P ≠ C c - X) →
    ∃ a : List F, a.Nodup ∧ (∀ x, x ∈ a) ∧
      ∀ y : F, y ∈ List.zipWith (fun x z => P.eval x + z) a (a.rotate 1)

'''.replace("¤", BT)
txt = hdr + block("F8", 8, 2, "Frobenius") + "\n\n" + block("F9", 9, 3, "Frobenius") + "\n\ntheorem not_conj_A229827_i : ¬ Conj_A229827_i := not_conj_A229827_i_F8\n\nend SunLean\n\n#print axioms SunLean.not_conj_A229827_i_F8\n#print axioms SunLean.not_conj_A229827_i_F9\n"
open("Riemann/Deep/SunPrim/Lean/A229827.lean","w").write(txt)
print("ok")

