
BT = chr(96)
def tables8():
    def mul(a,b):
        r=0
        for i in range(3):
            if b>>i&1: r^=a<<i
        for d in (4,3):
            if r>>d&1: r^=(0b1011)<<(d-3)
        return r
    add=[[a^b for b in range(8)] for a in range(8)]
    mu=[[mul(a,b) for b in range(8)] for a in range(8)]
    neg=list(range(8))
    inv=[0]+[next(b for b in range(1,8) if mu[a][b]==1) for a in range(1,8)]
    return add,mu,neg,inv
def tables9():
    def add(x,y): return ((x%3+y%3)%3)+3*(((x//3)+(y//3))%3)
    def mul(x,y):
        a,b=x%3,x//3; c,d=y%3,y//3
        return ((a*c-b*d)%3)+3*((a*d+b*c)%3)
    ad=[[add(a,b) for b in range(9)] for a in range(9)]
    mu=[[mul(a,b) for b in range(9)] for a in range(9)]
    neg=[((-(x%3))%3)+3*((-(x//3))%3) for x in range(9)]
    inv=[0]+[next(b for b in range(1,9) if mu[a][b]==1) for a in range(1,9)]
    return ad,mu,neg,inv
PRE="Riemann.Deep.SunPrim.Lean."
def w(fn, txt): open("Riemann/Deep/SunPrim/Lean/"+fn+".lean","w").write(txt)
def gen(name, q, desc, tabs, e, Pdesc, nsplit):
    add,mu,neg,inv = tabs
    enc2=lambda m: sum(m[a][b] << (4*(q*a+b)) for a in range(q) for b in range(q))
    enc1=lambda v: sum(v[a] << (4*a) for a in range(q))
    D=[]
    D += ["import Mathlib.Algebra.Field.Defs","import Mathlib.Data.Fintype.Card","import Mathlib.Data.List.Nodup","",
          "/-! The field %s as packed tables on %sFin %d%s (definitions only). -/" % (desc, BT, q, BT), "",
          "namespace SunLean","",
          "/-- %s. -/" % desc,
          "def %s : Type := Fin %d" % (name,q), "", "namespace %s" % name, "",
          "instance : DecidableEq %s := inferInstanceAs (DecidableEq (Fin %d))" % (name,q),
          "instance : Fintype %s := inferInstanceAs (Fintype (Fin %d))" % (name,q), "",
          "/-- An element given by its code. -/", "def ofFin (i : Fin %d) : %s := i" % (q,name), "",
          "/-- Tables packed into one natural number, 4 bits per entry (entry %s(a, b)%s at bit %s4 (%d a + b)%s). -/" % (BT,BT,BT,q,BT),
          "def ADD : ℕ := %d" % enc2(add), "def MUL : ℕ := %d" % enc2(mu), "def NEG : ℕ := %d" % enc1(neg), "def INV : ℕ := %d" % enc1(inv), ""]
    for nm, T, two in [("addT","ADD",True),("mulT","MUL",True),("negT","NEG",False),("invT","INV",False)]:
        if two: D.append("def %s (a b : Fin %d) : Fin %d := ⟨(%s >>> (4 * (%d * a.val + b.val))) %% 16 %% %d, Nat.mod_lt _ (by decide)⟩" % (nm,q,q,T,q,q))
        else: D.append("def %s (a : Fin %d) : Fin %d := ⟨(%s >>> (4 * a.val)) %% 16 %% %d, Nat.mod_lt _ (by decide)⟩" % (nm,q,q,T,q))
    D += ["", "instance : Zero %s := ⟨(0 : Fin %d)⟩" % (name,q), "instance : One %s := ⟨(1 : Fin %d)⟩" % (name,q),
          "instance : Add %s := ⟨addT⟩" % name, "instance : Mul %s := ⟨mulT⟩" % name, "instance : Neg %s := ⟨negT⟩" % name,
          "instance : Inv %s := ⟨invT⟩" % name, "", "end %s" % name, "", "end SunLean", ""]
    w(name+"Defs", "\n".join(D))
    groups = {
      "Add": [("add_assoc'", "∀ a b c : %s, a + b + c = a + (b + c)"), ("zero_add'", "∀ a : %s, 0 + a = a"),
              ("add_zero'", "∀ a : %s, a + 0 = a"), ("neg_add_cancel'", "∀ a : %s, -a + a = 0"),
              ("add_comm'", "∀ a b : %s, a + b = b + a"), ("zero_ne_one'", "(0 : %s) ≠ 1")],
      "Mul": [("mul_assoc'", "∀ a b c : %s, a * b * c = a * (b * c)"), ("zero_mul'", "∀ a : %s, 0 * a = 0"),
              ("mul_zero'", "∀ a : %s, a * 0 = 0"), ("one_mul'", "∀ a : %s, 1 * a = a"), ("mul_one'", "∀ a : %s, a * 1 = a"),
              ("mul_comm'", "∀ a b : %s, a * b = b * a"), ("mul_inv_cancel'", "∀ a : %s, a ≠ 0 → a * a⁻¹ = 1"),
              ("inv_zero'", "(0 : %s)⁻¹ = 0")],
      "DistL": [("left_distrib'", "∀ a b c : %s, a * (b + c) = a * b + a * c")],
      "DistR": [("right_distrib'", "∀ a b c : %s, (a + b) * c = a * c + b * c")],
    }
    for g, ax in groups.items():
        T = ["import "+PRE+name+"Defs", "", "/-! Field axioms of %s%s%s, part %s (%sdecide%s). -/" % (BT,name,BT,g,BT,BT), "",
             "namespace SunLean", "", "namespace %s" % name, ""]
        for n_, st in ax: T.append("theorem %s : %s := by decide +kernel" % (n_, st % name))
        T += ["", "end %s" % name, "", "end SunLean", ""]
        w(name+g, "\n".join(T))
    F = ["import "+PRE+name+"Add", "import "+PRE+name+"Mul", "import "+PRE+name+"DistL", "import "+PRE+name+"DistR", "",
         "/-! %s is a field with %d elements. -/" % (desc, q), "", "namespace SunLean", "", "namespace %s" % name, "",
         "instance instField : Field %s where" % name]
    for f in ["add := (· + ·)", "add_assoc := add_assoc'", "zero := 0", "zero_add := zero_add'", "add_zero := add_zero'",
              "nsmul := nsmulRec", "neg := Neg.neg", "zsmul := zsmulRec", "neg_add_cancel := neg_add_cancel'",
              "add_comm := add_comm'", "mul := (· * ·)", "left_distrib := left_distrib'", "right_distrib := right_distrib'",
              "zero_mul := zero_mul'", "mul_zero := mul_zero'", "mul_assoc := mul_assoc'", "one := 1", "one_mul := one_mul'",
              "mul_one := mul_one'", "mul_comm := mul_comm'", "inv := Inv.inv", "exists_pair_ne := ⟨0, 1, zero_ne_one'⟩",
              "mul_inv_cancel := mul_inv_cancel'", "inv_zero := inv_zero'", "nnqsmul := _",
              "nnqsmul_def := fun _ _ => rfl", "qsmul := _", "qsmul_def := fun _ _ => rfl"]:
        F.append("  " + f)
    F += ["", "theorem card_eq : Fintype.card %s = %d := by decide" % (name, q), "", "end %s" % name, "", "end SunLean", "",
          "#print axioms SunLean.%s.instField" % name, ""]
    w(name, "\n".join(F))

gen("F8", 8, "F₈ = F₂[t]/(t³ + t + 1), element b₀ + b₁ t + b₂ t² ↦ b₀ + 2 b₁ + 4 b₂", tables8(), 2, "P(x) = x²", 2)
gen("F9", 9, "F₉ = F₃[i]/(i² + 1), element a + b i ↦ a + 3 b", tables9(), 3, "P(x) = x³", 3)
print("ok")

