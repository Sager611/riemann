
import sys
BT = chr(96)
PRE="Riemann.Deep.SunPrim.Lean."
def w(fn, txt): open("Riemann/Deep/SunPrim/Lean/"+fn+".lean","w").write(txt)
def mul8(a,b):
    r=0
    for i in range(3):
        if b>>i&1: r^=a<<i
    for d in (4,3):
        if r>>d&1: r^=(0b1011)<<(d-3)
    return r
def add9(x,y): return ((x%3+y%3)%3)+3*(((x//3)+(y//3))%3)
def mul9(x,y):
    a,b=x%3,x//3; c,d=y%3,y//3
    return ((a*c-b*d)%3)+3*((a*d+b*c)%3)
def powf(mul, x, e):
    r=1
    for _ in range(e): r=mul(r,x)
    return r
def gen(name, q, e, mul, add, leafk, nmods, desc):
    ST=[[add(powf(mul,x,e),y) for y in range(q)] for x in range(q)]
    enc=sum(ST[a][b] << (4*(q*a+b)) for a in range(q) for b in range(q))
    xs=["x%d"%i for i in range(q)]
    s=lambda i,j: "sT %s %s" % (xs[i], xs[j])
    E=lambda v: "(ofFin ⟨%d, by decide⟩)" % v
    S=["import "+PRE+name, "", "/-!", "# Search tree for circular orders of %s with the sums %sx ^ %d + y%s" % (desc, BT, e, BT), "",
       "%sK k x0 … xk%s: after choosing %sx0, …, xk%s (guards: %sxk%s is new, its sum %sx(k-1) ^ %d + xk%s is new), every way" % (BT,BT,BT,BT,BT,BT,BT,e,BT),
       "to continue fails; at the last level the closing sum %sx%d ^ %d + x0%s repeats. %sK1%s with %sx0 = 0%s says no circular" % (BT,q-1,e,BT,BT,BT,BT,BT),
       "order starting at 0 has pairwise distinct sums. Each level has its own %sDecidable%s instance." % (BT,BT), "-/", "",
       "namespace SunLean", "", "namespace %s" % name, "",
       "/-- All elements. -/", "def LL : List %s := [%s]" % (name, ", ".join(E(i) for i in range(q))), "",
       "theorem mem_LL : ∀ x : %s, x ∈ LL := by decide" % name, "",
       "/-- The table of %sx ^ %d + y%s, packed as in %s%sDefs%s. -/" % (BT, e, BT, BT, name, BT),
       "def STAB : ℕ := %d" % enc, "",
       "/-- %sx ^ %d + y%s by table lookup. -/" % (BT, e, BT),
       "def sT (a b : %s) : %s := ofFin ⟨(STAB >>> (4 * (%d * (a : Fin %d).val + (b : Fin %d).val))) %% 16 %% %d, Nat.mod_lt _ (by decide)⟩" % (name, name, q, q, q, q), "",
       "theorem sT_eq : ∀ a b : %s, sT a b = a ^ %d + b := by decide +kernel" % (name, e), ""]
    for k in range(q-1,0,-1):
        args=" ".join(xs[:k+1])
        g="%s ∉ [%s] → " % (xs[k], ", ".join(xs[j] for j in range(k-1,-1,-1)))
        if k>=2: g+="%s ∉ [%s] → " % (s(k-1,k), ", ".join(s(j,j+1) for j in range(k-2,-1,-1)))
        if k==q-1: g+="%s ∉ [%s] → False" % (s(q-1,0), ", ".join(s(j,j+1) for j in range(q-2,-1,-1)))
        else: g+="∀ %s ∈ LL, K%d %s %s" % (xs[k+1], k+1, args, xs[k+1])
        S.append("/-- Level %d. -/" % k)
        S.append("def K%d (%s : %s) : Prop :=" % (k,args,name))
        S.append("  " + g)
        S.append("")
        S.append("instance (%s : %s) : Decidable (K%d %s) := by unfold K%d; infer_instance" % (args,name,k,args,k))
        S.append("")
    S += ["end %s" % name, "", "end SunLean", ""]
    w(name+"Search", "\n".join(S))
    # validity of prefixes (x0 = 0): own guard at level k passes
    import itertools
    def own_ok(t):
        path=(0,)+t; k=len(t)
        if path[k] in path[:k]: return False
        if k>=2:
            sums=[ST[path[j]][path[j+1]] for j in range(k)]
            if sums[-1] in sums[:-1]: return False
        return True
    def valid(t): return all(own_ok(t[:m]) for m in range(1,len(t)+1))
    lname=lambda t: "k%d_%s" % (len(t), "_".join(map(str,t)))
    decided=[]; assembled=[]
    for k in range(1, leafk+1):
        for t in itertools.product(range(q), repeat=k):
            if not valid(t[:-1]): continue
            if not own_ok(t) or k==leafk: decided.append(t)
            else: assembled.append(t)
    per=(len(decided)+nmods-1)//nmods
    for mi in range(nmods):
        ch=decided[mi*per:(mi+1)*per]
        M=["import "+PRE+name+"Search", "", "/-! Search leaves, part %d: %sK k 0 x1 … xk%s by %sdecide +kernel%s (full subtrees at level %d; vacuous ones where the last guard fails). -/" % (mi,BT,BT,BT,BT,leafk), "",
           "namespace SunLean", "", "namespace %s" % name, ""]
        for t in ch:
            M.append("theorem %s : K%d %s %s := by decide +kernel" % (lname(t), len(t), E(0), " ".join(E(v) for v in t)))
        M += ["", "end %s" % name, "", "end SunLean", ""]
        w(name+"Leaf%d" % mi, "\n".join(M))
    T=["\n".join("import "+PRE+name+"Leaf%d" % mi for mi in range(nmods)), "", "/-! Assembly of the leaves: no circular order of %s starting at 0 has pairwise distinct sums. -/" % desc, "",
       "namespace SunLean", "", "namespace %s" % name, ""]
    for t in sorted(assembled, key=lambda t: -len(t)):
        k=len(t)
        us = " ".join(["_"]*(1 if k==1 else 2))
        T.append("theorem %s : K%d %s %s := fun %s y _ =>" % (lname(t), k, E(0), " ".join(E(v) for v in t), us))
        T.append("  match (y : Fin %d) with" % q)
        for v in range(q): T.append("  | ⟨%d, _⟩ => %s" % (v, lname(t+(v,))))
        T.append("  | ⟨n + %d, h⟩ => absurd h (by omega)" % q)
        T.append("")
    T.append("theorem no_order_zero : ∀ x1 ∈ LL, K1 %s x1 := fun y _ =>" % E(0))
    T.append("  match (y : Fin %d) with" % q)
    for v in range(q): T.append("  | ⟨%d, _⟩ => k1_%d" % (v, v))
    T.append("  | ⟨n + %d, h⟩ => absurd h (by omega)" % q)
    T += ["", "end %s" % name, "", "end SunLean", "", "#print axioms SunLean.%s.no_order_zero" % name, ""]
    w(name+"Top", "\n".join(T))
    print(name, "decided", len(decided), "assembled", len(assembled))
which=sys.argv[1]
if which=="F8": gen("F8", 8, 2, mul8, lambda a,b:a^b, 2, 4, "F₈")
if which=="F9": gen("F9", 9, 3, mul9, add9, 3, 18, "F₉")
print("ok", which)

