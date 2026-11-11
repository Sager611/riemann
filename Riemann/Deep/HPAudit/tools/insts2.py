import os, re, sys, collections
exec(open(__import__('os').path.join(__import__('os').path.dirname(__file__), 'insts.py')).read().split("print('global instances:'")[0])
relevant = set('TopologicalSpace UniformSpace NormedAddCommGroup SeminormedAddCommGroup NormedSpace NormedField NontriviallyNormedField RCLike LT LE Preorder PartialOrder LinearOrder HPow Pow Div HDiv OfNat NatCast IntCast Field CommRing Ring Semiring AddCommGroup AddCommMonoid Module SMul Inv Neg Sub Add Mul One Zero Norm Dist MetricSpace PseudoMetricSpace EMetricSpace PseudoEMetricSpace Membership HasCompl Compl Bot Top Max Min Lattice CompleteLattice HasSum Summable T2Space ContinuousAdd IsTopologicalRing IsTopologicalAddGroup NumberField CharZero Algebra Coe CoeTail CoeHead CoeOut CoeSort CoeFun DecidableEq Decidable Nonempty Inhabited Subsingleton Unique OrderTopology MeasurableSpace'.split())
var_re = re.compile(r'^([A-Za-z_.]+)\s+(\S+)(\s+\S+)?$')
cnt = 0
for p, ln, concl, hdr in res:
    c = concl.strip()
    m = var_re.match(c)
    if not m: continue
    cls, a1 = m.group(1), m.group(2)
    if cls.split('.')[-1] not in relevant: continue
    args = [a1] + ([m.group(3).strip()] if m.group(3) else [])
    # a bound variable: appears in hdr binders like (a : Type ...) or {a : Type ...} or [..] — approximate: single-letter or Greek identifier
    if any(re.fullmatch(r'[α-ωA-Za-z]₀?|[α-ω][₀-₉]*|[A-Z][a-z]?', a) for a in args):
        cnt += 1
        print(p.split('/hpaudit/')[-1] + ':' + str(ln), '|', c, '||', hdr[:150])
print('generic-variable instances of relevant classes:', cnt)
