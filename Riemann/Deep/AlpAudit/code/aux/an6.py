import re, collections
D = '/Users/adriansager/playground/riemann/.scratch/alpaudit/'
r = re.compile(r"(?<![\w.'])_nested\w*")
hits = [l for l in open(D + 'cls_bad.tsv') if '\t_nested\t' in l and r.search(l.split('\t', 2)[2])]
print('_nested tokens at identifier start:', len(hits), hits[:3])
dn = re.compile(r'\b(def|theorem|lemma|abbrev|instance|structure|class|inductive|opaque|axiom|alias)\s+(_root_\.)?([A-Za-z0-9_.]*\.)?(RiemannHypothesis|riemannZeta|LFunction|IsPrimitive|NeZero|re|DirichletCharacter|Complex|RiemannZetaZerosOnLine|DirichletLNonvanishingRightOfHalf|DirichletLZerosOnLine|PrimitiveDirichletLZerosOnLine|rhOfAll|allChi|critAll|claimAll|EveryCharacter|everyCharacter|rh_of_every|critLine)\b')
names = collections.Counter(); rows = []; nstr = 0; nchar = 0; strfiles = set(); exs = []
for l in open(D + 'scan3.txt', encoding='utf-8', errors='replace'):
    path, ln, txt = l.rstrip('\n').split(':', 2)
    code = txt.split('--')[0]
    m = dn.search(code)
    if m:
        full = (m.group(2) or '') + (m.group(3) or '') + m.group(4)
        names[(m.group(1), full)] += 1; rows.append((path[4:] + ':' + ln, txt.strip()[:150]))
    if '"' in code: nstr += 1; strfiles.add(path)
    if re.search(r"(^|[^A-Za-z0-9_'])'([^'\\]|\\.)'", code): nchar += 1; exs.append(path[4:] + ':' + ln + ' ' + txt.strip()[:100])
print('decl-name hits', len(rows))
for (k, n), v in names.most_common(): print(v, k, n)
for r_ in rows[:40]: print(r_)
print('string-literal lines', nstr, 'files', len(strfiles), '| char-literal lines', nchar, exs[:5])
