import re, collections
D = '/Users/adriansager/playground/riemann/.scratch/alpaudit/'
S = '/Users/adriansager/playground/riemann/.scratch/alp/lean/src/'
ind = [l.rstrip('\n').split('\t', 1) for l in open(D + 'cls_ind.tsv')]
nested = []
for loc, txt in ind:
    fn, ln = loc.split(':')[:2]; ln = int(ln)
    lines = open(S + fn, encoding='utf-8').read().split('\n')
    m = re.match(r'^\s*(?:class\s+)?inductive\s+([^\s({\[:]+)', txt)
    if not m: continue
    name = m.group(1).split('.')[-1]
    block = [lines[ln - 1]]
    for x in lines[ln:ln + 60]:
        if re.match(r'^(theorem|lemma|def|instance|end|namespace|section|open|variable|noncomputable|abbrev|@\[|/--|structure|inductive|class|mutual)\b', x): break
        block.append(x)
    b = ' '.join(block)
    ctors = b.split('|')[1:]
    for c in ctors:
        head = c.split(':', 1)[-1]
        if re.search(r'\b(List|Array|Option|Multiset|Finset|Vector|Prod|Sum|Subtype|Sigma|PSigma|Set|Thunk)\b\s*\(?\s*(\(?\w+\s*:\s*[^)]*\)\s*→\s*)?' + re.escape(name) + r'\b', c) or re.search(r'(×|⊕)\s*' + re.escape(name) + r'\b|\b' + re.escape(name) + r'\b[^→]*\s(×|⊕)', c.split('→')[0] if '→' in c else ''):
            nested.append((loc, name, c.strip()[:150])); break
print('true-nested candidates', len(nested))
for n in nested: print(n)
print('== mutual block')
fn, ln = 'Development/ResolutionOfSingularities/EmbeddedResolution/IdealPrincipalization/PrincipalizationCertificate.lean', 10
L = open(S + fn, encoding='utf-8').read().split('\n')
for i in range(ln - 1, min(ln + 25, len(L))):
    if re.match(r'^\s*(mutual|inductive|structure|end|\|)', L[i]) or 'inductive' in L[i]: print(i + 1, L[i][:150])
