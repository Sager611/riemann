import os, re, sys, collections
ours_root = os.environ.get('OURS_OAI', os.path.join(os.path.dirname(__file__), '..', '..', 'QRH', 'OAI'))
trk_root = os.environ['TRK_OAI']  # <QuasiRiemannTracker clone>/submissions/2026-10-06-openai-baseline/src/Submissions/openai_baseline/OAI
def norm(s, kind):
    if kind == 'ours':
        s = s.replace('Riemann.Deep.QRH.OAI.', 'OAI.').replace('Riemann.Deep.QRH.Ext.', '')
    else:
        s = s.replace('Submissions.openai_baseline.OAI.', 'OAI.').replace('Submissions.openai_baseline.', '')
    s = re.sub(r'^\s*(public\s+)?import\s+.*$', '', s, flags=re.M)
    s = re.sub(r'^--.*$', '', s, flags=re.M)
    s = re.sub(r'[ \t]+\n', '\n', s)
    s = re.sub(r'\n{2,}', '\n', s)
    return s.strip()
def files(root):
    out = {}
    for dp, dn, fn in os.walk(root):
        for f in fn:
            if f.endswith('.lean'):
                p = os.path.join(dp, f); out[os.path.relpath(p, root)] = p
    return out
A = files(ours_root); B = files(trk_root)
print('ours', len(A), 'tracker', len(B), 'common', len(set(A) & set(B)))
print('only tracker:', sorted(set(B) - set(A)))
diff = []
for k in sorted(set(A) & set(B)):
    a = norm(open(A[k], encoding='utf-8').read(), 'ours'); b = norm(open(B[k], encoding='utf-8').read(), 'trk')
    if a != b: diff.append(k)
print('content differs after normalizing headers/imports:', len(diff))
for k in diff: print('  ', k)
