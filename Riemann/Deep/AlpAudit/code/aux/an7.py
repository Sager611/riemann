import re, collections
D = '/Users/adriansager/playground/riemann/.scratch/alpaudit/'
q = collections.Counter()
for l in open(D + 'cls_bad.tsv'):
    p = l.rstrip('\n').split('\t')
    if len(p) >= 3 and p[1] == 'Quot':
        for m in re.finditer(r'\bQuot\.(\w+)', p[2]): q['Quot.' + m.group(1)] += 1
print('Quot breakdown', q.most_common(12))
nums = [l.rstrip('\n').split('\t') for l in open(D + 'cls_num.tsv')]
nf = collections.Counter(n[0].split(':')[0] for n in nums)
mx = max(nums, key=lambda n: len(n[1]))
print('bignum occurrences', len(nums), 'files', len(nf), 'longest', mx[0], len(mx[1]), 'digits', '| top files', nf.most_common(3))
print('bignum examples', [(n[0], n[2][:90]) for n in nums[:3]])
dec = [l.rstrip('\n').split('\t') for l in open(D + 'cls_dec.tsv')]
kinds = collections.Counter()
for d in dec:
    t = d[1]
    if re.search(r'by\s+decide|:=\s*by\s+decide|<;>\s*decide|\bdecide\b\s*$', t): kinds['tactic decide'] += 1
    elif 'native_decide' in t: kinds['native'] += 1
    elif re.search(r'decide_eq_true|of_decide_eq_true', t): kinds['decide_eq_true lemma'] += 1
    elif re.search(r'simp|omega|norm_num', t): kinds['decide inside simp/other'] += 1
    else: kinds['other'] += 1
print('decide kinds', kinds.most_common())
print('decide examples', [(d[0], d[1][:80]) for d in dec[:2]])
wf = [l.rstrip('\n').split('\t') for l in open(D + 'cls_wf.tsv')]
print('wf files', sorted(set(w[0].split(':')[0] for w in wf))[:8])
