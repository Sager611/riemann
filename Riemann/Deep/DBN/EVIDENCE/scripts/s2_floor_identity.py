"""Step 2: the two dbn4 evidence source files are Planat's refined sources with kernels appended, nothing else.
- rebuild them from Planat's refined inputs + dbn4/catalogue.json exactly as dbn4/run_phase1.py does;
- their sha256 must equal the hashes recorded in both dbn4 history headers (what the auditors bind to);
- byte level: evidence = original[:a] + insertion + original[a:], the insertion sits right before the closing
  bracket of the "kernels" array and parses as exactly the 786 appended kernels;
- every top-level key except "kernels" is identical; the first 57 kernels are identical;
- every floor (608 cells x 10 shifts x 2 providers) and height is identical;
- every floor and height CONSUMED by the two dbn4 histories equals Planat's original value.
Output: EVIDENCE/02_floor_identity.txt"""
from fractions import Fraction as F
from common import *

guard(600)
log = Log(EVID / '02_floor_identity.txt')
log('# Floor identity: dbn4 evidence sources vs Planat refined sources (step s2_floor_identity.py)')
cat = load_json(DBN4 / 'catalogue.json')
log('catalogue.json sha256', sha(DBN4 / 'catalogue.json'), 'kernels', cat['count'])
new = [{'r': k['r'], 'rho': k['rho'], 'lambda': k['lam'], 'self_cap': k['self_cap'], 'dbn4_certificate': True,
        'leaves': k['leaves']} for k in cat['kernels']]
ev = WORK / 'evidence'; ev.mkdir(parents=True, exist_ok=True)
H1 = load_json(HIST_SINGLE); H2 = load_json(HIST_MIX)
O, E = [], []
for pi, n in enumerate(NAMES):
    ob = read_bytes(PLANAT / 'inputs/refined' / n)
    o = loads(ob)
    log('\n## provider', pi, n)
    log('Planat refined file bytes', len(ob), 'sha256', sha_bytes(ob), '(listed in SHA256SUMS.json:',
        sha_bytes(ob) == planat_sums()['inputs/refined/' + n]['sha256'], ')')
    log('Planat file is already canonical json (sorted keys, compact):', canon(o) == ob)
    e = loads(ob); e['kernels'] = e['kernels'] + new
    eb = canon(e)
    (ev / n).write_bytes(eb)
    he = sha_bytes(eb)
    log('rebuilt evidence bytes', len(eb), 'sha256', he)
    log('  == hash in classical_dbn4 header:', he == H1['source_hashes'][pi], '; == hash in classical_dbn4_mix header:',
        he == H2['source_hashes'][pi], '; == recorded EVIDENCE_SHA256:', he == EVIDENCE_SHA256[pi])
    assert he == H1['source_hashes'][pi] == H2['source_hashes'][pi] == EVIDENCE_SHA256[pi]
    old = REPO / '.scratch/dbn4/evidence' / n
    if old.exists():
        log('  == the file used in the original dbn4 run (.scratch/dbn4/evidence):', sha(old) == he)
    a = 0
    while a < len(ob) and ob[a] == eb[a]:
        a += 1
    s = 0
    while s < len(ob) - a and ob[-1 - s] == eb[-1 - s]:
        s += 1
    ins = eb[a:len(eb) - s]
    log('byte diff: common prefix', a, 'bytes, common suffix', s, 'bytes, prefix+suffix == whole original:', a + s == len(ob))
    log('  => evidence = original[:%d] + insertion(%d bytes) + original[%d:]; no original byte removed or changed' % (a, len(ins), a))
    start = ob.index(b'"kernels":[')
    depth, j = 0, start + len(b'"kernels":')
    while True:
        c = ob[j:j + 1]
        if c == b'[': depth += 1
        elif c == b']':
            depth -= 1
            if depth == 0: break
        j += 1
    log('  "kernels" array of the original closes at byte', j, '; insertion offset', a, '; equal:', j == a)
    log('  context around the insertion point:', repr(ob[a - 30:a + 15].decode()))
    parsed = loads('[' + ins[1:].decode() + ']')
    log('  insertion starts with ",":', ins[:1] == b',', '; parses as', len(parsed), 'kernels; equals the 786 appended:', parsed == new)
    assert a + s == len(ob) and j == a and parsed == new and ins[:1] == b','
    for key in sorted(o):
        if key == 'kernels':
            continue
        log('  key %-42s sha256(canonical) %s identical: %s' % (key, sha_bytes(canon(o[key]))[:16], canon(o[key]) == canon(e[key])))
        assert canon(o[key]) == canon(e[key])
    log('  kernels: original', len(o['kernels']), 'evidence', len(e['kernels']), '; first 57 identical:', e['kernels'][:57] == o['kernels'])
    modes = sorted(k for k in o['rows'][0] if k.startswith('r_'))
    nf = sum(1 for i in range(608) for m in modes if e['rows'][i][m] == o['rows'][i][m])
    nh = sum(1 for i in range(608) if e['plan']['rows'][i]['height'] == o['plan']['rows'][i]['height'])
    log('  floors identical:', nf, 'of', 608 * len(modes), '(608 cells x shifts', ','.join(m[2:].replace('_', '/') for m in modes) + '); heights identical:', nh, 'of 608')
    assert nf == 608 * len(modes) and nh == 608
    O.append(o); E.append(e)
log('\n## floors and heights consumed by the histories, compared with Planat\'s ORIGINAL refined files')
for label, H in (('classical_dbn4 (0.148409)', H1), ('classical_dbn4_mix (0.148369)', H2)):
    nf = nh = nc = bad = 0
    used = set()
    for rec in H['records']:
        i, pi = rec['cell'], rec['provider']
        if rec.get('type', 'single') == 'single':
            k = E[pi]['kernels'][rec['kernel']]
            mode = 'r_' + k['r'].replace('/', '_')
            bad += rec['floor'] != O[pi]['rows'][i][mode]; nf += 1
            bad += F(rec['self_cap']) != F(k['self_cap']); nc += 1
        else:
            m = H['mixed_kernels'][rec['mix']]
            if rec['mix'] not in used:
                used.add(rec['mix'])
                for t in ('1', '2'):
                    mode = 'r_' + m['r' + t].replace('/', '_')
                    bad += m['floor' + t] != O[pi]['rows'][m['cell']][mode]; nf += 1
                bad += F(m['height']) != F(O[pi]['plan']['rows'][m['cell']]['height']); nh += 1
        bad += F(rec['height']) != F(O[pi]['plan']['rows'][i]['height']); nh += 1
    log('%-32s records %5d: floors compared %5d, heights %5d, own-pair caps %5d, mismatches %d' % (label, len(H['records']), nf, nh, nc, bad))
    assert bad == 0
log('RESULT floor identity holds byte for byte; the only change is the appended kernel list')
log.save()
