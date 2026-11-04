#!/usr/bin/env python3
"""replay_kravatskiy.py -- optional second check (needs NumPy and a copy of Kravatskiy's 38D package).

Copies A. Kravatskiy's package verifications/improved/dim38-joint-latitudes (github.com/alexlegeartis/KissingNumbers,
commit b814435) to a temporary directory, replaces its equator file source/source/equator_lines_6855.txt by
equator_lines_13104.txt, updates the two certificate hashes and the claimed size, changes the hard-coded equator file
name and line count (6855 -> 13104) in refill/verify.py and refill/verify_coordinates.py, and runs both verifiers
unchanged otherwise.  Expected: 'VERIFIED: tau_38 >= 624020.' and 'VERIFIED independently: tau_38 >= 624020'.
Usage: python3 replay_kravatskiy.py /path/to/KissingNumbers/verifications/improved/dim38-joint-latitudes
"""
import sys, os, shutil, tempfile, subprocess, json, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
src = sys.argv[1]
tmp = tempfile.mkdtemp()
dst = os.path.join(tmp, 'pkg'); shutil.copytree(src, dst)
ss = os.path.join(dst, 'source', 'source')
os.remove(os.path.join(ss, 'equator_lines_6855.txt'))
shutil.copy(os.path.join(HERE, 'equator_lines_13104.txt'), os.path.join(ss, 'equator_lines_13104.txt'))
h = hashlib.sha256(open(os.path.join(ss, 'equator_lines_13104.txt'), 'rb').read()).hexdigest()
for f in ('refill/certificate.json', 'source/certificate.json'):
    p = os.path.join(dst, f); c = json.load(open(p))
    c['source_sha256'].pop('equator_lines_6855.txt'); c['source_sha256']['equator_lines_13104.txt'] = h
    if f.startswith('refill'): c['claimed_size'] = 611522 - 2 * 6855 + 2 * 13104
    json.dump(c, open(p, 'w'))
for f in ('refill/verify.py', 'refill/verify_coordinates.py'):
    p = os.path.join(dst, f); s = open(p).read()
    s = s.replace('equator_lines_6855.txt', 'equator_lines_13104.txt').replace('6855', '13104')
    open(p, 'w').write(s)
ok = True
for f in ('verify.py', 'verify_coordinates.py'):
    r = subprocess.run([sys.executable, f], cwd=os.path.join(dst, 'refill'), capture_output=True, text=True)
    print(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-500:])
    ok &= r.returncode == 0
shutil.rmtree(tmp)
sys.exit(0 if ok else 1)

