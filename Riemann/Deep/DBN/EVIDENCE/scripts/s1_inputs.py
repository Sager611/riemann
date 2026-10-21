"""Step 1: Planat's arXiv source is authentic and complete.
- sha256 of the arXiv e-print tarball (downloaded fresh from TARBALL_URL) equals TARBALL_SHA256;
- extract it to PLANAT's parent (if not yet extracted);
- every one of the 54 entries in Planat's own SHA256SUMS.json matches byte for byte. arXiv stores the eight
  *.json.gz members decompressed; they are re-compressed exactly as Planat's producer writes gzip (level 9,
  mtime 0, empty name) and must then match too.
Output: EVIDENCE/01_planat_inputs_sha256.txt"""
import gzip, io, json, tarfile
from common import *

guard(600)
log = Log(EVID / '01_planat_inputs_sha256.txt')
log('# Planat arXiv:2610.18214v1 source: authenticity of every file used (step s1_inputs.py)')
tb = read_bytes(TARBALL)
h = sha_bytes(tb)
log('tarball', TARBALL_URL, 'bytes', len(tb), 'sha256', h, 'expected', TARBALL_SHA256, 'MATCH' if h == TARBALL_SHA256 else 'MISMATCH')
assert h == TARBALL_SHA256
if not (PLANAT / 'SHA256SUMS.json').exists():
    with tarfile.open(TARBALL, 'r:gz') as t:
        t.extractall(PLANAT.parents[1], filter='data')
sums = planat_sums()
log('SHA256SUMS.json sha256', sha(PLANAT / 'SHA256SUMS.json'), 'entries', len(sums))
ok = 0
log('%-52s %9s %-64s %s' % ('file (Planat SHA256SUMS.json)', 'bytes', 'sha256', 'how checked'))
for rel, v in sums.items():
    p = PLANAT / rel
    if p.exists():
        b = read_bytes(p); how = 'file bytes'
    else:
        assert rel.endswith('.gz') and (PLANAT / rel[:-3]).exists(), rel
        buf = io.BytesIO()
        with gzip.GzipFile(filename='', fileobj=buf, mode='wb', mtime=0, compresslevel=9) as g:
            g.write(read_bytes(PLANAT / rel[:-3]))
        b = buf.getvalue(); how = 'gzip -9 of arXiv-decompressed member'
    good = sha_bytes(b) == v['sha256'] and len(b) == v['bytes']
    ok += good
    log('%-52s %9d %s %s %s' % (rel, len(b), sha_bytes(b), how, 'OK' if good else 'MISMATCH'))
log('RESULT', ok, 'of', len(sums), 'entries match' if ok == len(sums) else 'MISMATCH')
assert ok == len(sums)
log.save()
