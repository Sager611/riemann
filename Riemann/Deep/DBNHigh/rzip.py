#!/usr/bin/env python3
# remote zip reader: list / extract members of a Zenodo zip via HTTP Range (curl -r); stdlib only.
# usage: rzip.py list PATTERN | rzip.py get OUTDIR MEMBER...
import subprocess, struct, sys, zlib, os, re, json
URL = 'https://zenodo.org/api/records/23752103/files/Proof_and_Replay.zip/content'
SIZE = 201060812
def rng(a, b):
    out = subprocess.run(['perl', '-e', 'alarm 120; exec @ARGV', 'curl', '-s', '-L', '-r', '%d-%d' % (a, b), URL], capture_output=True).stdout
    assert len(out) == b - a + 1, (len(out), b - a + 1)
    return out
def central():
    cache = os.environ.get('RZIP_CACHE', '/tmp/dbnhigh_cd.bin')
    if os.path.exists(cache):
        return open(cache, 'rb').read()
    tail = rng(SIZE - 65536, SIZE - 1)
    i = tail.rfind(b'PK\x05\x06')
    _, _, _, _, n, cdsize, cdoff, _ = struct.unpack('<IHHHHIIH', tail[i:i+22])
    if cdoff == 0xFFFFFFFF or n == 0xFFFF:
        j = tail.rfind(b'PK\x06\x06')
        n, cdsize, cdoff = struct.unpack('<QQQ', tail[j+32:j+56])
    cd = rng(cdoff, cdoff + cdsize - 1)
    open(cache, 'wb').write(cd)
    return cd
def entries():
    cd = central(); p = 0; out = []
    while p < len(cd) and cd[p:p+4] == b'PK\x01\x02':
        (method, crc, csize, usize, nlen, elen, clen, lhoff) = struct.unpack('<6xH4xIIIHHH8xI', cd[p+4:p+46])
        name = cd[p+46:p+46+nlen].decode('utf8', 'replace')
        extra = cd[p+46+nlen:p+46+nlen+elen]
        if 0xFFFFFFFF in (csize, usize, lhoff):
            q = 0
            while q < len(extra):
                hid, hl = struct.unpack('<HH', extra[q:q+4])
                if hid == 1:
                    vals = list(struct.unpack('<%dQ' % (hl // 8), extra[q+4:q+4+hl])); k = 0
                    if usize == 0xFFFFFFFF: usize = vals[k]; k += 1
                    if csize == 0xFFFFFFFF: csize = vals[k]; k += 1
                    if lhoff == 0xFFFFFFFF: lhoff = vals[k]; k += 1
                q += 4 + hl
        out.append(dict(name=name, method=method, csize=csize, usize=usize, off=lhoff))
        p += 46 + nlen + elen + clen
    return out
def get(e, outdir):
    h = rng(e['off'], e['off'] + 29)
    nlen, elen = struct.unpack('<HH', h[26:30])
    start = e['off'] + 30 + nlen + elen
    data = rng(start, start + e['csize'] - 1) if e['csize'] else b''
    if e['method'] == 8: data = zlib.decompress(data, -15)
    assert len(data) == e['usize']
    f = os.path.join(outdir, e['name'].replace('/', '__'))
    open(f, 'wb').write(data)
    return f
if __name__ == '__main__':
    es = entries()
    if sys.argv[1] == 'list':
        pat = re.compile(sys.argv[2]) if len(sys.argv) > 2 else None
        print(len(es), 'entries')
        for e in es:
            if pat is None or pat.search(e['name']):
                print(e['usize'], e['csize'], e['name'])
    else:
        outdir = sys.argv[2]; os.makedirs(outdir, exist_ok=True)
        want = set(sys.argv[3:])
        for e in es:
            if e['name'] in want:
                assert e['usize'] <= 25e6 and e['csize'] <= 25e6, e
                print('got', get(e, outdir), e['usize'])
