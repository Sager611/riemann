#!/usr/bin/env python3
# alpdep name-level dependency graph builder (static over-approximation).
# pass1: parse every file (namespaces, sections, opens, decls, aliases) -> names.tsv + hash arrays
# pass2: re-parse, tokenize each decl body (statement+proof), resolve tokens through namespace/open
#        prefixes (+ dot-field fallback over the decl's resolved type-like names) -> edge files.
import sys, os, re, time, array
import numpy as np
SRC, OUT, MODE = sys.argv[1], sys.argv[2], sys.argv[3]
LIM = int(sys.argv[4]) if len(sys.argv) > 4 else 0
FILES = '/Users/adriansager/playground/riemann/.scratch/alpcore/files.txt'
MASK = (1 << 64) - 1; MUL = 1000003
def hc(c): return hash(c) & MASK
def fold(h, comps):
    for c in comps: h = (h * MUL + hc(c) + 1) & MASK
    return h
BLOCK = re.compile(r"/-(?:(?!/-|-/).)*?-/", re.S)
LINEC = re.compile(r"--[^\n]*")
def strip(t):
    while '/-' in t:
        t2 = BLOCK.sub(lambda m: '\n' * m.group().count('\n'), t)
        if t2 == t: break
        t = t2
    return LINEC.sub('', t)
MODS = r"(?:(?:@\[[^\]]*\]\s*)|(?:private|protected|noncomputable|partial|unsafe|nonrec|scoped|local)\s+)*"
DECL = re.compile(r"^" + MODS + r"(theorem|lemma|def|abbrev|structure|class\s+inductive|class|inductive|instance|opaque|axiom|example|irreducible_def|alias)\b(.*)$")
CMD = re.compile(r"^(?:@\[|#|(?:private|protected|noncomputable|partial|unsafe|nonrec|scoped|local)\s|(?:theorem|lemma|def|abbrev|structure|class|inductive|instance|opaque|axiom|example|irreducible_def|namespace|section|end|open|export|variable|universe|mutual|attribute|notation|infix|infixl|infixr|prefix|postfix|macro|macro_rules|syntax|elab|set_option|initialize|alias|import|run_cmd|assert_not_exists|suppress_compilation|omit|include|unseal|seal|deriving)\b)")
KC = {'theorem':'T','lemma':'T','def':'D','abbrev':'A','structure':'S','class':'C','class inductive':'I','inductive':'I','instance':'N','opaque':'O','axiom':'X','example':'E','irreducible_def':'D','alias':'L'}
L = r"A-Za-z_\u03b1-\u03ba\u03bc-\u03c9\u0391-\u039f\u03a1\u03a4-\u03a9\u03ca-\u03fb\u1f00-\u1ffe\u2100-\u214f"
R = L + r"0-9'!?\u2080-\u209c\u1d62-\u1d6a\u2c7c"
ATOM = r"(?:«[^»\n]*»|[" + L + r"][" + R + r"]*)"
TOK = re.compile(r"(?<![" + R + r".])(\.?)(" + ATOM + r"(?:\.(?:" + ATOM + r"|[0-9]+))*)")
NAMEAT = re.compile(r"\s*(«[^»]*»|[^\s:(){}\[\]⦃⦄,]+)")
PRIO = re.compile(r"^\s*\(priority\s*:=[^)]*\)")
CTOR = re.compile(r"^\s+([^\s:(]+)\s*::")
FIELD = re.compile(r"^\s+([" + L + r"][" + R + r"]*(?:\s+[" + L + r"][" + R + r"]*)*)\s*[:({\[]")
INDC = re.compile(r"\|\s*([" + L + r"][" + R + r"]*)")
WHEREB = re.compile(r"^\s+([" + L + r"][" + R + r"]*)\s*(?::=|:)")
EXT = re.compile(r"\bextends\s+(.*?)(?:\bwhere\b|:=|$)")
KW = frozenset('fun by have show from let in if then else match with do at using only calc exact intro intros apply refine rw simp omega linarith nlinarith positivity constructor rfl obtain rcases cases induction unfold norm_num field_simp ring ring_nf push_cast exact_mod_cast gcongr by_contra by_cases use exists left right trivial decide aesop tauto specialize set change convert congr ext funext filter_upwards subst split split_ifs contradiction exfalso assumption simpa dsimp suffices this Type Prop Sort fun_prop bound refine_lift nth_rewrite rwa simp_all norm_cast lift'.split())
def unq(s): return s.replace('«', '').replace('»', '')
def parse(text, fidx):
    lines = strip(text).split('\n')
    scopes = [['root', [], [], []]]   # kind, comps, opens, vars
    decls = []; cur = None; var_tgt = None
    def close(i):
        nonlocal cur
        if cur is not None:
            cur[3] = i - 1; decls.append(cur); cur = None
    for i, s in enumerate(lines):
        if s and s[0] not in ' \t' and CMD.match(s):
            close(i); var_tgt = None
            m = DECL.match(s)
            if m:
                kw = re.sub(r'\s+', ' ', m.group(1)); rest = m.group(2)
                if kw == 'instance': rest = PRIO.sub('', rest)
                nm = None
                if kw not in ('example',):
                    mm = NAMEAT.match(rest)
                    if mm and mm.group(1)[0] not in ':([{⦃' and mm.group(1) not in ('where', ':='):
                        nm = unq(mm.group(1))
                ns = [c for sc in scopes if sc[0] == 'ns' for c in sc[1]]
                opens = [o for sc in scopes for o in sc[2]]
                vars_ = [v for sc in scopes for v in sc[3]]
                if nm is None: full = '_anon_%d_%d' % (fidx, i)
                elif nm.startswith('_root_.'): full = nm[7:]
                else: full = '.'.join(ns + [nm]) if ns else nm
                cur = [full, KC[kw], i, i, ns, opens, vars_, []]
                continue
            w = s.split(None, 1); k0 = w[0]; arg = w[1].strip() if len(w) > 1 else ''
            if k0 == 'namespace': scopes.append(['ns', unq(arg.split()[0]).split('.') if arg else [], [], []])
            elif k0 == 'section' or (k0 == 'noncomputable' and arg.startswith('section')): scopes.append(['sec', [], [], []])
            elif k0 == 'mutual': scopes.append(['mut', [], [], []])
            elif k0 == 'end':
                if len(scopes) > 1: scopes.pop()
            elif k0 in ('open', 'export'):
                for t in TOK.finditer(arg):
                    o = unq(t.group(2))
                    if o not in ('scoped', 'hiding', 'renaming', 'in') and not t.group(1): scopes[-1][2].append(o)
            elif k0 == 'variable':
                scopes[-1][3].append(arg); var_tgt = scopes[-1][3]
        else:
            if var_tgt is not None and s.strip(): var_tgt.append(s)
    close(len(lines))
    out = []
    for d in decls:
        full, kind, st, en = d[0], d[1], d[2], d[3]
        body = lines[st:en + 1]
        al = []
        if kind in ('S', 'C'):
            ctor = None
            for ln in body[1:]:
                m = CTOR.match(ln)
                if m: ctor = unq(m.group(1)); continue
                m = FIELD.match(ln)
                if m: al.extend(m.group(1).split())
            al.append(ctor or 'mk')
            me = EXT.search(' '.join(body[:3]))
            if me:
                for p in me.group(1).split(','):
                    t = TOK.search(p)
                    if t: al.append('to' + unq(t.group(2)).split('.')[-1])
            al += ['rec', 'casesOn', 'recOn', 'noConfusion', 'ext', 'ext_iff']
        elif kind == 'I':
            for ln in body:
                al.extend(m.group(1) for m in INDC.finditer(ln))
            al += ['rec', 'casesOn', 'recOn', 'noConfusion']
        else:
            wh = False
            for ln in body:
                if wh:
                    m = WHEREB.match(ln)
                    if m: al.append(m.group(1))
                if re.search(r'\bwhere\s*$', ln): wh = True
        d[7] = al
        out.append((full, kind, st, en, d[4], d[5], d[6], al, body))
    return out
def files_iter():
    with open(FILES) as fh:
        for idx, ln in enumerate(fh):
            if LIM and idx >= LIM: break
            yield idx, ln.rstrip('\n')
t0 = time.time()
if MODE == 'pass1':
    H = array.array('Q'); I = array.array('i'); HID = array.array('Q'); K = bytearray(); LAST = array.array('Q'); NS = set()
    nper = array.array('i'); nid = 0
    nf = open(OUT + '/names.tsv', 'w')
    for fidx, p in files_iter():
        try: txt = open(os.path.join(SRC, p), encoding='utf-8', errors='replace').read()
        except Exception: nper.append(0); continue
        ds = parse(txt, fidx); nper.append(len(ds))
        for (full, kind, st, en, ns, opens, vars_, al, body) in ds:
            comps = full.split('.'); h = 0
            for c in comps[:-1]:
                h = (h * MUL + hc(c) + 1) & MASK; NS.add(h)
            h = fold(h, comps[-1:])
            H.append(h); I.append(nid); HID.append(h); K.append(ord(kind)); LAST.append(hc(comps[-1]))
            for a in al:
                H.append((h * MUL + hc(a) + 1) & MASK); I.append(nid); LAST.append(hc(a))
            nf.write('%d\t%s\t%s\t%d\t%d\t%d\n' % (nid, full, kind, fidx, st + 1, en + 1))
            nid += 1
        if fidx % 20000 == 0: print('pass1', fidx, nid, round(time.time() - t0), file=sys.stderr, flush=True)
    nf.close()
    Hn = np.frombuffer(H, dtype=np.uint64); In = np.frombuffer(I, dtype=np.int32)
    o = np.argsort(Hn, kind='stable')
    np.save(OUT + '/h_sorted.npy', Hn[o]); np.save(OUT + '/id_sorted.npy', In[o])
    np.save(OUT + '/hid.npy', np.frombuffer(HID, dtype=np.uint64)); np.save(OUT + '/kind.npy', np.frombuffer(bytes(K), dtype=np.uint8))
    np.save(OUT + '/last.npy', np.unique(np.frombuffer(LAST, dtype=np.uint64))); np.save(OUT + '/nper.npy', np.frombuffer(nper, dtype=np.int32))
    np.save(OUT + '/ns.npy', np.array(sorted(NS), dtype=np.uint64))
    print('pass1 done decls', nid, 'keys', len(H), 'secs', round(time.time() - t0), file=sys.stderr)
else:
    hs = np.load(OUT + '/h_sorted.npy'); ids = np.load(OUT + '/id_sorted.npy'); hid = np.load(OUT + '/hid.npy')
    kind = np.load(OUT + '/kind.npy'); nper = np.load(OUT + '/nper.npy'); last = np.load(OUT + '/last.npy')
    NSS = set(int(x) for x in np.load(OUT + '/ns.npy'))
    TYPEK = np.zeros(256, dtype=bool)
    for c in 'SCIDAO': TYPEK[ord(c)] = True
    nh = len(hs); nl = len(last)
    def look(cands):  # -> list of ids whose full-name hash is in cands
        if not cands: return []
        a = np.fromiter(cands, dtype=np.uint64, count=len(cands))
        lo = np.searchsorted(hs, a, 'left'); hi = np.searchsorted(hs, a, 'right')
        m = hi > lo
        if not m.any(): return []
        out = []
        for l_, h_ in zip(lo[m].tolist(), hi[m].tolist()): out.extend(ids[l_:h_].tolist())
        return out
    def inlast(strs):
        if not strs: return set()
        L_ = list(strs); a = np.fromiter((hc(s) for s in L_), dtype=np.uint64, count=len(L_))
        pos = np.searchsorted(last, a); pos[pos >= nl] = nl - 1
        ok = last[pos] == a
        return set(s for s, o in zip(L_, ok.tolist()) if o)
    phcache = {}
    print('loaded', round(time.time() - t0), file=sys.stderr, flush=True)
    cnt = open(OUT + '/ecount.bin', 'wb'); tgt = open(OUT + '/etarget.bin', 'wb')
    nid = 0; nfb = 0; nedge = 0
    for fidx, p in files_iter():
        try: txt = open(os.path.join(SRC, p), encoding='utf-8', errors='replace').read()
        except Exception: txt = ''
        ds = parse(txt, fidx) if txt else []
        assert len(ds) == nper[fidx], (p, len(ds), nper[fidx])
        ecnt = array.array('i'); etg = array.array('i')
        for (full, kind_, st, en, ns, opens, vars_, al, body) in ds:
            fc = full.split('.')
            ck = (tuple(fc[:-1]), tuple(ns), tuple(opens))
            ph = phcache.get(ck)
            if ph is None:
                ph = [fold(0, ns[:j]) for j in range(len(ns), -1, -1)]
                if len(fc) > len(ns) + 1:
                    for j in range(len(fc) - 1, len(ns), -1): ph.append(fold(0, fc[:j]))
                base = list(ph)
                for o in opens:
                    oc = o.split('.')
                    for b in base:
                        x = fold(b, oc)
                        if x in NSS: ph.append(x)
                ph = list(dict.fromkeys(ph))
                if len(phcache) > 2000: phcache.clear()
                phcache[ck] = ph
            text = '\n'.join(body)
            if vars_: text += '\n' + '\n'.join(vars_)
            toks = set()
            for m in TOK.finditer(text):
                toks.add((m.group(1), m.group(2)))
            tl = []; comps_all = set()
            for lead, nm in toks:
                comps = unq(nm).split('.')
                if not lead and len(comps) == 1 and comps[0] in KW: continue
                tl.append((lead, comps)); comps_all.update(c for c in comps if not c.isdigit())
            okc = inlast(comps_all)
            cands = set(); pend = set()
            for lead, comps in tl:
                if lead:
                    for c in comps:
                        if c in okc: pend.add(c)
                    continue
                if not any(c in okc for c in comps): continue
                hcs = [hc(c) + 1 for c in comps]; oks = [c in okc for c in comps]
                for b in ph:
                    h = b
                    for k in range(len(comps)):
                        h = (h * MUL + hcs[k]) & MASK
                        if oks[k]: cands.add(h)
                for c in comps[1:]:
                    if c in okc: pend.add(c)
            res = set(look(cands))
            if pend and res:
                th = [int(hid[x]) for x in res if TYPEK[kind[x]]]
                if th:
                    if len(th) * len(pend) <= 8000:
                        pf = [hc(f) + 1 for f in pend]
                        res.update(look([(hx * MUL + hf) & MASK for hx in th for hf in pf]))
                    else: nfb += 1
            res.discard(nid)
            r = sorted(res); etg.extend(r); ecnt.append(len(r))
            nedge += len(r); nid += 1
        tgt.write(etg.tobytes()); cnt.write(ecnt.tobytes())
        if fidx % 10000 == 0: print('pass2', fidx, nid, nedge, nfb, round(time.time() - t0), file=sys.stderr, flush=True)
    cnt.close(); tgt.close()
    print('pass2 done decls', nid, 'edges', nedge, 'fallback-skipped', nfb, 'secs', round(time.time() - t0), file=sys.stderr)
