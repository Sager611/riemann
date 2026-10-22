# pyright: reportMissingImports=false
"""merge_parts.py OUTPREFIX PART1,PART2,...  -- sum grams_c part files into one single-parity file per parity
(OUTPREFIX_p0.bin, OUTPREFIX_p1.bin), keeping the provenance header (tag, N, prec, npieces, pieces = union).
Same checks as cert4.grams_parity (same tag/N/prec, pieces cover 0..npieces-1 exactly once).  Exact ball sums."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import ctx
import bc, cert4

def main():
    out, parts = sys.argv[1], sys.argv[2].split(',')
    heads = [cert4.read_head(f) for f in parts]
    tag = heads[0]['tag'][1]; S = bc.Setup(tag); prec = int(heads[0]['prec'][1]); ctx.prec = prec
    for p in (0, 1):
        Ar, Br, N, cov = cert4.grams_parity(parts, heads, S, p, prec)
        npc = int(heads[0]['npieces'][1])
        res = {f'A0raw_{p}': Ar, f'B0raw_{p}': Br, 'N': N, 'prec': prec, 'part': (0, 1), 'tag': tag, 'npieces': npc,
               'pieces': list(range(npc)), 'mode': '+'.join(sorted({h['mode'][1] for h in heads})), 'parities': [p]}
        fn = f"{out}_p{p}.bin"
        bc.save_mats(fn, res)
        print(f"parity {p}: {cov} -> {fn} ({os.path.getsize(fn)/2**30:.2f} GB)", flush=True)
        del Ar, Br, res

if __name__ == "__main__":
    main()
