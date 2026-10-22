"""Extract every lstlisting block of the arXiv 2609.30306 TeX source into ex/bNN.txt.
Usage: python3 extract_appendix.py src/rhtm.tex  (b04-b06 source, b07 dispatch, b08 builder, b09-b13 tools, b20 table)."""
import sys
from pathlib import Path
src = sys.argv[1] if len(sys.argv) > 1 else 'src/rhtm.tex'
try:
    lines = Path(src).read_text().split('\n')
except OSError as e:
    sys.exit('cannot read %s: %s' % (src, e))
blocks = []
i = 0
while i < len(lines):
    if lines[i].strip().startswith('\\begin{lstlisting}'):
        j = i + 1
        while j < len(lines) and not lines[j].strip().startswith('\\end{lstlisting}'):
            j += 1
        blocks.append((i + 1, j + 1, lines[i + 1:j]))
        i = j + 1
    else:
        i += 1
out = Path('ex')
try:
    out.mkdir(exist_ok=True)
    for k, (a, b, body) in enumerate(blocks):
        (out / ('b%02d.txt' % k)).write_text('\n'.join(body) + '\n')
        print(k, a, b, len(body))
except OSError as e:
    sys.exit('cannot write ex/: %s' % e)
