#!/usr/bin/env python3
r"""Apply the post-RH draft (Riemann/Deep/PostRH/DRAFT.md) once ROOT (TREE.md, Alpöge's Lean RH + Dirichlet GRH) is verified.

  python3 Riemann/Deep/PostRH/apply.py                        dry run (default): checks everything, prints what it would change, writes nothing
  python3 Riemann/Deep/PostRH/apply.py --write [--date '21 Nov 2026']   applies
  python3 Riemann/Deep/PostRH/apply.py --preview DIR          dry run that also writes DIR/TREE.md and DIR/COROLLARIES.md (the would-be files) for diffing

What --write does (exact-substring edits only):
  TREE.md        §1 of DRAFT.md: "<sentence> " inserted right after each node anchor; on the ROOT line right after the
                 anchor's leading [ … ] tag, which is never edited (the orchestrator edits it by hand).
  COROLLARIES.md §2 of DRAFT.md appended under '## If ROOT holds (Alpöge, verified <date>)' (default date: today, UTC).
Markdown pipe escapes '\|' in the draft become '|' (TREE.md and COROLLARIES.md are lists, not tables).
Every assertion (anchor count == 1, sentence absent, no 'If ROOT holds' already on the node line, block absent, ROOT tag
intact, files unchanged since read) runs before
the first write; any failure exits with status 1 and writes nothing.
"""
import argparse
import datetime
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DRAFT = os.path.join(HERE, "DRAFT.md")
TREE = os.path.join(REPO, "TREE.md")
COR = os.path.join(REPO, "COROLLARIES.md")
ROOT_ANCHOR = "- ROOT (RH itself) — "
HEADING_PREFIX = "## If ROOT holds (Alpöge, verified "
ROW = re.compile(r"^\| \x60(.+?)\x60 \| (.+) \|$")


def fail(msg: str) -> None:
    print("FAIL: " + msg)
    print("Nothing written.")
    sys.exit(1)


def read(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        return f.read()


def sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def unescape(text: str) -> str:
    return text.replace("\\|", "|")


def section(text: str, start: str, end: str) -> str:
    if text.count(start) != 1 or text.count(end) != 1:
        fail(f"DRAFT.md: section markers {start!r}/{end!r} not unique")
    return text.split(start, 1)[1].split(end, 1)[0]


def parse_draft(text: str) -> tuple[list[tuple[str, str]], list[str]]:
    rows: list[tuple[str, str]] = []
    for line in section(text, "\n## 1.", "\n## 2.").split("\n"):
        m = ROW.match(line)
        if m:
            rows.append((m.group(1), unescape(m.group(2)).strip()))
    block = [unescape(l).rstrip() for l in section(text, "\n## 2.", "\n## 3.").split("\n") if l.startswith("- If ROOT holds")]
    return rows, block


def root_tag(tree: str) -> str:
    """The bracketed tag right after the ROOT anchor ('' if none), bracket-depth aware, on the same line."""
    i = tree.index(ROOT_ANCHOR) + len(ROOT_ANCHOR)
    if not tree.startswith("[", i):
        return ""
    depth = 0
    for j in range(i, len(tree)):
        c = tree[j]
        if c == "\n":
            break
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                return tree[i : j + 1]
    fail("ROOT tag: no closing ']' on the ROOT line")
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true", help="apply the edits (default: dry run)")
    ap.add_argument("--preview", metavar="DIR", help="dry run; also write the would-be TREE.md and COROLLARIES.md into DIR (never the repo root)")
    ap.add_argument("--date", default=datetime.datetime.now(datetime.timezone.utc).strftime("%-d %b %Y"), help="verification date for the COROLLARIES heading")
    args = ap.parse_args()

    draft, tree, cor = read(DRAFT), read(TREE), read(COR)
    tree_sha, cor_sha = sha(tree), sha(cor)
    rows, block = parse_draft(draft)
    mode = "WRITE" if args.write else "DRY RUN"
    print(f"apply.py {mode}  {datetime.datetime.now(datetime.timezone.utc):%Y-%m-%d %H:%M} UTC")
    print(f"DRAFT.md sha {sha(draft)}; TREE.md sha {tree_sha} ({tree.count(chr(10))} lines); COROLLARIES.md sha {cor_sha} ({cor.count(chr(10))} lines)")

    # ---- assertions, all before any write ----
    if not rows:
        fail("DRAFT.md §1: no rows parsed")
    if not block:
        fail("DRAFT.md §2: no '- If ROOT holds' lines parsed")
    anchors = [a for a, _ in rows]
    if len(set(anchors)) != len(anchors):
        fail("DRAFT.md §1: duplicate anchors")
    for a, s in rows:
        if "\n" in s or not s:
            fail(f"bad sentence for {a!r}")
        if tree.count(a) != 1:
            fail(f"TREE.md: anchor {a!r} occurs {tree.count(a)} times (need 1)")
        if s in tree:
            fail(f"TREE.md: sentence for {a!r} already present")
        node_line = tree[tree.index(a) :].split("\n", 1)[0]
        if "If ROOT holds" in node_line or s.split(":", 1)[0] + ": RH and Dirichlet GRH" in node_line:
            fail(f"TREE.md: line of {a!r} already carries a post-ROOT sentence (edited by hand?)")
    tag = root_tag(tree) if ROOT_ANCHOR in anchors else ""
    if ROOT_ANCHOR in anchors and tag:
        if not tree.startswith(" ", tree.index(ROOT_ANCHOR) + len(ROOT_ANCHOR) + len(tag)):
            fail("ROOT tag not followed by a space")
    if HEADING_PREFIX in cor:
        fail("COROLLARIES.md: an 'If ROOT holds (Alpöge, verified …)' heading is already present")
    for line in block:
        if line in cor:
            fail(f"COROLLARIES.md: block line already present: {line[:80]!r}")
    if len(set(block)) != len(block):
        fail("DRAFT.md §2: duplicate lines")

    # ---- compute the new texts in memory (exact-substring replacements, each old string unique) ----
    new_tree = tree
    plan = []
    for a, s in rows:
        old = a + tag + " " if (a == ROOT_ANCHOR and tag) else a
        if new_tree.count(old) != 1:
            fail(f"TREE.md: edit point {old[:60]!r}… occurs {new_tree.count(old)} times")
        before = len(new_tree)
        new_tree = new_tree.replace(old, old + s + " ", 1)
        if len(new_tree) != before + len(s) + 1:
            fail(f"TREE.md: unexpected length change at {a!r}")
        line_no = new_tree[: new_tree.index(a)].count("\n") + 1
        plan.append((line_no, a, old, s))
    if new_tree.count("\n") != tree.count("\n"):
        fail("TREE.md: line count would change")
    if tag and not new_tree.split(ROOT_ANCHOR, 1)[1].startswith(tag + " "):
        fail("ROOT tag would move or change")
    heading = HEADING_PREFIX + args.date + ")"
    sep = "" if cor.endswith("\n") else "\n"
    new_cor = cor + sep + "\n" + heading + "\n" + "\n".join(block) + "\n"

    # ---- report ----
    print(f"TREE.md: {len(plan)} insertions, +{len(new_tree) - len(tree)} chars, line count unchanged")
    for line_no, a, old, s in plan:
        where = f"after ROOT tag {tag[:40]!r}…" if old != a else "after anchor"
        print(f"  l.{line_no:<3} {a[2:].split(' ')[0]:<5} count 1, {where}: + {s[:80]!r}")
    print(f"COROLLARIES.md: append heading {heading!r} and {len(block)} lines (+{len(new_cor) - len(cor)} chars)")
    for line in block:
        print(f"  + {line[:80]!r}")
    print("ALL ASSERTIONS PASSED")

    if args.preview:
        out = os.path.abspath(args.preview)
        if out in (REPO, HERE) or args.write:
            fail("--preview DIR must be a scratch directory and excludes --write")
        os.makedirs(out, exist_ok=True)
        for name, text in (("TREE.md", new_tree), ("COROLLARIES.md", new_cor)):
            with open(os.path.join(out, name), "w", encoding="utf-8") as f:
                f.write(text)
        print(f"PREVIEW written to {out} (repo files untouched)")
    if not args.write:
        print("DRY RUN: nothing written. Apply with: python3 Riemann/Deep/PostRH/apply.py --write [--date '<d Mon YYYY>']")
        return
    if sha(read(TREE)) != tree_sha or sha(read(COR)) != cor_sha:
        fail("TREE.md or COROLLARIES.md changed while running; rerun")
    for path, text in ((TREE, new_tree), (COR, new_cor)):
        tmp = path + ".apply.tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            f.write(text)
        os.replace(tmp, path)
    print("WRITTEN: TREE.md and COROLLARIES.md (the ROOT tag is untouched; edit it by hand)")


if __name__ == "__main__":
    main()
