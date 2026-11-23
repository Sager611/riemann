"""Check PROOF.md against the Lean sources.

Checks: relative links resolve, Lean links land inside the cited declaration,
every goal line has exactly 7 words. Prints the top-level word count.
Run: python3 papers/02_gale-fails-every-dimension/check.py
"""

from __future__ import annotations

import pathlib
import re
import sys

BT = chr(96)
HERE = pathlib.Path(__file__).resolve().parent
DOC = HERE / "PROOF.md"
KEYWORDS = "theorem|lemma|def|abbrev|structure|inductive|instance"
MODIFIERS = "(noncomputable |private |protected )*"
TOP = re.compile(
    r"^(/--|/-!|noncomputable |private |protected |theorem |lemma |def |abbrev "
    r"|structure |inductive |instance |end( |$)|section|namespace |variable "
    r"|open |include |@\[)",
)
LINK = re.compile(r"\[([^\]]+)\]\(([^)#\s]+)(?:#L(\d+)(?:-L(\d+))?)?\)")
GOAL_WORDS = 7


def to_int(text: str | None) -> int | None:
    """Parse a line number, or return None."""
    if text is None:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def doc_start(lines: list[str], decl: int) -> int:
    """Return the first line of the docstring right above line index decl."""
    if decl == 0 or not lines[decl - 1].rstrip().endswith("-/"):
        return decl
    k = decl - 1
    while k >= 0 and not lines[k].lstrip().startswith("/--"):
        k -= 1
    return k if k >= 0 else decl


def decl_span(lines: list[str], name: str) -> tuple[int, int] | None:
    """Return the 1-based line span of declaration name, docstring included."""
    pat = re.compile(rf"^{MODIFIERS}({KEYWORDS})\s+{re.escape(name)}(\s|$|:)")
    hits = [i for i, line in enumerate(lines) if pat.match(line)]
    if len(hits) != 1:
        return None
    start = hits[0]
    end = start + 1
    while end < len(lines) and not TOP.match(lines[end]):
        end += 1
    end -= 1
    while end > start and not lines[end].strip():
        end -= 1
    return doc_start(lines, start) + 1, end + 1


def check_decl(
    label: str, path: pathlib.Path, lines: list[str], span: tuple[int, int]
) -> str:
    """Return an error if a Lean-identifier link misses its declaration."""
    if path.suffix != ".lean" or not label.startswith(BT) or " " in label:
        return ""
    name = label.strip(BT).split(".")[-1]
    decl = decl_span(lines, name)
    if decl is None:
        return f"no unique declaration {name} in {path.name}"
    if not decl[0] <= span[0] <= span[1] <= decl[1]:
        return f"{name}: L{span[0]}-L{span[1]} outside L{decl[0]}-L{decl[1]}"
    return ""


def check_link(label: str, target: str, first: str | None, last: str | None) -> str:
    """Return an error message for one link, or an empty string."""
    path = (HERE / target).resolve()
    if not path.exists():
        return f"missing file: {target}"
    a = to_int(first)
    if a is None:
        return ""
    b = to_int(last) if last else a
    lines = path.read_text(encoding="utf-8").split("\n")
    if b is None or not 1 <= a <= b <= len(lines):
        return f"bad range {target}#L{first}-L{last}"
    return check_decl(label, path, lines, (a, b))


def top_level_words(text: str) -> int:
    """Count words outside details blocks, code blocks, math and link targets."""
    top = re.sub(r"<details>.*?</details>", "", text, flags=re.DOTALL)
    top = re.sub(BT * 3 + r".*?" + BT * 3, "", top, flags=re.DOTALL)
    top = re.sub(r"\([^)]*\)|\$[^$]*\$", "", top)
    return len(re.findall(r"[A-Za-z]{2,}", top))


def main() -> int:
    """Run all checks; return the process exit code."""
    text = DOC.read_text(encoding="utf-8")
    links = [
        m.groups() for m in LINK.finditer(text) if not m.group(2).startswith("http")
    ]
    errors = [e for e in (check_link(*g) for g in links) if e]
    goals = re.findall(r"🎯 ((?:(?!<br>)[^|\n])+)", text)
    errors += [
        f"goal not 7 words: {g.strip()}" for g in goals if len(g.split()) != GOAL_WORDS
    ]
    spans = sum(1 for g in links if g[2])
    words = top_level_words(text)
    sys.stdout.write(
        f"line-span links: {spans}, goals: {len(goals)}, top-level words: {words}\n"
    )
    for err in errors:
        sys.stdout.write(f"ERROR {err}\n")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
