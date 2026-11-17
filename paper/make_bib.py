"""Generate paper/refs.bib from FINAL_SPEC.md §17 (single source of truth).

Each numbered reference line ``N. **[KEY]** text`` (bold optional) becomes
``@misc{KEY, note = {text}}`` so that bibtex (unsrtnat) prints exactly the
recorded bibliographic text, in order of first citation, and only for cited
keys.  Markdown italics become \\emph{}, bare URLs become \\url{}.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "docs" / "design" / "FINAL_SPEC.md"
OUT = ROOT / "paper" / "refs.bib"


def latexify(s: str) -> str:
    urls: list[str] = []

    def keep_url(m: re.Match) -> str:
        urls.append(m.group(0).rstrip(".,;)"))
        tail = m.group(0)[len(urls[-1]):]
        return f"\x00{len(urls) - 1}\x00{tail}"

    s = re.sub(r"https?://\S+", keep_url, s)
    s = s.replace("\\", "\\textbackslash{}")
    for ch in "&%$#_":
        s = s.replace(ch, "\\" + ch)
    s = re.sub(r"\*([^*]+)\*", r"\\emph{\1}", s)
    s = re.sub(r"`([^`]+)`", r"\\texttt{\1}", s)
    s = s.replace("–", "--").replace("—", "---").replace("’", "'")
    s = s.replace("“", "``").replace("”", "''")
    s = re.sub(r'"([^"]*)"', r"``\1''", s)
    s = s.replace("×", "$\\times$").replace("→", "$\\rightarrow$")
    s = re.sub(r"\x00(\d+)\x00", lambda m: "\\url{" + urls[int(m.group(1))] + "}", s)
    return s


def main() -> None:
    text = SPEC.read_text(encoding="utf-8")
    sec = text.split("## 17. References to cite", 1)[1].split("\n---", 1)[0]
    lines = sec.splitlines()
    entries: list[tuple[str, str]] = []
    cur: list[str] | None = None
    key = None
    for ln in lines:
        m = re.match(r"^\d+\.\s+(?:\*\*)?\[([A-Za-z0-9\-]+)\](?:\*\*)?\s+(.*)$", ln)
        if m:
            if cur is not None:
                entries.append((key, " ".join(cur)))
            key, cur = m.group(1), [m.group(2).strip()]
        elif cur is not None and ln.startswith(" ") and ln.strip():
            cur.append(ln.strip())
        elif cur is not None and not ln.strip():
            entries.append((key, " ".join(cur)))
            cur = None
    if cur is not None:
        entries.append((key, " ".join(cur)))
    out = []
    for k, body in entries:
        out.append("@misc{%s,\n  note = {%s}\n}\n" % (k, latexify(body)))
    OUT.write_text("\n".join(out), encoding="utf-8")
    print(f"wrote {len(entries)} entries to {OUT}")


if __name__ == "__main__":
    main()
