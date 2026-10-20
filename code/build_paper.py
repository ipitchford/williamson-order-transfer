"""Typeset one public manuscript without editing its scientific Markdown."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

p = argparse.ArgumentParser()
p.add_argument("repository", type=Path)
args = p.parse_args()
root = args.repository.resolve()
qa = root.parent / "evidence" / "pdf-qa"
qa.mkdir(parents=True, exist_ok=True)
source = (root / "paper.md").read_text()
first, body = source.split("\n", 1)
title = re.sub(r"^#+\s*", "", first).strip()
metadata = ["-M", "title=" + title, "-M", "author=Anonymous", "-M", "date=20 October 2026"]
if source.startswith("---\n"):
    body = source
    metadata = []
# Authorship is metadata, not a duplicate line under the title.
body = re.sub(r"^\s*(?:\*\*)?Anonymous\..*?\n", "\n", body, count=1)
# These source symbols occur in prose rather than mathematical delimiters.
body = body.replace("□", r"\(\square\)").replace("p≤6", r"\(p\leq6\)")
render = qa / "render.md"
render.write_text(body)
header = qa / "header.tex"
header.write_text(r"""
\usepackage{microtype}
\usepackage{xcolor}
\usepackage{needspace}
\definecolor{linkink}{HTML}{553653}
\AtBeginDocument{\hypersetup{colorlinks=true,linkcolor=linkink,urlcolor=linkink,citecolor=linkink}}
\setlength{\emergencystretch}{3em}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.42em}
\widowpenalty=10000
\clubpenalty=10000
\displaywidowpenalty=10000
\allowdisplaybreaks[1]
""")
layout = qa / "layout.lua"
layout.write_text(r"""
function Table(el)
  local n = 0
  for _, body in ipairs(el.bodies) do n = n + #body.body end
  return {pandoc.RawBlock("latex", "\\Needspace{" .. math.min(30, 3*n+4) .. "\\baselineskip}"), el}
end
""")
command = [
    "pandoc", str(render),
    "--from=markdown+tex_math_single_backslash+autolink_bare_uris",
    "--standalone", "--shift-heading-level-by=-1",
    "--resource-path=" + str(root),
    "--lua-filter=" + str(layout),
    "--include-in-header=" + str(header),
    *metadata,
    "-V", "documentclass=scrartcl", "-V", "papersize=a4",
    "-V", "fontsize=11pt", "-V", "geometry:top=23mm,bottom=23mm,left=24mm,right=24mm",
    "-V", "linestretch=1.04", "-o", str(root / "paper.tex"),
]
subprocess.run(command, check=True, cwd=root)
result = subprocess.run(["tectonic", "--keep-logs", str(root / "paper.tex")],
                        cwd=root, capture_output=True, text=True)
(qa / "build-output.log").write_text(result.stdout + result.stderr)
if result.returncode:
    raise SystemExit(result.returncode)
(qa / "build-receipt.json").write_text(json.dumps({
    "paper_md_sha256": hashlib.sha256(source.encode()).hexdigest(),
    "paper_pdf_sha256": hashlib.sha256((root/"paper.pdf").read_bytes()).hexdigest(),
    "command": command,
    "engine": "tectonic",
}, indent=2) + "\n")
print(root / "paper.pdf")
