# Rebuilding the manuscript

The scientific source is `paper.md`. The accompanying `paper.tex` and `paper.pdf` are generated from it.

Install Pandoc and Tectonic, then run from the repository root:

```sh
python3 code/build_paper.py .
```

The build uses bundled Latin Modern fonts and produces `paper.tex` and `paper.pdf`. It writes diagnostic material to `../evidence/pdf-qa/`; this material is not part of the manuscript. Font packages may be fetched by Tectonic on its first run. Compiling the supplied TeX directly with Tectonic is also supported.

The conversion normalises title/authorship metadata and prose mathematical glyphs for typesetting. It does not edit `paper.md` or run research computations.
