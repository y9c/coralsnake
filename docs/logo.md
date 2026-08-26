---
layout: default
title: Logo
nav_order: 6
---

# Logo

`coralsnake logo` plots a DNA/RNA sequence logo from a set of motif sequences.
The scoring engine is pure numpy; the renderer needs matplotlib (the `plot`
extra).

## Usage

```bash
coralsnake logo -m ACGT -m ACGG -m CCGT -o logo.png
```

Per-motif weights from a file (`seq<TAB>count`):

```bash
cat motifs.tsv
# ACGT  5
# ACGG  3
coralsnake logo -i motifs.tsv -o logo.svg
```

### Options

| Option       | Description                                        |
| ------------ | -------------------------------------------------- |
| `-m, --motifs` | Motif sequence(s). Repeatable or comma-separated. |
| `-i, --input`  | Input file, one motif per line (`seq<TAB>count`). |
| `-o, --output` | Output image (png/svg). **required**.             |
| `-w, --weights`| Comma-separated per-motif weights.                |
| `--t2u/--no-t2u`| Convert T to U (default on).                     |
| `--2bit/--no-2bit`| Use 2-bit information logo (default on).        |
| `--normed`    | Normalize letter heights to sum to 1.              |

## Python API

```python
from coralsnake import Mlogo

m = Mlogo(motifs=["ACGT", "ACGG", "CCGT"], to2bit=True)
m.plot(ax)  # requires matplotlib (plot extra)
```
