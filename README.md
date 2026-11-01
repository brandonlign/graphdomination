# Mean order and counts of dominating sets in graphs of small degree

Brandon Li

A dominating set in a graph is a set of vertices such that every vertex is either in the set or next to a vertex in it. This project asks how few dominating sets a graph can have when its vertices have degree 2, 3 or 4, and how large those sets are on average. In both questions, disjoint 4-cycles turn out to be the extreme case.

The repository contains the paper and the code that checks its computer-assisted steps.

## Where to start

1. Read the introduction of the [paper](docs/manuscript.pdf). Section 1.2 outlines the proof.
2. To check the computations, follow "Run the checks" below. The [verification guide](docs/verification_guide.pdf) explains the data files.
3. For the separate result on subdivided graphs, see the [note on subdivision counts](docs/subdivision_supplement.pdf).

## Main result

Let D_G(t) be the domination polynomial of G, and let n be its number of vertices.

- If every degree of G is 2, 3 or 4, then D_G(t) ≥ (t⁴ + 4t³ + 6t²)^(n/4) for every t ≥ 1. At t = 1, this says that G has at least 11^(n/4) dominating sets.
- The mean size of a dominating set, weighted by t, is at most n(t² + 3t + 3)/(t² + 4t + 6) when all degrees are 2 or 3, and when all degrees are 2, 3 or 4 and every edge has a degree-two endpoint.
- Equality holds only for disjoint unions of 4-cycles.

## Run the checks

Use Python 3.12 and SymPy 1.14.0. From the repository root:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python supplement/replay_proof.py --destination .replay
```

For the subdivision note:

```sh
cd subdivision_checks
python middle_interval_check/check_rational_independent.py
python upper_interval_check/check_independent.py
```

[`supplement/README.md`](supplement/README.md) explains which folder checks which result.

## Build the PDFs

With pdfLaTeX installed:

```sh
sh manuscript/build.sh
```

## Contents

| Folder | Contents |
|---|---|
| `docs/` | The three PDFs |
| `manuscript/` | LaTeX sources and build script |
| `supplement/` | Certificates, view data and checking programs for the main theorem |
| `subdivision_checks/` | Coefficients and checking programs for the subdivision note |

`SHA256SUMS` lists a checksum for every file. An archived copy is on [Zenodo](https://doi.org/10.5281/zenodo.23844235).

Code and data are licensed as described in [`supplement/LICENSE`](supplement/LICENSE).
