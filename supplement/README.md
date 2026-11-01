# Exact proof checks

These files support the article *Mean order and counts of dominating sets in graphs of small degree* by Brandon Li. They contain the rational certificates, the finite view data and the programs that check them. The arguments that connect these checks to all graphs are in the article.

## Words used here

- **View:** what one chosen vertex can see of the graph around it.
- **Column:** a group of views that give exactly the same polynomials, so they only need to be checked once.
- **Certificate:** a list of exact numbers. If the checks on those numbers pass, the inequality in the paper follows.

## Run the checks

Use Python 3.12 and SymPy 1.14.0. From this directory, run:

    python replay_proof.py --destination NEW_DIRECTORY

The destination must not already exist, and Python must not be run with `-O`. The program checks every input against `PROOF_MANIFEST.sha256`, copies the files to the destination and runs 42 groups of checks. Results are written to `REPLAY_REPORT.json` and to one log per group. It takes about four minutes, needs no network connection, and blocks numerical optimizers. Some groups feed in invalid data on purpose and pass when that data is rejected.

## What each part checks

| Article result | Folder | What is checked |
|---|---|---|
| Certificate for graphs with degrees 2 and 3 (Proposition 4.1) | `interval_certificate/`, `examples/` | All 231,816 views in 9,867 columns, and 167,739 coefficient signs |
| Certificate for the reduced class (Proposition 4.2) | `degree_four/structural/`, `degree_four/interval/` | All 17,553 views in 1,722 columns, 29,274 coefficient signs and 162 multiplier signs |
| Equality for 0 < z < 1 (Section 5.2) | `equality/`, `degree_four/equality/` | The views with zero slack, listed exactly |
| Equality at z = 1 (Section 5.3) | `endpoint/`, `degree_four/equality/` | The endpoint certificate on 9,853 columns, and the 370 degree-four views |
| Cycles (Lemma 5.1) | `cycles/` | The recurrence on cycles of length 3 to 40 |
| Comparisons (Appendix A) | `degree_four/comparisons/` | The exact identities and integer comparisons |

Most results are checked by two separately written programs, and the programs compare their full output. The `audit/` and `independent/` folders hold the second programs.

The certificates for the two classes are `interval_certificate/subcubic_certificate.json` and `degree_four/interval/certificate.json`.

The programs verify the finite calculations; the general arguments are presented in the article.

Licence terms are in `LICENSE`.
