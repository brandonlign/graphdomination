# Checks for the subdivision note

These files support the note *Scalar inequalities for subdivision counts*. Each folder holds the coefficients of one polynomial, and each `_check` folder holds a program that rebuilds that polynomial using exact integer arithmetic.

| Folder | Interval | Coefficients |
|---|---|---|
| `middle_interval/` | 1/4 ≤ q < 1/2 | 192 |
| `upper_interval/` | 1/2 ≤ q < 1 | 308 |

From this directory, run:

    python middle_interval_check/check_rational_independent.py
    python upper_interval_check/check_independent.py

The programs need only the Python standard library. Do not run Python with `-O`.
