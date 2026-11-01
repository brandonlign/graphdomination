#!/usr/bin/env python3
"""Independent ordinary-set cycle exposure and LOW_BITS slack checks.

Standard library only; imports no local-view generator, certificate checker,
optimizer, or symbolic algebra library.  Enumerates actual graphs C3,...,C12,
not the full mixed-degree local relaxation.  The latter remains a dependency.
All polynomial arrays are in increasing power order.  Failures are explicit,
and optimized Python is refused as an additional defensive check.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path

if not __debug__:
    raise RuntimeError("Cycle exposure audit refuses optimized Python")

HERE = Path(__file__).resolve().parent
CERTIFICATE_SHA256 = "c6f20f9d56dd89013c8e1b186e1f0614368b7585c3243fa12ddc8b64fbfa0205"
ZERO_REPORT_SHA256 = "24ee52f4c265929a93d5dc2ac4b9f06dd3266d0e7cc72f869f55305b875d2e45"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def powerset(vertices):
    vertices = sorted(vertices)
    for size in range(len(vertices) + 1):
        for subset in combinations(vertices, size):
            yield frozenset(subset)


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return tuple(poly) or (0,)


def add(a, b):
    return trim((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                for i in range(max(len(a), len(b))))


def scale(poly, scalar):
    return trim(scalar * x for x in poly)


def multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def shift(poly, exponent):
    return trim([0] * exponent + list(poly))


def coefficient_counts(sets, selected=None):
    out = [0] * (1 + max(map(len, sets)))
    for subset in sets:
        if selected is None or selected in subset:
            out[len(subset)] += 1
    return trim(out)


def cycle(n):
    adjacency = {v: set() for v in range(n)}
    for v in range(n):
        w = (v + 1) % n
        adjacency[v].add(w)
        adjacency[w].add(v)
    require(all(len(neighbors) == 2 for neighbors in adjacency.values()), "Not a simple cycle")
    return adjacency


def feasible(adjacency, occupied):
    return all(not ({v} | neighbors) <= occupied for v, neighbors in adjacency.items())


def view(adjacency, root, exterior):
    neighbors = sorted(adjacency[root])
    labels = {u: i + 1 for i, u in enumerate(neighbors)}
    region = {root} | set(neighbors)
    boundary = set().union(*(adjacency[u] for u in neighbors)) - region
    blocks = []
    for w in boundary:
        incidence = tuple(sorted(labels[u] for u in neighbors if w in adjacency[u]))
        require(incidence, "Empty distance-two block")
        outside_closed = {w} | (adjacency[w] - region)
        blocks.append(dict(T=incidence, degree=len(adjacency[w]), x=int(w in exterior),
                           b=int(outside_closed <= exterior)))
    blocks.sort(key=lambda block: block["T"])
    return dict(d=len(neighbors), neighbor_degrees=[len(adjacency[u]) for u in neighbors],
                neighbor_edges=[(labels[u], labels[w]) for u, w in combinations(neighbors, 2)
                                if w in adjacency[u]], blocks=blocks)


def structurally_zero(record):
    if record["d"] != 2 or record["neighbor_degrees"] != [2, 2] or record["neighbor_edges"]:
        return False
    if all(block["x"] == block["b"] == 0 for block in record["blocks"]):
        return True
    blocks = record["blocks"]
    return (len(blocks) == 1 and tuple(blocks[0]["T"]) == (1, 2)
            and blocks[0]["degree"] == 2 and blocks[0]["x"] == blocks[0]["b"] == 1)


def local_rules(adjacency, root, exterior, inside):
    region = {root} | adjacency[root]
    if region <= inside:
        return False
    for u in adjacency[root]:
        closed = {u} | adjacency[u]
        if closed - region <= exterior and closed & region <= inside:
            return False
    boundary = set().union(*(adjacency[u] for u in adjacency[root])) - region
    for w in boundary:
        if {w} | (adjacency[w] - region) <= exterior and adjacency[w] & region <= inside:
            return False
    return True


def symbol(adjacency, vertex, occupied):
    return 2 * (len(adjacency[vertex]) - 2) + int(vertex in occupied)


def six_times_selected_transports(adjacency, root, occupied, selected):
    """Count actual oriented walks and stars, without a local-view engine."""
    values = Counter()
    a = symbol(adjacency, root, occupied)
    for u in adjacency[root]:
        b = symbol(adjacency, u, occupied)
        for w in adjacency[u] - {root}:
            c = symbol(adjacency, w, occupied)
            if a != c:
                name = f"R_{int(w in adjacency[root])}_{min(a, c)}_{b}_{max(a, c)}"
                if name in selected:
                    values[name] += 6 if a < c else -6
    for vertex, amount in [(root, 6)] + [(u, -3) for u in adjacency[root]]:
        neighbor_symbols = "".join(str(x) for x in sorted(
            symbol(adjacency, w, occupied) for w in adjacency[vertex]))
        name = f"S_2_{int(vertex in occupied)}_{neighbor_symbols}"
        if name in selected:
            values[name] += amount
    return values


def bernstein16(power):
    require(len(power) <= 17, "Unexpected slack degree")
    return tuple(sum(F(comb(k, j), comb(16, j)) * value
                     for j, value in enumerate(power) if j <= k) for k in range(17))


def main():
    path = HERE / "inputs/subcubic_certificate.json"
    require(hashlib.sha256(path.read_bytes()).hexdigest() == CERTIFICATE_SHA256, "Certificate hash mismatch")
    certificate = json.loads(path.read_text())
    require((certificate["a"], certificate["b"], certificate["degree"], certificate["bernstein_degree"])
            == ("0", "1", 2, 16), "Unexpected certificate format")
    support = certificate["support"]
    beta = [[F(x) for x in row] for row in certificate["bernstein_polynomials"]]
    require(len(support) == len(beta) == 74 and len(set(support)) == 74, "Support shape mismatch")
    require(all(len(row) == 3 for row in beta), "Polynomial degree mismatch")
    coefficients = {name: (b[0], 2 * (b[1] - b[0]), b[0] - 2 * b[1] + b[2])
                    for name, b in zip(support, beta)}
    inequalities = [name for name in support if name.startswith(("X_", "B_", "L_"))]
    selected = set(support) - set(inequalities)
    require(len(inequalities) == 9 and len(selected) == 65, "Unexpected support families")
    require(all(x >= 0 for name, row in zip(support, beta) if name in inequalities for x in row),
            "Negative inequality multiplier")
    require(all(name.startswith(("R_", "S_")) for name in selected), "Unexpected transport family")
    # All selected exterior rows have root degree 3. All selected association
    # patterns differ from (2,0), the pattern of every root of an ordinary cycle.
    require(all(name.split("_")[1] == "3" for name in inequalities if name.startswith(("X_", "B_"))),
            "Nonzero cycle exterior feature requires implementation")
    require("L_2_0" not in inequalities, "Nonzero cycle association feature requires implementation")

    zero_path = HERE / "inputs/zero_structure_full_replay.json"
    require(hashlib.sha256(zero_path.read_bytes()).hexdigest() == ZERO_REPORT_SHA256,
            "Full-universe zero-support report hash mismatch")
    source_zero = json.loads(zero_path.read_text())
    require(source_zero["all_passed"] and source_zero["certificate_sha256"] == CERTIFICATE_SHA256,
            "Unmatched full-universe zero-support dependency")
    require(source_zero["independently_regenerated_full_views"] == 231816
            and source_zero["independent_symbolic_columns"] == 9867
            and source_zero["zero_full_views"] == 7 and len(source_zero["zero_columns"]) == 6,
            "Unexpected full-universe zero-support summary")
    require(all(structurally_zero(record) for record in source_zero["zero_views"]),
            "Declared zero full view does not have the stated structure")

    cache = {}
    summaries = []
    for n in range(3, 13):
        adjacency = cycle(n)
        global_sets = [occupied for occupied in powerset(adjacency) if feasible(adjacency, occupied)]
        partition = coefficient_counts(global_sets)
        total_z = (0,)
        total_o = (0,)
        total_h = (0,)
        total_transport = {name: (0,) for name in selected}
        exposure_count = zero_count = 0
        witness = None
        for root in range(n):
            region = {root} | adjacency[root]
            groups = defaultdict(list)
            for occupied in global_sets:
                groups[occupied - region].append(occupied & region)
            for exterior, inside_sets in groups.items():
                require(len(set(inside_sets)) == len(inside_sets), "Duplicate inside set")
                actual = {inside for inside in powerset(region) if feasible(adjacency, exterior | inside)}
                local = {inside for inside in powerset(region) if local_rules(adjacency, root, exterior, inside)}
                require(actual == local == set(inside_sets), "Conditional family or local-rule mismatch")
                require(frozenset() in actual, "Missing empty completion")
                record = view(adjacency, root, exterior)
                zpoly = coefficient_counts(inside_sets)
                opoly = coefficient_counts(inside_sets, root)
                features = {name: [0] * 4 for name in selected}
                for inside in inside_sets:
                    raw = six_times_selected_transports(adjacency, root, exterior | inside, selected)
                    for name, value in raw.items():
                        features[name][len(inside)] += value
                features = {name: trim(poly) for name, poly in features.items()}
                key = json.dumps(record, sort_keys=True, separators=(",", ":"))
                if key not in cache:
                    target = scale(add(multiply((1, 4, 6), opoly), scale(multiply((0, 1, 3), zpoly), -1)), 6)
                    h = target
                    for name in selected:
                        h = add(h, scale(multiply(coefficients[name], features[name]), -1))
                    bernstein = bernstein16(h)
                    require(all(value >= 0 for value in bernstein), "Negative exact cycle slack coefficient")
                    require((not any(h)) == structurally_zero(record), "Cycle zero-support iff mismatch")
                    cache[key] = dict(view=record, Z=zpoly, O=opoly, transports=features,
                                      H=h, bernstein16=bernstein, occurrences=0)
                entry = cache[key]
                require((zpoly, opoly, features) == (entry["Z"], entry["O"], entry["transports"]),
                        "Local record lost information")
                entry["occurrences"] += 1
                exposure_count += 1
                zero_count += structurally_zero(record)
                total_z = add(total_z, shift(zpoly, len(exterior)))
                total_o = add(total_o, shift(opoly, len(exterior)))
                total_h = add(total_h, shift(entry["H"], len(exterior)))
                for name in selected:
                    total_transport[name] = add(total_transport[name], shift(features[name], len(exterior)))
                chosen = frozenset() if n == 3 else frozenset({2})
                if root == 0 and exterior == chosen:
                    witness = dict(root=root, exterior=sorted(exterior), view=record,
                                   feasible_inside=[sorted(inside) for inside in sorted(inside_sets, key=lambda s: (len(s), sorted(s)))],
                                   Z=zpoly, O=opoly, exposure_mass_numerator=shift(zpoly, len(exterior)),
                                   exposure_mass_denominator="n Z_Cn(z)", zero_slack=structurally_zero(record),
                                   cleared_slack_H=entry["H"], bernstein16=entry["bernstein16"])
        z_derivative = tuple(i * value for i, value in enumerate(partition))
        gap = add(multiply((1, 4, 6), z_derivative), scale(multiply((0, 1, 3), partition), -n))
        require(total_z == scale(partition, n), "Exposure masses do not sum to one")
        require(total_o == trim(z_derivative), "Root occupancy identity failed")
        require(all(not any(poly) for poly in total_transport.values()), "Averaged transport identity failed")
        require(total_h == scale(gap, 6), "Exact averaged slack identity failed")
        require(witness is not None, "Missing requested witness")
        if n == 4:
            require(exposure_count == zero_count and partition == (1, 4, 6), "C4 equality failed")
        else:
            require(not witness["zero_slack"] and any(witness["cleared_slack_H"]), "Witness is not strict in the interior")
        if n >= 5:
            require(witness["Z"] == (1, 3, 2) and witness["O"] == (0, 1, 1), "Singleton witness polynomial mismatch")
            require(any(tuple(block["T"]) == (1,) and (block["x"], block["b"]) == (1, 0)
                        for block in witness["view"]["blocks"]), "Missing singleton block with flag (1,0)")
        summaries.append(dict(cycle_length=n, feasible_global_sets=len(global_sets), partition=partition,
                              root_exterior_pairs=exposure_count, zero_slack_exposures=zero_count,
                              strict_interior_exposures=exposure_count-zero_count,
                              mean_gap_numerator=gap, mean_gap_denominator="n A(z) Z_Cn(z)",
                              all_65_selected_transport_identities_exact=True, witness=witness))

    views = []
    for key, entry in sorted(cache.items()):
        views.append({name: value for name, value in entry.items() if name != "transports"})
    output = dict(all_passed=True, certificate_sha256=CERTIFICATE_SHA256,
                  implementation="ordinary Python sets and exact rational arithmetic; no project imports",
                  scope="Actual cycle exposures for C3,...,C12 only; full 231816-view zero support is a pinned dependency, not rerun here",
                  optimizer_used=False, source_code_imports_used=False, evaluated_activities=[],
                  certificate_slack_checks="exact power polynomials and all 17 degree-16 Bernstein coefficients",
                  cycles=summaries, distinct_actual_cycle_views=len(views), local_views=views)
    def encode(value):
        if isinstance(value, F):
            return str(value)
        raise TypeError(type(value).__name__)
    (HERE / "CYCLE_EXPOSURES.json").write_text(json.dumps(output, default=encode, indent=2) + "\n")
    print(json.dumps(dict(all_passed=True, cycles=list(range(3, 13)), distinct_actual_cycle_views=len(views),
                          root_exterior_pairs=sum(row["root_exterior_pairs"] for row in summaries),
                          graph_transport_polynomial_identities=65 * len(summaries),
                          certificate_sha256=CERTIFICATE_SHA256), indent=2))


if __name__ == "__main__":
    main()
