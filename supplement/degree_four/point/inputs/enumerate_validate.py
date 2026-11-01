#!/usr/bin/env python3
"""Exact finite grammar screen. No numerical activities or solver imports."""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Exact validation refuses optimized Python (-O/-OO)")

from collections import Counter, defaultdict
from itertools import combinations, product
from math import comb, factorial
from pathlib import Path
import gzip
import hashlib
import json
import time

HERE = Path(__file__).resolve().parent


def dump(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def partitions(stubs):
    """Set partitions of distinct stub indices, canonical first-appearance order."""
    blocks = []

    def rec(k):
        if k == len(stubs):
            yield tuple(tuple(B) for B in blocks)
            return
        for B in blocks:
            if len(B) < 4 and all(stubs[j] != stubs[k] for j in B):
                B.append(k)
                yield from rec(k + 1)
                B.pop()
        blocks.append([k])
        yield from rec(k + 1)
        blocks.pop()

    yield from rec(0)


def boundary_choices(T, degrees):
    for h in (2, 3, 4):
        if h < len(T) or (h > 2 and any(degrees[i - 1] > 2 for i in T)):
            continue
        for x, b in ((0, 0), (1, 0), (1, 1)):
            if h == len(T) and x != b:
                continue
            yield h, x, b


def masks_from_clauses(d, edges, blocks, choices):
    forbidden = [(1 << (d + 1)) - 1]
    for i in range(1, d + 1):
        if all(x for T, (_, x, _) in zip(blocks, choices) if i in T):
            U = {0, i}
            for a, b in edges:
                if a == i:
                    U.add(b)
                if b == i:
                    U.add(a)
            forbidden.append(sum(1 << j for j in U))
    for T, (_, _, b) in zip(blocks, choices):
        if b:
            forbidden.append(sum(1 << i for i in T))
    return tuple(S for S in range(1 << (d + 1))
                 if all(S & Q != Q for Q in forbidden))


def key_of(view):
    return (view['d'], tuple(view['degrees']), tuple(map(tuple, view['edges'])),
            tuple(map(tuple, view['parts'])), tuple(map(tuple, view['choices'])))


def all_views():
    for d in (2, 3, 4):
        pairs = tuple(combinations(range(1, d + 1), 2))
        for degrees in product((2, 3, 4), repeat=d):
            if d > 2 and any(h > 2 for h in degrees):
                continue
            for e in range(1 << len(pairs)):
                edges = tuple(p for j, p in enumerate(pairs) if e >> j & 1)
                if any(degrees[a - 1] > 2 and degrees[b - 1] > 2 for a, b in edges):
                    continue
                internal = [sum(i in edge for edge in edges) for i in range(1, d + 1)]
                remaining = [h - 1 - r for h, r in zip(degrees, internal)]
                if any(r < 0 for r in remaining):
                    continue
                stubs = tuple(i for i, r in enumerate(remaining, 1) for _ in range(r))
                for parts in partitions(stubs):
                    blocks = tuple(tuple(stubs[j] for j in B) for B in parts)
                    for choices in product(*(tuple(boundary_choices(T, degrees)) for T in blocks)):
                        masks = masks_from_clauses(d, edges, blocks, choices)
                        yield dict(d=d, degrees=degrees, edges=edges, stubs=stubs,
                                   parts=parts, blocks=blocks, choices=choices,
                                   masks=masks)


def closed_form_counts():
    """Independent weighted matching/Bell-polynomial count, no view generation."""
    degree_two_table = {}
    for p in (2, 3, 4):
        for q in range(p, 5):
            a, b = p - 1, q - 1
            singleton_p = 9 if p == 2 else 3
            singleton_q = 9 if q == 2 else 3
            paired = 8 if p == q == 2 else 2
            count = sum(comb(a, k) * comb(b, k) * factorial(k)
                        * singleton_p ** (a - k) * singleton_q ** (b - k)
                        * paired ** k for k in range(min(a, b) + 1))
            degree_two_table[f'{p}{q}'] = count
    root_two_noedge = sum(v * (1 if k[0] == k[1] else 2)
                          for k, v in degree_two_table.items())
    root_two_edge = 1 + 2 * 3 + 2 * 3**2
    # Weighted set partition counts with block weights a1=9,a2=8,a3=5,a4=2.
    a = [0, 9, 8, 5, 2]
    B = [1]
    for n in range(1, 5):
        B.append(sum(comb(n - 1, k - 1) * a[k] * B[n - k]
                     for k in range(1, n + 1)))
    high = {}
    for d in (3, 4):
        high[d] = [factorial(d) // (factorial(d - 2 * m) * 2**m * factorial(m))
                   * B[d - 2 * m] for m in range(d // 2 + 1)]
    return dict(root_two_noedge_table=degree_two_table,
                root_two_edge=root_two_edge, weighted_partition_counts=B,
                by_root={'2': root_two_noedge + root_two_edge,
                         '3': sum(high[3]), '4': sum(high[4])},
                high_matching_contributions={str(k): v for k, v in high.items()})


def enumerate_grammar():
    start = time.monotonic()
    keys = set()
    h = hashlib.sha256()
    by_d, by_tuple, by_edges, by_stubs, by_polynomial = (Counter() for _ in range(5))
    assignment_upper = 0
    feasible_assignments = 0
    with gzip.open(HERE / 'views.jsonl.gz', 'wt') as f:
        for view in all_views():
            key = key_of(view)
            assert key not in keys
            keys.add(key)
            raw = json.dumps(view, sort_keys=True, separators=(',', ':')) + '\n'
            h.update(raw.encode())
            f.write(raw)
            d = view['d']
            assert len(view['stubs']) <= (6 if d == 2 else d)
            assert view['masks'] and view['masks'][0] == 0
            by_d[str(d)] += 1
            by_tuple[f"{d}:" + ''.join(map(str, view['degrees']))] += 1
            by_edges[f"{d}:{len(view['edges'])}"] += 1
            by_stubs[f"{d}:{len(view['stubs'])}"] += 1
            poly = tuple(sum(S.bit_count() == k for S in view['masks']) for k in range(d + 2))
            by_polynomial[f'{d}:{poly}'] += 1
            assignment_upper += 1 << (d + 1)
            feasible_assignments += len(view['masks'])
    independent = closed_form_counts()
    assert dict(by_d) == independent['by_root']
    assert by_edges['2:1'] == independent['root_two_edge']
    for d in (3, 4):
        for m, value in enumerate(independent['high_matching_contributions'][str(d)]):
            assert by_edges[f'{d}:{m}'] == value
    summary = dict(all_passed=True, total=len(keys), by_root=dict(by_d),
                   by_degree_tuple=dict(by_tuple), by_neighbor_edge_count=dict(by_edges),
                   by_stub_count=dict(by_stubs), by_local_polynomial=dict(by_polynomial),
                   generation_sha256=h.hexdigest(), closed_form=independent,
                   all_inside_mask_checks=assignment_upper,
                   feasible_inside_assignments=feasible_assignments,
                   elapsed_seconds=time.monotonic() - start)
    dump('enumeration.json', summary)
    return keys, summary


def cycle(n):
    return n, {(i, (i + 1) % n) for i in range(n)}


def complete_bipartite(a, b):
    return a + b, {(i, a + j) for i in range(a) for j in range(b)}


def fixtures():
    out = {f'C{n}': cycle(n) for n in (3, 4, 5)}
    out.update({f'K2,{b}': complete_bipartite(2, b) for b in (3, 4)})
    out['two_triangles_one_root'] = (5, {(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)})
    out['triangle_and_square_one_root'] = (6, {(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 5), (5, 0)})
    edges = set()
    for w, (a, b) in enumerate(combinations(range(4), 2), 4):
        edges.update({(a, w), (w, b)})
    out['once_subdivision_K4'] = (10, edges)
    edges = {(0, w) for w in (2, 3, 4)} | {(1, w) for w in (2, 3, 4)}
    edges |= {(0, 5), (5, 6), (6, 1)}
    out['three_once_one_twice'] = (7, edges)
    edges = {(0, 2), (2, 1), (0, 3), (3, 1), (0, 4), (4, 5), (5, 0), (1, 6), (6, 7), (7, 1)}
    out['two_once_two_returning_triangles'] = (8, edges)
    return out


def actual_view(adj, v, E):
    nbs = sorted(adj[v])
    d = len(nbs)
    R = {v, *nbs}
    local = {u: i for i, u in enumerate([v] + nbs)}
    degrees = tuple(len(adj[u]) for u in nbs)
    edges = tuple((i + 1, j + 1) for i in range(d) for j in range(i + 1, d)
                  if nbs[j] in adj[nbs[i]])
    stubs, endpoints = [], []
    for i, u in enumerate(nbs, 1):
        for w in sorted(adj[u] - R):
            stubs.append(i)
            endpoints.append(w)
    index, parts, ws = {}, [], []
    for j, w in enumerate(endpoints):
        if w not in index:
            index[w] = len(parts)
            parts.append([])
            ws.append(w)
        parts[index[w]].append(j)
    blocks = tuple(tuple(stubs[j] for j in B) for B in parts)
    choices = tuple((len(adj[w]), int(w in E), int(({w} | (adj[w] - R)) <= E)) for w in ws)
    view = dict(d=d, degrees=degrees, edges=edges, stubs=tuple(stubs), parts=tuple(map(tuple, parts)),
                blocks=blocks, choices=choices)
    view['masks'] = masks_from_clauses(d, edges, blocks, choices)
    return view, local


def validate_actual(keys):
    reports = {}
    for name, (n, edge_list) in fixtures().items():
        adj = [set() for _ in range(n)]
        for a, b in edge_list:
            assert a != b
            adj[a].add(b)
            adj[b].add(a)
        assert all(2 <= len(A) <= 4 for A in adj)
        assert all(len(adj[a]) == 2 or len(adj[b]) == 2 for a, b in edge_list)
        closed = [{u} | A for u, A in enumerate(adj)]
        full_family = [S for S in range(1 << n)
                       if all(not all(S >> u & 1 for u in C) for C in closed)]
        whole_poly = [sum(S.bit_count() == k for S in full_family) for k in range(n + 1)]
        exposures = 0
        repeated = 0
        boundary_classes = Counter()
        for v in range(n):
            R = {v} | adj[v]
            outer = sorted(set(range(n)) - R)
            root_mass = [0] * (n + 1)
            for B in range(1 << len(outer)):
                E = {u for j, u in enumerate(outer) if B >> j & 1}
                actual = []
                inside = [v] + sorted(adj[v])
                for S in range(1 << len(inside)):
                    I = E | {u for j, u in enumerate(inside) if S >> j & 1}
                    if all(not C <= I for C in closed):
                        actual.append(S)
                if not actual:
                    continue
                view, local = actual_view(adj, v, E)
                assert key_of(view) in keys
                assert tuple(actual) == view['masks']
                exposures += 1
                if len(set(view['blocks'])) < len(view['blocks']):
                    repeated += 1
                for T, (h, x, b) in zip(view['blocks'], view['choices']):
                    boundary_classes[f'{h},{len(T)},{x},{b}'] += 1
                for S in actual:
                    root_mass[len(E) + S.bit_count()] += 1
            assert root_mass == whole_poly
        reports[name] = dict(n=n, edges=sorted(sorted(e) for e in edge_list),
                             degrees=list(map(len, adj)), feasible_full_sets=len(full_family),
                             complement_polynomial=whole_poly, exposures=exposures,
                             repeated_boundary_exposures=repeated,
                             observed_boundary_classes=dict(boundary_classes))
    # Deliberate invalid-grammar controls, independent of every positive result.
    assert ((0, 1),) not in set(partitions((1, 1)))
    assert tuple(boundary_choices((1, 2, 3), (3, 2, 2))) == ()
    assert (2, 1, 0) not in set(boundary_choices((1, 2), (2, 2)))
    assert all(not (v['d'] > 2 and any(h > 2 for h in v['degrees'])) for v in all_views())
    summary = dict(all_passed=True, fixtures=reports,
                   total_exposures=sum(r['exposures'] for r in reports.values()),
                   negative_controls=['repeated neighbor in block', 'degree below incidence',
                                      'high-high edge', 'exhausted boundary flags disagree'])
    dump('validation.json', summary)
    return summary


if __name__ == '__main__':
    keys, enum = enumerate_grammar()
    val = validate_actual(keys)
    print(json.dumps(dict(total=enum['total'], by_root=enum['by_root'],
                          generation_sha256=enum['generation_sha256'],
                          all_inside_mask_checks=enum['all_inside_mask_checks'],
                          feasible_inside_assignments=enum['feasible_inside_assignments'],
                          fixture_exposures=val['total_exposures'], all_passed=True), indent=2))
