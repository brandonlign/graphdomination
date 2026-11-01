# Reduced-class view enumeration

`enumerate_validate.py` lists every local view for graphs with degrees in {2,3,4} in which each edge has a degree-two endpoint. The rules are those of Section 2 of the article:

- The root has degree 2, 3 or 4, and every displayed edge has a degree-two endpoint.
- Remaining stubs are partitioned into blocks, one block for each boundary vertex.
- Each boundary vertex has a degree and flags (x, b) with values (0,0), (1,0) or (1,1).

`audit/grammar_count_independent.py` counts the same views by a separate method.
