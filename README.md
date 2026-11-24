# ComponentSearch2

Graph compaction and component identification for
[MATAM](https://github.com/xvazquezc/matam). This is a fork of
[bonsai-team/ComponentSearch2](https://github.com/bonsai-team/ComponentSearch2)
(f705924).

It reads the read overlap graph written by `ovgraphbuild`, groups reads into
meta nodes, filters and splices the meta graph, and writes the components
that MATAM assembles one by one.

## Usage

```
componentsearch -n nodes.csv -e edges.csv -b basename [-N INT] [-E INT] [-r SEED]
```

| Option | Meaning | Default |
|---|---|---|
| `-n` | nodes CSV (`;`-separated, one header line, one read per line) | `nodes.csv` |
| `-e` | edges CSV (`;`-separated; `Source` and `Target` columns found by name) | `edges.csv` |
| `-b` | basename of the outputs | empty |
| `-N` | meta nodes with at most this many reads are removed | 0 |
| `-E` | meta edges with at most this many read links are removed (0: none) | 0 |
| `-r` | seed of the random BFS starts | current time |
| `-o` | accepted, does nothing (as upstream) | |

Outputs: `basename.components.csv` (the nodes file plus `metanode` and
`component` columns, `-1` for reads filtered out), `basename.metaNodes.csv`
and `basename.metaEdges.csv`.

Build with `make` (C++11).

## Differences from upstream

- **Memory and time.** The read graph is held once, as a compressed adjacency
  array (8 bytes per edge instead of about 50), and the meta graph is changed
  in place. For a 294 M-edge graph: 14.8 GB and 9 min upstream, 2.3 GB and
  2 min here.
- **Same output without the edge filter.** With `-E 0`, the outputs are byte
  for byte those of upstream with its edge filter off, for the same seed
  (glibc `rand`); the BFS start is still the k-th remaining read for
  `rand() % remaining`. Upstream's edge filter rarely removed anything: on
  two production graphs at seed 42, `-E 0` here and upstream `-N 1 -E 1`
  give identical files.
- **Edge filter fixed.** Upstream passed `-N` instead of `-E` to the edge
  filter, counted the links of a meta node per outside read instead of per
  neighboring meta node, and removed the neighbor whose meta node number
  equalled such a read number, on one side only. On production graphs it
  removed none of the edges of weight 1. Here the weight of an edge is the
  number of read links between the two meta nodes, and edges of weight at
  most `-E` are removed from both sides.
- **Input errors.** An edge to a read missing from the nodes file is an error
  (exit code 1) instead of undefined behaviour.
