# Morphology QC Hub

A browser tool for visually quality-checking labeled single-neuron reconstructions of
hippocampal neurons. It loads the labeled SWCs into the Allen mouse brain in 3D and lets
you inspect each neuron's compartment labels one at a time — catching an axon branch
tagged as dendrite or a soma in the wrong place — then flag status, add notes, and export.

It currently holds **9,007 labeled neurons** across eight areas of the hippocampal
formation:

| Area | Neurons |
|------|--------:|
| CA1  | 3,690 |
| DG   | 3,059 |
| CA3  | 896 |
| SUB  | 689 |
| PRE  | 329 |
| POST | 168 |
| PARA | 129 |
| CA2  | 47 |

## Run it

Live: https://center-for-integrative-connectomics-cic.github.io/morphology-qc-hub/

Or serve the folder locally:

```sh
python3 -m http.server 8000
```

## Using it

- **Pick a neuron** — click a soma in the 3D view, or a row in the left list. The full
  arbor loads, colour-coded by compartment (axon blue, basal red, apical yellow). Once a
  neuron is selected, **↑ / ↓** steps through the rest of its group without going back
  to the list.
- **Areas and subgroups** — each area in the list expands into named subgroups that can
  be opened one at a time, all at once, or neuron by neuron. SUB is grouped by the
  Saustad et al. scheme (column, layer, projection class or cell type) and coloured with
  the paper's own palette; the other areas split by reconstruction type.
- **Compartment toggles** — hide/show axon, basal, apical or soma, for the selected
  neuron and every loaded group at once. Isolated mislabels (a node whose type disagrees
  with all its neighbours) are absorbed into the surrounding compartment for display,
  and marker nodes follow the branch they annotate, so neither clutters the view — the
  underlying SWC is never modified.
- **CCF structures** — ten Allen CCFv3 surfaces (CA1, CA2, CA3, DG, SUB, ProS, POST,
  PRE, PAR, HATA), each shown or hidden on its own. The ◇ beside an area shows that
  area's surface by itself.
- **QC** — set pass / flag / redo and write notes per neuron (kept in the browser via
  `localStorage`). "Export QC" dumps everything to CSV; "SWC" downloads the original file.
- **Batches** — group neurons into a named collection to review, bulk-set status,
  export, or download together.

## Rebuilding the data

`build_hub.py` reads the labeled SWCs and writes the binaries the viewer fetches
(`data/neurons.json`, `data/brain.bin`, `data/swc/*.bin`) plus download copies in
`data/swc_src/`. Point it at wherever the labeled sets live:

```sh
LABELED_SETS=/path/to/labled_sets python3 build_hub.py
```

`build_ccf_meshes.py` bakes the CCF surfaces and the soma-to-structure lookup from the
Allen atlas, in the same display frame as the neurons:

```sh
python3 build_ccf_meshes.py         # needs brainglobe-atlasapi, trimesh, scikit-image
```

## Layout

```
index.html            the app (no framework, three.js vendored in lib/)
build_hub.py          SWC -> binary baker
build_ccf_meshes.py   Allen CCFv3 -> surfaces + soma lookup
data/
  neurons.json        per-neuron metadata + QC stats
  brain.bin           Allen brain-shell mesh
  _root_dec.obj       source mesh for brain.bin
  ccf_regions.bin     CCF surfaces in the viewer's display frame
  ccf_regions.json    per-structure colour, volume and byte offsets
  soma_ccf.json       neuron id -> CCF structure containing its soma
  id_migration.json   old id -> new id, so saved QC survives a rebuild
  swc/                per-neuron geometry (one .bin each)
  swc_src/            original SWCs for download (gitignored)
```
