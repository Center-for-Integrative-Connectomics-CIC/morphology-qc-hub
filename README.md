# Morphology QC Hub

A browser tool for visually quality-checking labeled single-neuron reconstructions of
hippocampal neurons. It loads the labeled SWCs into the Allen mouse brain in 3D and lets
you inspect each neuron's compartment labels one at a time — catching an axon branch
tagged as dendrite or a soma in the wrong place — then flag status, add notes, and export.

It currently holds **2,622 labeled neurons** across six hippocampal areas:

| Area | Neurons |
|------|--------:|
| CA1  | 951 |
| DG   | 999 |
| PRE  | 328 |
| POST | 168 |
| PARA | 129 |
| CA2  | 47  |

## Run it

Live: https://center-for-integrative-connectomics-cic.github.io/morphology-qc-hub/

Or serve the folder locally:

```sh
python3 -m http.server 8000
```

## Using it

- **Pick a neuron** — click a soma in the 3D view, or a row in the left list. The full
  arbor loads, colour-coded by compartment (axon blue, basal red, apical yellow).
- **Compartment toggles** — hide/show axon, basal, or apical to inspect the labeling.
  Isolated mislabels (a node whose type disagrees with all its neighbours) are absorbed
  into the surrounding compartment for display, so they don't clutter the view — the
  underlying SWC is never modified.
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

Most areas are read as `<AREA>/**/*_labeled.swc`; CA1 sits in `hpf_ca1/` with a different
convention (`*_edited.swc`, deduplicated by neuron id). Coordinates are Allen CCFv3
microns; the build recentres them on the brain-shell centroid for display.

## Layout

```
index.html            the app (no framework, three.js via CDN)
build_hub.py          SWC -> binary baker
data/
  neurons.json        per-neuron metadata + QC stats
  brain.bin           Allen brain-shell mesh
  _root_dec.obj       source mesh for brain.bin
  swc/                per-neuron geometry (one .bin each)
  swc_src/            original SWCs for download (gitignored)
```
