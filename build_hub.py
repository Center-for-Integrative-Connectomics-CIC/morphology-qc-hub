#!/usr/bin/env python3
"""Bake labeled single-neuron SWC reconstructions into the compact binaries the
viewer (index.html) fetches. Run once whenever the labeled sets change.

    LABELED_SETS=/path/to/labled_sets python3 build_hub.py

Areas: CA1, CA2, DG, POST, PRE, PARA. Most use "<AREA>/**/*_labeled.swc"; CA1 lives
in "hpf_ca1/" and its curated labels are "*_edited.swc" (deduplicated by neuron id).

Coordinates are Allen CCFv3 microns (AP, DV, ML). Display space is X=ML, Y=-DV (up),
Z=AP, centred on the brain-shell centroid and scaled to mm. SWC type column:
1 soma, 2 axon, 3 basal dendrite, 4 apical, 0 undefined, 10/11/12 custom.
"""
import os, glob, json, struct
import numpy as np
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
LSET = os.environ.get("LABELED_SETS",
                      os.path.join(HERE, "..", "clarity-template", "labled_sets"))
DATA = os.path.join(HERE, "data")
SWCOUT = os.path.join(DATA, "swc")
SRCOUT = os.path.join(DATA, "swc_src")

AREAS = ["CA1", "CA2", "DG", "POST", "PRE", "PARA"]
AREA_COLOR = {"CA1": "#f2762e", "CA2": "#e0403a", "DG": "#f2b705",
              "POST": "#2e9e57", "PRE": "#2aa4e0", "PARA": "#c840c8"}
AREA_DIR = {"CA1": "hpf_ca1"}
AREA_GLOB = {"CA1": "*_edited.swc"}
TYPE_NAME = {1: "soma", 2: "axon", 3: "dendrite", 4: "apical", 0: "undef", 10: "c10", 11: "c11", 12: "c12"}
EPS_UM = 6.0


def load_obj(fn):
    V, F = [], []
    for L in open(fn):
        if L.startswith("v "): V.append([float(x) for x in L.split()[1:4]])
        elif L.startswith("f "):
            idx = [int(p.split("/")[0]) - 1 for p in L.split()[1:]]
            for i in range(1, len(idx) - 1): F.append([idx[0], idx[i], idx[i + 1]])
    return np.array(V, float), np.array(F, np.int64)


_root_V, _root_F = load_obj(os.path.join(DATA, "_root_dec.obj"))
CENTER = _root_V.mean(0)
SCALE = 1000.0


def to_display(P):
    P = np.asarray(P, float).reshape(-1, 3)
    ap, dv, ml = P[:, 0], P[:, 1], P[:, 2]
    return np.column_stack([(ml - CENTER[2]) / SCALE, -(dv - CENTER[1]) / SCALE, (ap - CENTER[0]) / SCALE])


def parse_swc(fn):
    ids, typ, xyz, par = [], [], [], []
    for L in open(fn):
        if L.startswith("#"): continue
        p = L.split()
        if len(p) < 7: continue
        ids.append(int(float(p[0]))); typ.append(int(float(p[1])))
        xyz.append((float(p[2]), float(p[3]), float(p[4]))); par.append(int(float(p[6])))
    return np.array(ids), np.array(typ), np.array(xyz, float), np.array(par)


def rdp(points, eps):
    """Ramer-Douglas-Peucker polyline simplification, iterative."""
    n = len(points)
    if n < 3: return list(range(n))
    keep = [False] * n; keep[0] = keep[-1] = True; stack = [(0, n - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1: continue
        p0, p1 = points[a], points[b]; d = p1 - p0; L = np.linalg.norm(d); seg = points[a + 1:b] - p0
        dist = np.linalg.norm(seg, axis=1) if L < 1e-9 else np.linalg.norm(np.cross(seg, d / L), axis=1)
        k = int(np.argmax(dist))
        if dist[k] > eps:
            idx = a + 1 + k; keep[idx] = True; stack.append((a, idx)); stack.append((idx, b))
    return [i for i in range(n) if keep[i]]


def process(fn):
    ids, typ, xyz, par = parse_swc(fn)
    id2i = {int(i): k for k, i in enumerate(ids)}
    children = defaultdict(list); roots = []
    for k, pa in enumerate(par):
        if pa == -1 or int(pa) not in id2i: roots.append(k)
        else: children[id2i[int(pa)]].append(k)

    tcount = {int(t): int(c) for t, c in zip(*np.unique(typ, return_counts=True))}
    seglen = 0.0; axonlen = 0.0
    for k, pa in enumerate(par):
        if pa == -1 or int(pa) not in id2i: continue
        d = np.linalg.norm(xyz[k] - xyz[id2i[int(pa)]]); seglen += d
        if typ[k] == 2: axonlen += d

    pts = []; ptyp = []; edges = []; cache = {}

    def pidx(k):
        if k not in cache: cache[k] = len(pts); pts.append(xyz[k]); ptyp.append(int(typ[k]))
        return cache[k]

    for r in roots:
        stack = [r]
        while stack:
            b = stack.pop()
            for c in children[b]:
                path = [b, c]; cur = c
                while len(children[cur]) == 1: path.append(children[cur][0]); cur = children[cur][0]
                poly = np.array([xyz[k] for k in path]); kept = [path[i] for i in rdp(poly, EPS_UM)]
                for i in range(len(kept) - 1): edges.append((pidx(kept[i]), pidx(kept[i + 1])))
                if len(children[cur]) >= 1: stack.append(cur)

    P = to_display(np.array(pts)).astype(np.float32)
    E = np.array(edges, np.uint32)
    T = np.array(ptyp, np.uint8)
    soma_raw = xyz[roots[0]] if roots else xyz[0]
    return P, E, T, to_display([soma_raw])[0], dict(
        n_nodes=int(len(ids)), types={TYPE_NAME.get(t, str(t)): c for t, c in tcount.items()},
        len_um=round(float(seglen), 1), axon_um=round(float(axonlen), 1),
        n_branches=int(sum(1 for k in range(len(ids)) if len(children[k]) > 1)),
        n_tips=int(sum(1 for k in range(len(ids)) if len(children[k]) == 0)),
        n_roots=int(len(roots)))


def write_mesh(fn, P, E):
    # header <II (n_points, n_tris) + Float32 points + Uint32 triangle index (stride 3)
    with open(fn, "wb") as o:
        o.write(struct.pack("<II", len(P), len(E))); o.write(P.tobytes()); o.write(E.tobytes())


def write_arbor(fn, P, E, T):
    # header <II (n_points, n_edges) + Float32 points + Uint32 line index (stride 2) + Uint8 SWC type per point
    with open(fn, "wb") as o:
        o.write(struct.pack("<II", len(P), len(E))); o.write(P.tobytes()); o.write(E.tobytes()); o.write(T.tobytes())


def link_source(fn, nid):
    """Expose the original SWC for download under data/swc_src/ (hardlink, copy fallback)."""
    dst = os.path.join(SRCOUT, nid + ".swc")
    if os.path.exists(dst): return
    try:
        os.link(fn, dst)
    except OSError:
        import shutil; shutil.copy2(fn, dst)


def main():
    os.makedirs(SWCOUT, exist_ok=True)
    os.makedirs(SRCOUT, exist_ok=True)
    write_mesh(os.path.join(DATA, "brain.bin"),
               to_display(_root_V).astype(np.float32), _root_F.astype(np.uint32))

    neurons = {}; tot = 0
    for area in AREAS:
        folder = AREA_DIR.get(area, area); pattern = AREA_GLOB.get(area, "*_labeled.swc")
        files = sorted(glob.glob(os.path.join(LSET, folder, "**", pattern), recursive=True))
        print(f"[{area}] {len(files)} {pattern} in {folder}/")
        for i, fn in enumerate(files):
            nid = os.path.basename(fn)[:-4]
            if nid in neurons:
                if neurons[nid]["area"] != area: print("  collision across areas:", nid)
                continue
            try:
                P, E, T, soma, st = process(fn)
            except Exception as e:
                print("  fail", nid, repr(e)); continue
            write_arbor(os.path.join(SWCOUT, nid + ".bin"), P, E, T)
            link_source(fn, nid)
            tot += os.path.getsize(os.path.join(SWCOUT, nid + ".bin"))
            neurons[nid] = dict(id=nid, area=area, set=os.path.basename(os.path.dirname(fn)),
                                src=f"data/swc_src/{nid}.swc",
                                soma=[round(float(v), 4) for v in soma],
                                np=int(len(P)), ne=int(len(E)), **st)
            if i % 100 == 0: print(f"   {i}/{len(files)} {nid} ({len(P)} pts)")

    out = dict(areaColor=AREA_COLOR, areas=AREAS,
               counts={a: sum(1 for n in neurons.values() if n["area"] == a) for a in AREAS},
               neurons=neurons)
    json.dump(out, open(os.path.join(DATA, "neurons.json"), "w"))
    print(f"\nwrote {len(neurons)} neuron bins ({tot/1e6:.1f} MB), brain.bin, neurons.json")
    print("counts:", out["counts"], "total", len(neurons))


if __name__ == "__main__":
    main()
