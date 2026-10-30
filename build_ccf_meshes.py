#!/usr/bin/env python3
"""Bake the Allen CCFv3 hippocampal surfaces the QC hub draws, and assign each
soma to the structure that contains it.

Writes three files the viewer loads directly:

  data/ccf_regions.bin    float32 positions, float32 normals, uint32 indices,
                          concatenated per structure
  data/ccf_regions.json   one entry per structure with its colour, volume and
                          byte offsets into the .bin
  data/soma_ccf.json      neuron id -> CCF acronym containing its soma

Coordinates come out in the hub's display frame: the vertices go through
``build_hub.to_display`` so they land in the same space as the neurons with no
transform applied in the browser.

    annotation voxel (i,j,k) -> (i,j,k) * 25 um = CCF (AP, DV, ML) microns
                             -> build_hub.to_display -> hub display mm

The middle step is a plain scale with no axis reversal, which the soma check at
the end verifies end to end rather than assuming.

Needs brainglobe-atlasapi, trimesh, scikit-image and scipy; the atlas itself is
downloaded and cached by brainglobe on first run.
"""
import json
import os
import struct
import sys
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import trimesh
from scipy import ndimage
from skimage import measure

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
sys.path.insert(0, HERE)
import build_hub                                   # the hub's own transform

ATLAS = "allen_mouse_25um"
VOX_MM3 = (25 ** 3) / 1e9
FACE_BUDGET = 40000

# Every CCF structure the project's areas map onto, each drawn as its own mesh.
# Allen gives CA1/CA2/CA3/DG the same green and SUB a near-identical one, which
# cannot be told apart on screen, so the hub uses its own palette.  ProS and
# HATA carry no hub area of their own but frame the subiculum.
CCF_COLOR = {"CA1": "#4a90d9", "CA2": "#d96fa8", "CA3": "#5ec8c0",
             "DG": "#d9a441", "SUB": "#9b7fd4", "ProS": "#b06bff",
             "POST": "#35c46a", "PRE": "#2aa4e0", "PAR": "#e86fb0",
             "HATA": "#8a9bad"}
ORDER = ("CA1", "CA2", "CA3", "DG", "SUB", "ProS", "POST", "PRE", "PAR", "HATA")
# hub area -> CCF acronym (PARA is PAR in the Allen nomenclature)
AREA_TO_CCF = {"CA1": "CA1", "CA2": "CA2", "CA3": "CA3", "DG": "DG",
               "SUB": "SUB", "POST": "POST", "PRE": "PRE", "PARA": "PAR"}


def atlas():
    from brainglobe_atlasapi import BrainGlobeAtlas
    return BrainGlobeAtlas(ATLAS)


def has(a, acronym):
    try:
        a.structures[acronym]
        return True
    except KeyError:
        return False


def descendants(a, acronym):
    """Every structure id at or below an acronym - CA1 covers CA1so/sp/sr/slm."""
    root = a.structures[acronym]["id"]
    return np.fromiter((s["id"] for s in a.structures.values()
                        if root in s["structure_id_path"]), dtype=np.uint32)


def surface(mask):
    """Marching cubes -> drop specks -> smooth -> decimate -> display frame."""
    pad = np.pad(mask, 1).astype(np.uint8)
    v, f, _, _ = measure.marching_cubes(pad, 0.5, step_size=2)
    v = (v - 1.0) * 25.0                           # un-pad, voxels -> microns
    m = trimesh.Trimesh(vertices=v, faces=f, process=True, validate=True)
    if m.body_count > 1:
        # stray single-voxel components survive smoothing as spikes, so they go
        # first rather than after
        comps = m.split(only_watertight=False)
        if comps:
            big = max(len(c.faces) for c in comps)
            keep = [c for c in comps if len(c.faces) >= max(24, 0.02 * big)]
            if keep and len(keep) < len(comps):
                m = trimesh.util.concatenate(keep)
    try:
        m = trimesh.smoothing.filter_taubin(m, iterations=10)
    except Exception:
        pass
    if len(m.faces) > FACE_BUDGET:
        try:
            cand = m.simplify_quadric_decimation(face_count=FACE_BUDGET)
            cand.update_faces(cand.nondegenerate_faces())
            if len(cand.faces):
                m = cand
        except Exception:
            pass
    m.update_faces(m.nondegenerate_faces())
    m.remove_unreferenced_vertices()
    m.fix_normals()
    V = build_hub.to_display(np.asarray(m.vertices)).astype(np.float32)
    F = np.asarray(m.faces, dtype=np.uint32)
    N = np.asarray(trimesh.Trimesh(vertices=V, faces=F, process=False).vertex_normals,
                   dtype=np.float32)
    return V, N, F


def bake(a):
    entries, blobs, off = [], [], 0
    print(f"hub frame: CENTER {np.round(build_hub.CENTER, 1)} um, SCALE {build_hub.SCALE}\n")
    for acr in ORDER:
        if not has(a, acr):
            print(f"  {acr:5} not in {ATLAS}, skipped")
            continue
        m = np.isin(a.annotation, descendants(a, acr))
        if not m.any():
            continue
        V, N, F = surface(m)
        pos, nrm, idx = V.tobytes(), N.tobytes(), F.tobytes()
        blobs += [pos, nrm, idx]
        entries.append(dict(name=acr, kind="ccf", parent="CCF",
                            color=CCF_COLOR.get(acr, "#9aa7b4"),
                            voxels=int(m.sum()), mm3=round(int(m.sum()) * VOX_MM3, 4),
                            nv=int(len(V)), nt=int(len(F)), posOff=off,
                            normOff=off + len(pos), idxOff=off + len(pos) + len(nrm)))
        off += len(pos) + len(nrm) + len(idx)
        print(f"  {acr:5} {len(F):7,d} tris  {entries[-1]['mm3']:7.3f} mm3")

    with open(os.path.join(DATA, "ccf_regions.bin"), "wb") as fh:
        for b in blobs:
            fh.write(b)
    allV = np.concatenate([np.frombuffer(blobs[i * 3], dtype=np.float32).reshape(-1, 3)
                           for i in range(len(entries))])
    json.dump(dict(regions=entries,
                   bbox=dict(min=allV.min(0).tolist(), max=allV.max(0).tolist()),
                   frame="hub display (build_hub.to_display)"),
              open(os.path.join(DATA, "ccf_regions.json"), "w"))
    print(f"\n{len(entries)} surfaces, {sum(e['nt'] for e in entries):,} tris, "
          f"{off / 1e6:.1f} MB -> data/ccf_regions.bin")
    return entries, allV


def soma_structures(a):
    """id -> CCF acronym containing that soma, read from the source SWCs."""
    hub = json.load(open(os.path.join(DATA, "neurons.json")))["neurons"]
    res = np.asarray(a.resolution, dtype=float)
    shape = np.asarray(a.annotation.shape)
    sets = {acr: set(descendants(a, acr).tolist()) for acr in ORDER if has(a, acr)}

    ids, pts = [], []
    for nid, n in hub.items():
        src = os.path.join(HERE, n.get("src", ""))
        if not n.get("src") or not os.path.exists(src):
            continue
        xyz = _soma_um(src)
        if xyz is None:
            continue
        ids.append(nid)
        pts.append(xyz)
    if not ids:
        print("\nno source SWCs on disk, skipping data/soma_ccf.json")
        return None

    V = np.rint(np.asarray(pts) / res).astype(np.int64)
    ok = np.all((V >= 0) & (V < shape), axis=1)
    out = {}
    sid = np.zeros(len(V), dtype=np.uint32)
    sid[ok] = a.annotation[V[ok, 0], V[ok, 1], V[ok, 2]]
    for nid, s in zip(ids, sid.tolist()):
        for acr, members in sets.items():
            if s in members:
                out[nid] = acr
                break
    json.dump(dict(soma=out), open(os.path.join(DATA, "soma_ccf.json"), "w"))

    print(f"\nsoma -> CCF structure: {len(out):,} of {len(ids):,} somas placed")
    hits = tot = 0
    for nid in ids:
        want = AREA_TO_CCF.get(hub[nid]["area"])
        if want:
            tot += 1
            hits += out.get(nid) == want
    print(f"agreement with the curator's area label: {hits / max(tot, 1):.1%} "
          f"({hits:,}/{tot:,})")
    return out


def _soma_um(swc_path):
    """First type-1 node as (AP, DV, ML) microns, else the root."""
    root = None
    with open(swc_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 7:
                continue
            xyz = (float(p[2]), float(p[3]), float(p[4]))
            if int(float(p[1])) == 1:
                return xyz
            if root is None and int(float(p[6])) == -1:
                root = xyz
    return root


def main():
    a = atlas()
    entries, allV = bake(a)
    soma_structures(a)

    # the surfaces have to land where the neurons already are
    hub = json.load(open(os.path.join(DATA, "neurons.json")))["neurons"]
    stored = np.array([n["soma"] for n in hub.values()], dtype=float)
    lo, hi = allV.min(0) - 0.5, allV.max(0) + 0.5
    inb = np.all((stored >= lo) & (stored <= hi), axis=1)
    print(f"somas inside the surface bbox: {inb.mean():.1%}")


if __name__ == "__main__":
    main()
