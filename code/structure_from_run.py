#!/usr/bin/env python
"""
structure_from_run.py

Turn one run of the forward generator into a real chemical structure.

The generator (mini_forward_generator.py) only tracks switches and counts bonds; it never
builds atoms. This script runs it unchanged, records every coupling in the order it fired
(donor first, as _fire_coupling and _nucleate_one pass them to _bond), and then assembles
the molecule with RDKit from real monolignols (coniferyl G, sinapyl S, p-coumaryl H),
applying the same edits the rulebook describes:

  beta-O-4 : donor Cb-O4 acceptor; donor Ca=Cb -> single; water adds OH at donor Ca
  beta-5   : donor Cb-C5 acceptor; acceptor O4 closes on donor Ca (phenylcoumaran)
  beta-1   : donor Cb-C1 acceptor; acceptor ring -> cyclohexadienone (C4=O); the acceptor's
             alpha-OH (from its own earlier beta-O-4) closes on donor Ca (spirodienone)
  beta-beta: Cb-Cb; each gamma-O closes on the partner's Ca (resinol)
  5-5      : C5-C5 (biphenyl)
  4-O-5    : first unit's O4 - second unit's C5 (diaryl ether)

Units whose beta was never used (joined only through C5 or O4, or never joined) keep their
cinnamyl-alcohol side chain, as the generator's summary() reports.

Outputs (results/structures/):
  <tag>.sdf        full structure with hydrogens and 3D coordinates (ETKDGv3 + MMFF94)
  <tag>.json       timeline for the tutorial viewer: atoms with 2D/3D coordinates and the event at
                   which each appears, bonds with their order history, units, events
  <tag>.log        run summary

Stereocentres are not assigned; the 3D model is one diastereomer, a geometry illustration.

Usage (RDKit needed; the lignin env's RDKit is broken by NumPy 2, use base conda):
    ~/miniconda3/bin/python code/structure_from_run.py --sg 1.0 --seed 7000 \
        --params results/fitted_params_ref_v11e_sw.json
"""
import argparse
import json
import os
import sys
import time

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdDepictor, rdMolDescriptors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mini_forward_generator as mfg          # noqa: E402
from repo_paths import rel                    # noqa: E402

OUT = os.path.join(HERE, os.pardir, "results", "structures")
SMI = {"G": "OC/C=C/c1ccc(O)c(OC)c1", "S": "OC/C=C/c1cc(OC)c(O)c(OC)c1", "H": "OC/C=C/c1ccc(O)cc1"}


class RecordingKMC(mfg.MiniForwardKMC):
    """The generator, unchanged, plus a record of units and couplings in firing order."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.events = []          # (kind, donor uid, acceptor uid)
        self.unit_step = {}       # uid -> number of events already fired when the unit appeared

    def _new_monomer(self):
        u = super()._new_monomer()
        self.unit_step[u.uid] = len(self.events)
        return u

    def _bond(self, kind, a, b):
        super()._bond(kind, a, b)
        self.events.append((kind, a.uid, b.uid))


def positions(m):
    ring = m.GetRingInfo().AtomRings()[0]
    c1 = 4
    nb = sorted(n.GetIdx() for n in m.GetAtomWithIdx(c1).GetNeighbors() if n.GetIdx() in ring)
    order = [c1, nb[0]]
    while len(order) < 6:
        order.append([n.GetIdx() for n in m.GetAtomWithIdx(order[-1]).GetNeighbors()
                      if n.GetIdx() in ring and n.GetIdx() != order[-2]][0])
    c1, c6, c5, c4, c3, c2 = order
    o4 = [n.GetIdx() for n in m.GetAtomWithIdx(c4).GetNeighbors() if n.GetSymbol() == "O"][0]
    return dict(Og=0, Cg=1, Cb=2, Ca=3, C1=c1, C2=c2, C3=c3, C4=c4, C5=c5, C6=c6, O4=o4)


def build(gen):
    """Assemble the RDKit molecule. Returns (mol, atom_step, atom_unit, bond_hist)."""
    rw = Chem.RWMol()
    P, atom_step, atom_unit = {}, {}, {}
    bond_hist = {}                                    # (i, j) sorted -> [[step, order], ...]
    bond_kind = {}                                    # (i, j) sorted -> coupling that made it
    uids = sorted(gen.units)
    for uid in uids:
        u = gen.units[uid]
        if u.mtype not in SMI:
            raise SystemExit("unit type %s not supported (grass extension off only)" % u.mtype)
        m = Chem.MolFromSmiles(SMI[u.mtype])
        Chem.Kekulize(m, clearAromaticFlags=True)
        pos = positions(m)
        off = rw.GetNumAtoms()
        for a in m.GetAtoms():
            na = Chem.Atom(a.GetAtomicNum())
            i = rw.AddAtom(na)
            atom_step[i] = gen.unit_step[uid]; atom_unit[i] = uid
        for b in m.GetBonds():
            i, j = b.GetBeginAtomIdx() + off, b.GetEndAtomIdx() + off
            rw.AddBond(i, j, b.GetBondType())
            bond_hist[tuple(sorted((i, j)))] = [[gen.unit_step[uid], 2 if b.GetBondType() == Chem.BondType.DOUBLE else 1]]
        P[uid] = {k: v + off for k, v in pos.items()}
    alpha_O = {}                                      # uid -> index of its alpha-OH oxygen

    def setb(i, j, order, step):
        b = rw.GetBondBetweenAtoms(i, j)
        b.SetBondType(Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
        bond_hist[tuple(sorted((i, j)))].append([step, order])

    def addb(i, j, step, order=1):
        rw.AddBond(i, j, Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
        bond_hist[tuple(sorted((i, j)))] = [[step, order]]
        bond_kind[tuple(sorted((i, j)))] = gen.events[step - 1][0]

    for step, (kind, d, a) in enumerate(gen.events, start=1):
        D, A = P[d], P[a]
        if kind in ("bo4", "b5", "b1", "bb"):
            setb(D["Ca"], D["Cb"], 1, step)
        if kind == "bo4":
            addb(D["Cb"], A["O4"], step)
            o = rw.AddAtom(Chem.Atom(8)); atom_step[o] = step; atom_unit[o] = d
            addb(D["Ca"], o, step)
            alpha_O[d] = o
        elif kind == "b5":
            addb(D["Cb"], A["C5"], step)
            addb(A["O4"], D["Ca"], step)
        elif kind == "b1":
            if a not in alpha_O:
                raise RuntimeError("beta-1 acceptor %d has no alpha-OH (never donated a beta-O-4)" % a)
            addb(D["Cb"], A["C1"], step)
            for x, y, o in [("C1", "C2", 1), ("C2", "C3", 2), ("C3", "C4", 1), ("C4", "C5", 1),
                            ("C5", "C6", 2), ("C6", "C1", 1)]:
                setb(A[x], A[y], o, step)
            setb(A["C4"], A["O4"], 2, step)
            addb(alpha_O[a], D["Ca"], step)
        elif kind == "bb":
            setb(A["Ca"], A["Cb"], 1, step)
            addb(D["Cb"], A["Cb"], step)
            addb(D["Og"], A["Ca"], step)
            addb(A["Og"], D["Ca"], step)
        elif kind == "55":
            addb(D["C5"], A["C5"], step)
        elif kind == "4o5":
            addb(D["O4"], A["C5"], step)
        else:
            raise RuntimeError("event kind %s not handled (optional channels must be off)" % kind)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    # the untouched side chains stay E (all monolignols enter as E)
    for b in mol.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and not b.IsInRing() \
                and b.GetBeginAtom().GetSymbol() == "C" and b.GetEndAtom().GetSymbol() == "C":
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            ni = [n.GetIdx() for n in b.GetBeginAtom().GetNeighbors() if n.GetIdx() != j]
            nj = [n.GetIdx() for n in b.GetEndAtom().GetNeighbors() if n.GetIdx() != i]
            if ni and nj:
                b.SetStereoAtoms(ni[0], nj[0]); b.SetStereo(Chem.BondStereo.STEREOTRANS)
    return mol, atom_step, atom_unit, bond_hist, bond_kind


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sg", type=float, default=1.0)
    ap.add_argument("--pct-h", type=float, default=0.0)
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--seed", type=int, default=7000)
    ap.add_argument("--params", default=os.path.join("results", "fitted_params_ref_v11e_sw.json"))
    ap.add_argument("--tag", default=None)
    ap.add_argument("--scale-nuclei", action="store_true",
                    help="scale n_nuclei by n/120, so a small molecule keeps the calibrated seed proportion")
    ap.add_argument("--no-3d", action="store_true")
    ap.add_argument("--mmff-iters", type=int, default=3000)
    args = ap.parse_args()

    params = {k: v for k, v in json.load(open(args.params)).items() if not k.startswith("_")}
    params["pct_h"] = args.pct_h
    if args.scale_nuclei:
        params["n_nuclei"] = params["n_nuclei"] * args.n / 120.0
    gen = RecordingKMC(args.sg, n_target=args.n, params=params, seed=args.seed)
    gen.run()
    tag = args.tag or "sg%.1f_n%d_seed%d" % (args.sg, args.n, args.seed)
    os.makedirs(OUT, exist_ok=True)

    t0 = time.time()
    mol, atom_step, atom_unit, bond_hist, bond_kind = build(gen)
    frags = Chem.GetMolFrags(mol)
    rdDepictor.SetPreferCoordGen(False)   # CoordGen returns all-zero coordinates at this size
    rdDepictor.Compute2DCoords(mol)
    c2 = mol.GetConformer()
    xy = [(round(c2.GetAtomPosition(i).x, 3), round(c2.GetAtomPosition(i).y, 3)) for i in range(mol.GetNumAtoms())]
    t2d = time.time() - t0

    xyz, sdf_path, t3d, method = None, None, 0.0, "none"
    if not args.no_3d:
        t0 = time.time()
        # ETKDGv3 on the heavy atoms only: distance geometry costs ~N^3 (triangle smoothing), so
        # leaving the ~1800 hydrogens of a 120-unit lignin out of it cuts the cost ~8x. Hydrogens
        # are then added with coordinates and the whole structure is relaxed with MMFF94.
        method = "ETKDGv3 (heavy atoms) + MMFF94"
        mh0 = Chem.Mol(mol)
        ps = AllChem.ETKDGv3(); ps.randomSeed = args.seed; ps.useRandomCoords = True
        if AllChem.EmbedMolecule(mh0, ps) != 0:
            raise SystemExit("3D embedding failed")
        mh = Chem.AddHs(mh0, addCoords=True)
        mp = AllChem.MMFFGetMoleculeProperties(mh)
        ff = AllChem.MMFFGetMoleculeForceField(mh, mp, nonBondedThresh=9.0)
        ff.Minimize(maxIts=args.mmff_iters)
        heavy = mol.GetNumAtoms()
        conf = mh.GetConformer()
        xyz = [[round(v, 3) for v in conf.GetAtomPosition(i)] for i in range(mh.GetNumAtoms())]
        # hydrogens appear with their heavy atom
        for h in mh.GetAtoms():
            if h.GetIdx() >= heavy:
                p = h.GetNeighbors()[0].GetIdx()
                atom_step[h.GetIdx()] = atom_step[p]; atom_unit[h.GetIdx()] = atom_unit[p]
        mh.SetProp("_Name", "forward generator %s" % tag)
        mh.SetProp("generator", "mini_forward_generator.py v1.1.1")
        mh.SetProp("params", rel(args.params))
        mh.SetProp("geometry", method)
        mh.SetProp("sg", str(args.sg)); mh.SetProp("n_units", str(args.n)); mh.SetProp("seed", str(args.seed))
        sdf_path = os.path.join(OUT, tag + ".sdf")
        w = Chem.SDWriter(sdf_path); w.write(mh); w.close()
        molblock = Chem.MolToMolBlock(mh)
        t3d = time.time() - t0
    else:
        molblock = None

    counts = {}
    for k, _, _ in gen.events:
        counts[k] = counts.get(k, 0) + 1
    timeline = {
        "tag": tag, "sg": args.sg, "pct_h": args.pct_h, "n": args.n, "seed": args.seed,
        "params": rel(args.params), "n_nuclei": params["n_nuclei"], "geometry": method, "smiles": Chem.MolToSmiles(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol), "mw": round(Descriptors.MolWt(mol), 1),
        "fragments": len(frags), "counts": counts,
        "units": [gen.units[u].mtype for u in sorted(gen.units)],
        "unit_step": [gen.unit_step[u] for u in sorted(gen.units)],
        "events": [[k, d, a] for k, d, a in gen.events],
        "atoms": [[mol.GetAtomWithIdx(i).GetSymbol(), xy[i][0], xy[i][1], atom_step[i], atom_unit[i]]
                  for i in range(mol.GetNumAtoms())],
        "bonds": [[i, j, h, bond_kind.get((i, j), "")] for (i, j), h in sorted(bond_hist.items())],
        "atom_step_h": [atom_step[i] for i in range(len(xyz))] if xyz else None,
        "molblock": molblock,
    }
    json_path = os.path.join(OUT, tag + ".json")
    json.dump(timeline, open(json_path, "w"), separators=(",", ":"))
    lines = ["structure_from_run.py  %s" % time.strftime("%Y-%m-%d %H:%M"),
             "params %s  S/G %.2f  H %.2f  N %d  seed %d" % (rel(args.params), args.sg, args.pct_h, args.n, args.seed),
             "events %d  bonds %s  fragments %d" % (len(gen.events), counts, len(frags)),
             "formula %s  MW %.1f  heavy atoms %d" % (timeline["formula"], timeline["mw"], mol.GetNumAtoms()),
             "2D %.1f s  3D %.1f s (%s)" % (t2d, t3d, method),
             "wrote %s%s" % (rel(json_path), ("  " + rel(sdf_path)) if sdf_path else "")]
    open(os.path.join(OUT, tag + ".log"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
