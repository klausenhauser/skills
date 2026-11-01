"""
Split an intact fitted ligand into its component residues, preserving coordinates.

REFERENCE IMPLEMENTATION -- adapt the CONFIGURE block below for each project.
Everything under "matching logic" is general and should not need changing.

The idea: rather than rebuilding a multi-residue ligand from monomers and
trying to place them (which gives junctions ~1 A too long that no program will
link), carve the residues out of the ligand you have already fitted. Atoms do
not move, so every junction comes out at a real bond length by construction.

How it works:
  1. build each component's heavy-atom graph from its restraint CIF, MINUS its
     leaving atoms, so the graph is the residue as it exists inside the ligand
  2. perceive connectivity on the fitted ligand from coordinates alone
     (DetermineConnectivity) -- no bond orders, so no aromatic-perception traps
  3. match with bond orders stripped from both sides: element + topology only
  4. backtrack to an EXACT COVER -- placements mutually disjoint and together
     covering every heavy atom

Step 4 is the safeguard. An exact cover with unique per-component matches is a
proof the assignment is right, not a guess. No cover found means the component
list or the leaving atoms are wrong, and you get an error instead of a
plausible-looking wrong answer.

Always check the junction distances afterwards.
"""
import gemmi
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from itertools import product

# ---------------------------------------------------------------- CONFIGURE
# Everything in this block is project-specific. Change all of it.
#
#   MODEL       the model with the full ligand fitted as one residue
#   LIGAND      (chain, residue name) of that fitted ligand -- the copy to split
#   OUT_CHAIN   chain name for the split residues in the output
#   CIFDIR      where the post-processed component CIFs live
#   COMPONENTS  (id, leaving atoms to drop) per component that has a CIF;
#               standard CCD components with no Grade2 output are built by
#               hand with manual() further down
#   NUM         component id -> residue number in the output
#   JUNCTIONS   inter-residue bonds to check afterwards, with ideal lengths
#   ANCHOR      the covalent bond to the protein to check, or None
#
# The values below are placeholders showing the format: a four-residue
# peptide RES1-RES4 with a norleucine and C-terminal amide cap on chain C,
# where RES3 carries the warhead (leaving atom F1) bonded to TYR 123 OH on
# protein chain A.
#
MODEL = "/path/to/model_with_fitted_ligand.cif"
LIGAND = ("C", "LIG")
OUT_CHAIN = "C"
CIFDIR = "/path/to/Step 2"

def comp_from_cif(path, comp_id, drop=()):
    """Build an RDKit mol (heavy atoms) from a CCP4 restraint CIF, carrying atom names."""
    doc = gemmi.cif.read_file(path)
    names, elems, bonds = [], [], []
    for blk in doc:
        t = blk.find("_chem_comp_atom.", ["atom_id", "type_symbol"])
        for row in t:
            n = row.str(0)
            e = row.str(1)
            if e == "H" or n in drop:
                continue
            names.append(n); elems.append(e)
        t = blk.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2"])
        for row in t:
            a, b = row.str(0), row.str(1)
            if a in names and b in names:
                bonds.append((names.index(a), names.index(b)))
    m = Chem.RWMol()
    for e in elems:
        m.AddAtom(Chem.Atom(e.capitalize()))  # Grade2 writes CL/BR in capitals; RDKit needs Cl/Br
    for a, b in bonds:
        m.AddBond(a, b, Chem.BondType.SINGLE)
    return m.GetMol(), names

def manual(atoms, bonds):
    m = Chem.RWMol()
    for n, e in atoms:
        m.AddAtom(Chem.Atom(e.capitalize()))  # Grade2 writes CL/BR in capitals; RDKit needs Cl/Br
    names = [n for n, _ in atoms]
    for a, b in bonds:
        m.AddBond(names.index(a), names.index(b), Chem.BondType.SINGLE)
    return m.GetMol(), names

# ---- the components, as they appear INSIDE the peptide (leaving atoms dropped)
COMPONENTS = []
for cid, drop in [("RES1", ("OXT",)), ("RES2", ("OXT",)),
                  ("RES3", ("OXT", "F1")), ("RES4", ("OXT",))]:
    mol, names = comp_from_cif(f"{CIFDIR}/{cid}.restraints.fixed.cif", cid, drop)
    COMPONENTS.append((cid, mol, names))
COMPONENTS.append(("NLE", *manual(
    [("N","N"),("CA","C"),("C","C"),("O","O"),("CB","C"),("CG","C"),("CD","C"),("CE","C")],
    [("N","CA"),("CA","C"),("C","O"),("CA","CB"),("CB","CG"),("CG","CD"),("CD","CE")])))
COMPONENTS.append(("NH2", *manual([("N","N")], [])))

NUM = {"RES1": 1, "RES2": 2, "RES3": 3, "RES4": 4, "NLE": 5, "NH2": 6}

# (residue 1, atom 1, residue 2, atom 2, ideal length in A)
JUNCTIONS = [(1,'C',2,'N',1.34), (2,'C',3,'N',1.33), (3,'C',4,'N',1.33),
             (4,'C',5,'N',1.33), (5,'C',6,'N',1.33)]

# (protein chain, resname, resnum, atom, ligand residue number, ligand atom, ideal length in A)
ANCHOR = ("A", "TYR", 123, "OH", 3, "S1", 1.60)

# --------------------------------------------------------- end CONFIGURE
# ---- the fitted ligand, connectivity perceived from coordinates
st = gemmi.read_structure(MODEL)
lig_chain, lig_name = LIGAND
lig = [r for ch in st[0] if ch.name == lig_chain for r in ch if r.name == lig_name]
if not lig:
    raise SystemExit(f"no residue {lig_name} in chain {lig_chain} of {MODEL} -- check LIGAND")
lig = lig[0]
heavy = [a for a in lig if a.element.name != "H"]
print(f"{lig_name} (chain {lig_chain}) heavy atoms: {len(heavy)}")

xyz = f"{len(heavy)}\n\n" + "".join(
    f"{a.element.name} {a.pos.x:.4f} {a.pos.y:.4f} {a.pos.z:.4f}\n" for a in heavy)
raw = Chem.MolFromXYZBlock(xyz)
rdDetermineBonds.DetermineConnectivity(raw, charge=0)
print(f"perceived bonds: {raw.GetNumBonds()}")

# ================= matching logic: general, leave unchanged =================
# ---- match each component; bond orders ignored, elements enforced
def generic(q):
    qq = Chem.RWMol(q)
    for b in qq.GetBonds():
        b.SetBondType(Chem.BondType.UNSPECIFIED)
    p = Chem.AdjustQueryParameters.NoAdjustments()
    p.makeBondsGeneric = True
    return Chem.AdjustQueryProperties(qq.GetMol(), p)

target = generic(raw)
cands = {}
for cid, mol, names in COMPONENTS:
    ms = target.GetSubstructMatches(generic(mol), uniquify=True, maxMatches=5000)
    cands[cid] = ms
    print(f"  {cid:6s} {len(names):2d} atoms -> {len(ms)} candidate placements")

# ---- exact cover: disjoint placements covering every heavy atom once
order = sorted(cands, key=lambda c: len(cands[c]))
solution = {}
def solve(i, used):
    if i == len(order):
        return len(used) == len(heavy)
    cid = order[i]
    for m in cands[cid]:
        s = set(m)
        if len(s) == len(m) and not (s & used):
            solution[cid] = m
            if solve(i + 1, used | s):
                return True
            del solution[cid]
    return False

if not solve(0, set()):
    raise SystemExit("no disjoint full-coverage partition found")
print("\nexact cover found: every heavy atom assigned exactly once")

# ---- emit the assignment
assign = {}
for cid, mol, names in COMPONENTS:
    for qi, ti in enumerate(solution[cid]):
        assign[ti] = (NUM[cid], cid, names[qi])
import json
json.dump({str(k): v for k, v in assign.items()}, open("assignment.json", "w"), indent=0)
for cid in sorted(NUM, key=NUM.get):
    got = [v[2] for v in assign.values() if v[1] == cid]
    print(f"  {NUM[cid]} {cid:6s}: {' '.join(sorted(got))}")

# ---- build the split ligand, coordinates untouched
out = gemmi.Structure()
out.spacegroup_hm = st.spacegroup_hm
out.cell = st.cell
mdl = gemmi.Model("1")
ch = gemmi.Chain(OUT_CHAIN)
byres = {}
for i, a in enumerate(heavy):
    num, cid, nm = assign[i]
    byres.setdefault((num, cid), []).append((nm, a))
for (num, cid) in sorted(byres):
    r = gemmi.Residue()
    r.name = cid
    r.seqid = gemmi.SeqId(num, ' ')
    r.het_flag = 'H'
    for nm, a in byres[(num, cid)]:
        na = gemmi.Atom()
        na.name = nm
        na.element = a.element
        na.pos = a.pos
        na.occ = a.occ
        na.b_iso = a.b_iso
        r.add_atom(na)
    ch.add_residue(r)
mdl.add_chain(ch)
out.add_model(mdl)
out.setup_entities()
out.make_mmcif_document().write_file("ligand_split.cif")

def A(num, nm):
    for r in out[0][OUT_CHAIN]:
        if r.seqid.num == num:
            for a in r:
                if a.name == nm: return a
print("\n=== junction distances in the split model ===")
for n1, a1, n2, a2, ideal in JUNCTIONS:
    d = A(n1,a1).pos.dist(A(n2,a2).pos)
    flag = "ok" if abs(d-ideal) < 0.15 else "CHECK"
    print(f"  {n1}{a1} -> {n2}{a2}: {d:.3f} A  (ideal {ideal})  {flag}")
if ANCHOR:
    pchain, pname, pnum, patom, lnum, latom, ideal = ANCHOR
    pres = [r for c in st[0] if c.name==pchain for r in c if r.name==pname and r.seqid.num==pnum][0]
    pa = [a for a in pres if a.name==patom][0]
    print(f"  {pname}{pnum}({pchain}) {patom} -> {lnum}{latom}: {pa.pos.dist(A(lnum,latom).pos):.3f} A  (ideal {ideal})")
