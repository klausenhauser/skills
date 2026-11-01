#!/usr/bin/env python3
"""
make_deposition_cif.py  <refined.cif>  [out.cif]

Prepare a Phenix-refined model for OneDep:
  1. drop the comp_*/link_*/mod_* restraint blocks Phenix embeds in its output
  2. record, in _struct_conn, the covalent bonds that standard polymer
     connectivity does not imply -- the protein-ligand warhead bond and the
     inter-residue bonds of a split multi-residue ligand

PHENIX DOES WRITE _struct_conn. It finds bonds on its own -- notably
crystallographic-symmetry disulfides -- and writes them into its output. So
this script READS the existing loop, KEEPS every row, and APPENDS only what is
missing. It never rebuilds the loop.

That distinction matters: an earlier version of this script used
block.init_loop(), which REPLACES the loop, and silently dropped a
symmetry disulfide (a cysteine to its own symmetry mate) from a finished
deposition file. The row was only noticed later, by accident.

The reason that version was written at all is worth knowing, because the same
mistake is easy to repeat: the loop was queried with a tag list that included
_struct_conn.pdbx_dist_value, which Phenix does NOT write. Asking gemmi's
find() for several tags at once returns an EMPTY table if any single one is
missing -- which reads as "there are no rows" and is not. Hence read_existing()
below probes one tag at a time and tolerates Phenix's own tag set.

CAVEAT worth passing to the annotator: Phenix writes no _struct_conn.ptnr2_symmetry,
so a symmetry bond arrives looking like a chemically impossible self-bond
(same chain, same residue, same atom). This script flags such rows. Either add
the operator in wwPDB notation or mention it in the deposition notes.

EDIT THE LINKS TABLE BELOW for your structure. Distances are measured from the
coordinates, never hardcoded.
"""
import sys
import gemmi

# (id, chain_1, resnum_1, comp_1, atom_1, chain_2, resnum_2, comp_2, atom_2,
#  leaving_atom_flag, details)
# Placeholder example: a four-residue peptide RES1-RES4 plus NLE and an amide
# cap on chain C, with RES3's sulfonyl fluoride warhead bonded to TYR 123 OH.
LINKS = [
    ("covale1", "A", 123, "TYR",  "OH", "C", 3, "RES3", "S1", "both",
     "SuFEx: tyrosine O-sulfonylation; Tyr HH and sulfonyl F eliminated"),
    ("covale2", "C", 1,   "RES1", "C",  "C", 2, "RES2", "N",  "both", "N-terminal acyl cap"),
    ("covale3", "C", 2,   "RES2", "C",  "C", 3, "RES3", "N",  "both", "peptide bond"),
    ("covale4", "C", 3,   "RES3", "C",  "C", 4, "RES4", "N",  "both", "peptide bond"),
    ("covale5", "C", 4,   "RES4", "C",  "C", 5, "NLE",  "N",  "both", "peptide bond"),
    ("covale6", "C", 5,   "NLE",  "C",  "C", 6, "NH2",  "N",  "both", "C-terminal amide cap"),
]

# Tags this script writes. Phenix's own rows may carry fewer; missing values
# are written as "?" so the loop stays rectangular.
TAGS = ["id", "conn_type_id",
        "ptnr1_auth_asym_id", "ptnr1_auth_seq_id", "ptnr1_auth_comp_id",
        "ptnr1_label_atom_id", "ptnr1_symmetry",
        "ptnr2_auth_asym_id", "ptnr2_auth_seq_id", "ptnr2_auth_comp_id",
        "ptnr2_label_atom_id", "ptnr2_symmetry",
        "pdbx_dist_value", "pdbx_leaving_atom_flag", "details"]


def model_block(doc):
    for block in doc:
        if block.find_values("_atom_site.Cartn_x"):
            return block
    sys.exit("no block with coordinates found")


def read_existing(block):
    """Existing _struct_conn rows as dicts keyed by tag. Tolerates Phenix's
    tag set by probing one tag at a time -- see the note in the docstring."""
    cols, n = {}, 0
    for tag in set(TAGS) | {"pdbx_role", "pdbx_PDB_id"}:
        vals = block.find_values("_struct_conn." + tag)
        if vals:
            cols[tag] = [v.strip() for v in vals]
            n = max(n, len(vals))
    if not n:
        return []
    return [{t: cols[t][i] for t in cols} for i in range(n)]


def identity(row):
    """What makes a row the same bond, for duplicate detection."""
    def side(p):
        return (row.get(f"ptnr{p}_auth_asym_id", "?"),
                row.get(f"ptnr{p}_auth_seq_id", "?"),
                row.get(f"ptnr{p}_label_atom_id", "?"))
    return tuple(sorted([side(1), side(2)]))


def main(src, dst=None):
    dst = dst or src.replace(".cif", "_deposit.cif")

    st = gemmi.read_structure(src)
    st.setup_entities()
    model = st[0]
    n_before = sum(len(r) for c in model for r in c)

    def atom(ch, num, comp, name):
        for c in model:
            if c.name != ch:
                continue
            for r in c:
                if r.seqid.num == num and r.name == comp:
                    for a in r:
                        if a.name == name:
                            return a
        sys.exit(f"not found: {ch}/{num} {comp} {name}")

    doc = gemmi.cif.read_file(src)
    removed = []
    for i in range(len(doc) - 1, -1, -1):
        if doc[i].name.startswith(("comp_", "link_", "mod_")):
            removed.append(doc[i].name)
            del doc[i]
    if removed:
        print(f"removed {len(removed)} restraint block(s)")
    block = model_block(doc)

    existing = read_existing(block)
    print(f"existing _struct_conn rows: {len(existing)}")

    rows, seen = [], set()
    for row in existing:
        rows.append([row.get(t, "?") or "?" for t in TAGS])
        seen.add(identity(row))
        label = (f"{row.get('ptnr1_auth_comp_id','?')} "
                 f"{row.get('ptnr1_auth_asym_id','?')}"
                 f"{row.get('ptnr1_auth_seq_id','?')} "
                 f"{row.get('ptnr1_label_atom_id','?')} -- "
                 f"{row.get('ptnr2_auth_comp_id','?')} "
                 f"{row.get('ptnr2_auth_asym_id','?')}"
                 f"{row.get('ptnr2_auth_seq_id','?')} "
                 f"{row.get('ptnr2_label_atom_id','?')}")
        print(f"  kept {row.get('id','?')} ({row.get('conn_type_id','?')}): {label}")
        s1, s2 = (row.get("ptnr1_symmetry", "?"), row.get("ptnr2_symmetry", "?"))
        if identity(row)[0] == identity(row)[1] and s2 in ("?", ".", "", "1_555"):
            print("       ** self-bond with no symmetry operator. This is a "
                  "crystallographic-symmetry bond that Phenix recorded without "
                  "ptnr2_symmetry. Add the operator or flag it in the notes. **")

    added = 0
    for cid, c1, n1, k1, a1, c2, n2, k2, a2, leave, note in LINKS:
        cand = {"ptnr1_auth_asym_id": c1, "ptnr1_auth_seq_id": str(n1),
                "ptnr1_label_atom_id": a1,
                "ptnr2_auth_asym_id": c2, "ptnr2_auth_seq_id": str(n2),
                "ptnr2_label_atom_id": a2}
        if identity(cand) in seen:
            print(f"  skipped {cid}: already recorded")
            continue
        d = atom(c1, n1, k1, a1).pos.dist(atom(c2, n2, k2, a2).pos)
        vals = {"id": cid, "conn_type_id": "covale",
                "ptnr1_auth_asym_id": c1, "ptnr1_auth_seq_id": str(n1),
                "ptnr1_auth_comp_id": k1, "ptnr1_label_atom_id": a1,
                "ptnr1_symmetry": "1_555",
                "ptnr2_auth_asym_id": c2, "ptnr2_auth_seq_id": str(n2),
                "ptnr2_auth_comp_id": k2, "ptnr2_label_atom_id": a2,
                "ptnr2_symmetry": "1_555",
                "pdbx_dist_value": f"{d:.2f}",
                "pdbx_leaving_atom_flag": leave,
                "details": gemmi.cif.quote(note)}
        rows.append([vals[t] for t in TAGS])
        seen.add(identity(cand))
        added += 1
        print(f"  added {cid}: {c1}/{n1} {k1} {a1} -- {c2}/{n2} {k2} {a2}   {d:.2f} A")

    loop = block.init_loop("_struct_conn.", TAGS)
    for r in rows:
        loop.add_row(r)

    types = sorted({r[TAGS.index("conn_type_id")] for r in rows})
    ct = block.init_loop("_struct_conn_type.", ["id", "criteria", "reference"])
    for t in types:
        ct.add_row([t, "?", "?"])

    doc.write_file(dst)

    # verify
    st2 = gemmi.read_structure(dst)
    st2.setup_entities()
    n_after = sum(len(r) for c in st2[0] for r in c)
    pa = [(a.pos.x, a.pos.y, a.pos.z) for c in model for r in c for a in r]
    pb = [(a.pos.x, a.pos.y, a.pos.z) for c in st2[0] for r in c for a in r]
    dev = max((max(abs(x - y) for x, y in zip(p, q)) for p, q in zip(pa, pb)),
              default=0.0)
    print(f"\nwrote {dst}")
    print(f"  _struct_conn: {len(existing)} preserved + {added} added "
          f"= {len(rows)} rows")
    print(f"  atoms: {n_before} -> {n_after}   max coordinate deviation: {dev:.6f} A")
    if n_before != n_after or dev > 1e-4:
        sys.exit("  ** model changed -- do not deposit this file **")


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    main(*sys.argv[1:])
