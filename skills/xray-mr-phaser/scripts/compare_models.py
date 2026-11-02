# usage: phenix.python compare_models.py a.pdb chainA b.pdb chainB [offset]
# CA RMSD after superposition on residues common to both (b renumbered by +offset)
import iotbx.pdb, sys
from scitbx.array_family import flex
from scitbx.math import superpose
def ca(f, ch, off=0):
    d = {}
    for c in iotbx.pdb.input(f).construct_hierarchy().models()[0].chains():
        if ch in ('*', c.id):
            for r in c.residues():
                for at in r.atoms():
                    if at.name == ' CA ': d[r.resseq_as_int() + off] = at.xyz
    return d
a = ca(sys.argv[1], sys.argv[2]); b = ca(sys.argv[3], sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 0)
k = sorted(set(a) & set(b)); A = flex.vec3_double([a[i] for i in k]); B = flex.vec3_double([b[i] for i in k])
s = superpose.least_squares_fit(A, B); print(len(k), 'CA, rmsd', A.rms_difference(s.other_sites_best_fit()))
