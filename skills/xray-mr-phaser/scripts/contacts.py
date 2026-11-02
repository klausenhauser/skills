# usage: phenix.python contacts.py model.pdb   -- CA-CA pairs <8 A between chains incl. symmetry mates
import iotbx.pdb, sys, collections
from scitbx.array_family import flex
pdb = iotbx.pdb.input(sys.argv[1]); h = pdb.construct_hierarchy()
sel = [i for i, a in enumerate(h.atoms()) if a.name == ' CA ']
ch = [h.atoms()[i].parent().parent().parent().id for i in sel]
asu = pdb.xray_structure_simple().select(flex.size_t(sel)).pair_asu_table(distance_cutoff=8.0)
cnt = collections.Counter()
for i, d in enumerate(asu.table()):
    for j, syms in d.items(): cnt[(ch[i], ch[j])] += len(syms)
for c in sorted(set(ch)):
    print(c, {b: cnt[(c, b)] for b in sorted(set(ch)) if b != c and cnt[(c, b)]})
