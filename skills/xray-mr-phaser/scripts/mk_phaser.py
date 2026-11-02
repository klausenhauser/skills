#!/usr/bin/env python3
"""Generate a Phaser MR_AUTO keyword file.
usage: mk_phaser.py RUN MTZ F SIGF HIRES "SPACEGROUP" --seq a.fa:N --seq b.fa:N \
         --ens name=path.pdb:RMS [--ens ...] [--sgalt] SEARCH [SEARCH ...]
Each SEARCH is one placement step; use a|b for "a OR b". Writes RUN/phaser.in.
Run: (cd RUN && nohup phenix.phaser < phaser.in > phaser.log 2>&1 &)"""
import sys, os, argparse
p = argparse.ArgumentParser()
p.add_argument('run'); p.add_argument('mtz'); p.add_argument('f'); p.add_argument('sig')
p.add_argument('hires'); p.add_argument('sg')
p.add_argument('--seq', action='append', default=[], help='fasta:copies')
p.add_argument('--ens', action='append', default=[], help='name=pdb:rms')
p.add_argument('--sgalt', action='store_true'); p.add_argument('--kill', default='150')
p.add_argument('search', nargs='+')
a = p.parse_args()
os.makedirs(a.run, exist_ok=True)
ens = {}
for e in a.ens:
    n, rest = e.split('=', 1); path, rms = rest.rsplit(':', 1); ens[n] = (os.path.abspath(path), rms)
t = f"TITLE {a.run}\nMODE MR_AUTO\nHKLIN {os.path.abspath(a.mtz)}\nLABIN F={a.f} SIGF={a.sig}\nSPACEGROUP {a.sg}\nRESOLUTION {a.hires}\nCOMPOSITION BY ASU\n"
for s in a.seq:
    f, n = s.rsplit(':', 1); t += f"COMPOSITION PROTEIN SEQUENCE {os.path.abspath(f)} NUMBER {n}\n"
used = {w for s in a.search for w in s.replace('|', ' ').split()}
for k in sorted(used):
    path, rms = ens[k]; t += f"ENSEMBLE {k} PDB {path} RMS {rms}\n"
for s in a.search:
    t += 'SEARCH ' + ' OR '.join('ENSEMBLE ' + x for x in s.split('|')) + ' NUM 1\n'
t += "PACK SELECT PERCENT\nPACK CUTOFF 10\nTNCS RLIST ADD ON\n"
if a.sgalt: t += "SGALT SELECT ALL\n"
t += f"JOBS 1\nKILL TIME {a.kill}\nROOT {a.run}\nXYZOUT ON\nHKLOUT OFF\n"
open(f'{a.run}/phaser.in', 'w').write(t); print('wrote', f'{a.run}/phaser.in')
