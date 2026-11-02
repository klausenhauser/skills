#!/bin/zsh
# summarize Phaser runs: ./summ.sh RUN1 RUN2 ...   (each dir has phaser.log)
for r in "$@"; do
  L=$r/phaser.log
  echo "### $r  [finished=$(grep -c 'EXIT STATUS: SUCCESS' $L)] $(grep -oE 'TOP LLG = [0-9.]+|New best solution LLG = [0-9.]+' $L | tail -1)"
  grep -E "SOLU 6DIM" $r/$r.sol 2>/dev/null | awk '{print "  ",$4,$NF}' | tail -12
done
