#!/bin/zsh
# usage: round8/run_tri.sh <key> <tic> <config-name> [extra run_tri.py args...]
# Runs TRICERATOPS (round6/run_tri.py) 3 x 1e6 draws on the round-8 fold and apertures.
cd "$(dirname "$0")/.."
key=$1; tic=$2; name=$3; shift 3
T=runs/round8_p1/tri
lc=cache/r8/raw/$tic
[[ -f $T/${key}_apertures.npz ]] && lc=$T/${key}_apertures.npz
cut=$(ls -d cache/tesscut/* | while read d; do ls $d/*.fits 2>/dev/null | head -1 | grep -q . && echo $d; done | head -0)
cut=$(.venv/bin/python -c "import sys,json; sys.path.insert(0,'.'); from planet_hunt import pixelvet as pv; v=json.load(open('runs/round8_p1/inventory.json'))['$key']; print(pv.cut_dir(v['ra'],v['dec']))" 2>/dev/null)
nice -n 15 .trivenv/bin/python round6/run_tri.py --tic $tic --fold $T/${key}_fold.npz --lc-dir $lc --runs 3 --n 1000000 \
  --max-sectors 12 --out $T --name tri_${key}_${name} --local-cutouts $cut "$@" > $T/tri_${key}_${name}.log 2>&1
echo "done $key $name: $(grep -E 'FPP_mean|NFPP_mean' $T/tri_${key}_${name}.json 2>/dev/null | tr -d '\n ')"
