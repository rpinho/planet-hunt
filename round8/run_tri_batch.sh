#!/bin/zsh
# TRICERATOPS (round6/run_tri.py) on round-8 candidates; fold/apertures from round8/trifold.py
cd "$(dirname $0)/.."
T=runs/round8_crowd/tri
run() {  # id tic runs n
  local cut=$(ls -d cache/tesscut/$(python3 -c "import json;print([f\"{s['ra']:.5f}\" for s in json.load(open('runs/round8_crowd/deepvet_targets.json'))['targets'] if s['name']=='$1'][0])")_* | head -1)
  nice -n 15 .trivenv/bin/python round6/run_tri.py --tic $2 --fold $T/$1_fold.npz --lc-dir $T/$1_apertures.npz --runs $3 --n $4 \
     --out $T --name triceratops_$1 --local-cutouts $cut > $T/tri_$1.log 2>&1
}
run E1_101653309_0.7034 101653309 5 1000000 &
run E1_262843927_24.1087 262843927 3 1000000 &
wait
run E1_410418450_6.8202 410418450 3 1000000 &
run E1_317300777_2.8211 317300777 3 500000 &
wait
run E1_348760393_12.3887 348760393 3 500000
echo TRI_BATCH_DONE >> $T/tri_batch.done
