#!/bin/zsh
cd ~/.worktrees/planet-hunt-r6-toi4342
T=runs/round6_toi4342/tri
for cfg in "$@"; do
  set -- ${=cfg}; pl=$1; mode=$2; nruns=$3
  extra=(--contrast-curve $T/zorro832.csv --filt TESS)
  if [[ $mode == ccclear ]]; then
    if [[ $pl == d18 ]]; then extra+=(--clear-beyond 12); else extra+=(--clear-beyond 30); fi
  fi
  nice -n 20 .trivenv/bin/python round6/run_tri.py --tic 354944123 --fold $T/fold_$pl.npz --lc-dir cache/r6/raw --runs $nruns --n 1000000 \
     --max-sectors 9 --out $T --name tri_${pl}_${mode} --local-cutouts cache/tesscut/324.38942_-77.97915 "${extra[@]}" > $T/tri_${pl}_${mode}.log 2>&1
done
