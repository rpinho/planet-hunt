#!/bin/zsh
# keep running the hold-out on newly fetched stars while the fetch is alive, then a final pass
cd "$(dirname $0)/.."
while pgrep -f "round8/holdout.py fetch" >/dev/null; do
  nice -n 15 .venv/bin/python round8/holdout.py run --workers 5 --fetched >> runs/round8_crowd/run.log 2>&1
  sleep 90
done
nice -n 15 .venv/bin/python round8/holdout.py run --workers 6 --fetched >> runs/round8_crowd/run.log 2>&1
echo LOOP_DONE2 >> runs/round8_crowd/run.log
