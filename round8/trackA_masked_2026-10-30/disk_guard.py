"""Disk-only helper extracted from the source job; no network code."""
from pathlib import Path
import shutil
ROOT=Path.cwd()

def space():
 free=shutil.disk_usage(ROOT).free
 size=sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file())
 if free<30_000_000_000 or size>4_800_000_000:raise RuntimeError('Disk stop: '+str((free,size)))
 return free,size
