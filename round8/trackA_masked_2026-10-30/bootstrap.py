# Read-only dependency reuse; the specified shared Python remains the interpreter.
import sys, os
from pathlib import Path
sys.path.append('.venv/lib/python3.12/site-packages')
os.environ['PYTHONDONTWRITEBYTECODE']='1'
os.environ['XDG_CACHE_HOME']=str(Path.cwd()/'runtime/cache')
os.environ['XDG_CONFIG_HOME']=str(Path.cwd()/'runtime/config')
os.environ['NUMBA_CACHE_DIR']=str(Path.cwd()/'runtime/numba')
os.environ['MPLCONFIGDIR']=str(Path.cwd()/'runtime/matplotlib')
