
import os
import sys
import importlib
import ROOT as R

from . import file
from . import fit
from . import nuclei
from . import particle
from . import plot
from . import bwplot
from . import results
from . import settings
from . import smartphones

def __build_libraries__():
    R.gInterpreter.ProcessLine('.L {}/run.h+'.format(__directory__))
    R.gInterpreter.ProcessLine('.L {}/sim.h+'.format(__directory__))

def reset():
    for f in os.listdir(__directory__):
        if (f.endswith('.d') or f.endswith('.so') or f.endswith('.pcm')):
            os.remove(os.path.join(__directory__, f))
    print()
    print('>> please exit, and restart python (and import tools) to rebuild libraries')

def reload():
    for k, v in sys.modules.items():
        if (k.startswith('tools')):
            importlib.reload(v)

__directory__ = os.path.dirname(os.path.realpath(__file__))
__build_libraries__()
