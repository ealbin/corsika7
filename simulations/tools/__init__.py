
import os
import ROOT as R

from . import file
from . import nuclei
from . import particle
from . import results
from . import settings

def __build_libraries__():
    R.gInterpreter.ProcessLine('.L {}/run.h+'.format(__directory__))
    R.gInterpreter.ProcessLine('.L {}/sim.h+'.format(__directory__))

def reset():
    for f in os.listdir(__directory__):
        if (f.endswith('.d') or f.endswith('.so') or f.endswith('.pcm')):
            os.remove(os.path.join(__directory__, f))
    print()
    print('>> please exit, and restart python (and import tools) to rebuild libraries')

__directory__ = os.path.dirname(os.path.realpath(__file__))
__build_libraries__()
