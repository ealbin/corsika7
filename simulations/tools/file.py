
import ROOT as R

def load(filelist):
    print(' [DISCLAIMER]: assuming files in filelist have the same simulation parameters (energy, observation level, primary particle..)')
    
    run = R.TChain('run')
    sim = R.TChain('sim')

    for f in filelist:
        run.Add(f)
        sim.Add(f)

    return run, sim

