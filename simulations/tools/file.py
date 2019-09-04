
import ROOT as R

from . import results

def load(filelist):
    print(' [DISCLAIMER]: assuming files in filelist have the same simulation parameters (energy, observation level, primary particle..)')
    
    run = R.TChain('run')
    sim = R.TChain('sim')

    for f in filelist:
        run.Add(f)
        sim.Add(f)

    return run, sim

def batch(file):
    R.gROOT.SetBatch(R.kTRUE)

    f = R.TFile(file, 'update')
    r = results.Results()

    nlevels = r.run.run_ObservationLevel.size()
    for _ in range(-1, nlevels):
        print('Working on observation index = {} of {}'.format(_, nlevels))
        r.plot_density(obs_level=_,  save=True)
        r.plot_slice(obs_level=_,    save=True)
        r.plot_spectrum(obs_level=_, save=True)
        r.plot_content(obs_level=_,  save=True)
    print('Working on impact')
    r.plot_impact(save=True)
    print('Working on efficiency')
    r.plot_efficiency(save=True)
    
    print('Finished')
    f.Close()

