
import os
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


def collect_histograms(filename, filelist, verbose=False):

    def mkdir(name):
        if (R.gDirectory.GetDirectory(name) == None):
            if (verbose):
                print(' [info] making subdir: {}'.format(name))
            R.gDirectory.mkdir(name)
        else:    
            if (verbose):
                print(' [info] subdir exists: {}'.format(name))
        subdir = R.gDirectory.GetDirectory(name)
        R.gDirectory.cd(name)
        return subdir

    search_keys = ['density', 'spectrum', 'efficiency', 'content']

    out = R.TFile(filename, 'recreate')

    nfiles = len(filelist)
    for i, f in enumerate(filelist):
        print(' [info] working on file {} of {}: {}'.format(i+1, nfiles, f))
        tag     = ( '_'.join(os.path.basename(f).split('_')[1:]) ).strip('.root')
        dirname = os.path.dirname(f).lstrip('./')
        subdirs = dirname.split('/')
        out.cd()
        subdir  = mkdir(subdirs[0])
        for _ in subdirs[1:]:
            subdir = mkdir(_)
        data = R.TFile(f)
        subdir.cd()
        if (verbose):
            print(' [info] writing in directory:')
            R.gDirectory.pwd()
        for key in data.GetListOfKeys():
            if any(_ in key.GetName() for _ in search_keys):
                hist = data.Get(key.GetName())
                hist.SetDirectory(subdir)
                hist.Write(tag + key.GetName(), R.TObject.kOverwrite)
                if (verbose):
                    print(' [info]   - wrote {} as {} with title: {}'.format(hist.GetName(), tag + key.GetName(), hist.GetTitle()))
        data.Close()

    out.Close()
    print('finished!')

