
import os
import numpy as np
import ROOT as R

from . import nuclei
from . import particle
from . import settings

__directory__ = os.path.dirname(os.path.realpath(__file__))
R.gInterpreter.ProcessLine('.L {}/run.h+'.format(__directory__))
R.gInterpreter.ProcessLine('.L {}/sim.h+'.format(__directory__))

__data__ = R.TFile('./proton-1e14/DAT000001.root')


class Results:
    
    def __init__(self):
        filename = R.gDirectory.GetName()
        is_data_file  = filename.endswith('.root')
        has_run_tree  = False
        has_sim_tree  = False

        if (is_data_file):
            n_keys = R.gDirectory.GetNkeys()
            keys   = R.gDirectory.GetListOfKeys()
            for _ in range(n_keys):
                key_name = keys.At(_).GetName()

                if (key_name == 'run'):
                    self.run = R.Run(R.gDirectory.Get('run'))
                    self.run.GetEntry(0)
                    has_run_tree = True
                    print(" [info] found run tree")
                    continue

                if (key_name == 'sim'):
                    self.sim = R.Sim(R.gDirectory.Get('sim'))
                    self.sim.GetEntry(0)
                    has_sim_tree = True
                    print(" [info] found sim tree")
                    continue

        self.all_good = has_run_tree and has_sim_tree
        if (self.all_good):
            print(" [info] saul goodman")
            self.__initialize_histograms__();
        else:
            print(" [ERROR] snot goodman")
            if (not is_data_file):
                print("         - no root file!")
            else:
                if (not has_run_tree):
                    print("         - no run tree!")
                if (not has_sim_tree):
                    print("         - no sim tree!")


    ##########################################################################

    def observation_levels(self):
        if (not self.all_good):
            self.__print_std_err__()
            return

        print()
        print('index     altitude [km]')
        print('-----------------------')

        nlevels = self.run.run_ObservationLevel.size()
        for _ in range(nlevels):
            index = _ + 1
            if (_ == 9):
                index = 0
            altitude = self.run.run_ObservationLevel.at(_) / 1e5
            print('{: >3d}{}{: >7.2f}'.format(index, ' '*7, altitude))

    ##########################################################################

    def plot_density(self, obs_level=-1):
        if (not self.all_good):
            self.__print_std_err__()
            return

        maxlevel = self.run.run_ObservationLevel.size() - 1
        if (obs_level < -1 or obs_level > maxlevel):
            print(' [ERROR] Invalid observation level index')
            self.observation_levels()
            return

        title = 'Composit of All Observation Levels'
        if (obs_level != -1):
            index = obs_level - 1
            if (obs_level == 0):
                index = 9
            title = 'Observation Level = {:.2f} [km]'.format(
                    self.run.run_ObservationLevel[index]/1e5)

        for _ in settings.Catagories:
            self.density[_].SetTitle(title)
            self.density[_].Reset('ICEM')

        for _ in range(self.sim.particle__):
            
            if (obs_level != -1):
                if (self.sim.particle__ObservationLevel[_] != obs_level):
                    continue

            pid = self.sim.particle__ParticleID[_]
            x   = self.sim.particle__x[_] / 1e5 # [km]
            y   = self.sim.particle__y[_] / 1e5 # [km]

            if (str(pid) in particle.ID):
                catagory = particle.ID[str(pid)]['catagory']
            elif (str(pid)[-2:] in nuclei.ID):
                catagory = nuclei.ID[str(pid)[-2:]]['catagory']
            else:
                print(' [ERROR] Unknown particle ID: {}'.format(pid))
            
            r = np.sqrt(x*x + y*y)
            self.density[catagory].Fill(r)

        xaxis = self.density[settings.Catagories[0]].GetXaxis()
        for _ in range(xaxis.GetNbins()):
            i     = _ + 1
            low   = xaxis.GetBinLowEdge(i)
            width = xaxis.GetBinWidth(i)
            high  = low + width
            area  = np.pi * (high*high - low*low) # [km^2]
            area  = area * (1e5)*(1e5) # [cm^2]

            for cat in settings.Catagories:
                raw = self.density[cat].GetBinContent(_)
                self.density[cat].SetBinContent(_, raw / area)

                raw = self.density[cat].GetBinError(_)
                self.density[cat].SetBinError(_, raw / area)

        try:
            self.density_canvas
        except:
            pass
        else:
            del self.density_canvas

        name   = 'density_canvas'
        width  = settings.Density['width']
        height = settings.Density['height']
        left   = settings.Margins['left']
        right  = settings.Margins['right']
        bottom = settings.Margins['bottom']
        top    = settings.Margins['top']
        self.density_canvas = R.TCanvas(name, title, width, height)
        self.density_canvas.SetFillStyle(4000)
        self.density_canvas.SetMargin(left, right, bottom, top)
        self.density_canvas.SetLogx()
        self.density_canvas.SetLogy()

        for _ in settings.Catagories:
            self.density[_].SetMinimum(np.power(10., settings.Density['ymin']))
            self.density[_].SetMaximum(np.power(10., settings.Density['ymax']))
            self.density[_].Draw('hist pe same')

        legend = R.TLegend(.6, .89, .94, .7)
        legend.SetFillStyle(4000)
        legend.SetBorderSize(0)
        legend.SetNColumns(2)
        legend.AddEntry(self.density['photon'],   '#gamma',      'pe')
        legend.AddEntry(self.density['proton'],   'p, #bar{p}',  'pe')
        legend.AddEntry(self.density['electron'], 'e^{#pm}',     'pe')
        legend.AddEntry(self.density['neutron'],  'n, #bar{n}',  'pe')
        legend.AddEntry(self.density['muon'],     '#mu^{#pm}',   'pe')
        legend.AddEntry(self.density['ocharged'], 'other^{#pm}', 'pe')
        legend.AddEntry(self.density['nuclei'],   'nuclei',      'pe')
        legend.AddEntry(self.density['oneutral'], 'other^{0}',   'pe')
        legend.Draw()
        self.density_canvas.Update()

    ##########################################################################

    def plot_slice(self, obs_level=-1):
        if (not self.all_good):
            self.__print_std_err__()
            return

        maxlevel = self.run.run_ObservationLevel.size() - 1
        if (obs_level < -1 or obs_level > maxlevel):
            print(' [ERROR] Invalid observation level index')
            self.observation_levels()
            return

        title = 'Composit of All Observation Levels'
        if (obs_level != -1):
            index = obs_level - 1
            if (obs_level == 0):
                index = 9
            title = 'Observation Level = {:.2f} [km]'.format(
                    self.run.run_ObservationLevel[index]/1e5)

        for _ in settings.Catagories:
            self.slice[_].SetTitle(title)
            self.slice[_].Reset('ICEM')

        for _ in range(self.sim.particle__):
            
            if (obs_level != -1):
                if (self.sim.particle__ObservationLevel[_] != obs_level):
                    continue

            pid = self.sim.particle__ParticleID[_]
            x   = self.sim.particle__x[_] / 1e5 # [km] North
            y   = self.sim.particle__y[_] / 1e5 # [km] West

            if (str(pid) in particle.ID):
                catagory = particle.ID[str(pid)]['catagory']
            elif (str(pid)[-2:] in nuclei.ID):
                catagory = nuclei.ID[str(pid)[-2:]]['catagory']
            else:
                print(' [ERROR] Unknown particle ID: {}'.format(pid))
            
            # "x" is West-to-East (i.e. -y)
            # "y" is South-to-North (i.e. x)
            self.slice[catagory].Fill(-y, x)

        try:
            self.slice_canvas
        except:
            pass
        else:
            del self.slice_canvas

        name   = 'slice_canvas'
        width  = settings.Slice['width']
        height = settings.Slice['height']
        left   = settings.Margins['left']
        right  = .17 #settings.Margins['right']
        bottom = settings.Margins['bottom']
        top    = settings.Margins['top']
        self.slice_canvases = {}
        for _ in settings.Catagories:
            
            # skip empty plots
            if (self.slice[_].GetEntries() == 0.):
                continue

            canvas = R.TCanvas(name + '_' + _, title, width, height)
            canvas.SetFillStyle(4000)
            canvas.SetMargin(left, right, bottom, top)
            canvas.SetLogz()
            self.slice_canvases[_] = canvas

            self.slice[_].SetMinimum(np.power(10., settings.Slice['zmin']))
            self.slice[_].SetMaximum(np.power(10., settings.Slice['zmax']))
            self.slice[_].Draw('colz')
            
            legend = R.TLegend(.6, .89, .8, .8)
            legend.SetFillStyle(4000)
            legend.SetBorderSize(0)
            legend.AddEntry(self.slice[_], settings.Symbols[_], '')
            legend.Draw()
            canvas.Update()

    ##########################################################################

    def __log_binning__(self, axis):
        nbins = axis.GetNbins()
        axmin = axis.GetXmin()
        axmax = axis.GetXmax()
        width = (axmax - axmin) / float(nbins)

        logbins = []
        for _ in range(nbins + 1):
            logbins.append( np.power(10., axmin + _*width) )
        axis.Set(nbins, np.asarray(logbins, dtype=np.float))
   
    ##########################################################################

    def __initialize_histograms__(self):
        
        # Longitudinal Density
        #---------------------------------------------------------------------
        name   = settings.Density['name']
        title  = settings.Density['title']
        xtitle = settings.Density['xtitle']
        ytitle = settings.Density['ytitle']
        title  = '{};{};{}'.format(title, xtitle, ytitle)
        stats  = settings.Density['stats']
        xbins  = settings.Density['xbins']
        xmin   = settings.Density['xmin']
        xmax   = settings.Density['xmax']

        self.density = {}
        for _ in settings.Catagories:
            self.density[_] = R.TH1D(name + _, title, xbins, xmin, xmax)
            self.density[_].SetStats(stats)
            self.density[_].SetMarkerColor(settings.Colors[_])
            self.density[_].SetMarkerStyle(R.kFullCircle)
            self.density[_].SetLineColor(settings.Colors[_])
            self.density[_].GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
            self.density[_].GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
            self.__log_binning__(self.density[_].GetXaxis())


        # X-Y Slice
        #---------------------------------------------------------------------
        name   = settings.Slice['name']
        title  = settings.Slice['title']
        xtitle = settings.Slice['xtitle']
        ytitle = settings.Slice['ytitle']
        ztitle = settings.Slice['ztitle']
        title  = '{};{};{};{}'.format(title, xtitle, ytitle, ztitle)
        stats  = settings.Slice['stats']
        xbins  = settings.Slice['xbins']
        ybins  = settings.Slice['ybins']
        xmin   = settings.Slice['xmin']
        ymin   = settings.Slice['ymin']
        xmax   = settings.Slice['xmax']
        ymax   = settings.Slice['ymax']

        self.slice = {}
        for _ in settings.Catagories:
            self.slice[_] = R.TH2D(name + _, _ + title, 
                                   xbins, xmin, xmax, ybins, ymin, ymax)
            self.slice[_].SetStats(stats)
            self.slice[_].GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
            self.slice[_].GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
            self.slice[_].GetZaxis().SetTitleOffset(settings.TitleOffset['z'])


        # Energy Spectrum
        #---------------------------------------------------------------------
        name   = settings.Spectrum['name']
        title  = settings.Spectrum['title']
        xtitle = settings.Spectrum['xtitle']
        ytitle = settings.Spectrum['ytitle']
        title  = '{};{};{}'.format(title, xtitle, ytitle)
        stats  = settings.Spectrum['stats']
        xbins  = settings.Spectrum['xbins']
        xmin   = settings.Spectrum['xmin']
        xmax   = settings.Spectrum['xmax']

        self.spectrum = {}
        for _ in settings.Catagories:
            self.spectrum[_] = R.TH1D(name + _, title, xbins, xmin, xmax)
            self.spectrum[_].SetStats(stats)
            self.spectrum[_].SetLineColor(settings.Colors[_])
            self.spectrum[_].GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
            self.spectrum[_].GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
            self.__log_binning__(self.spectrum[_].GetXaxis())


        # Particle Content
        #---------------------------------------------------------------------
        name   = settings.Content['name']
        title  = settings.Content['title']
        ytitle = settings.Content['ytitle']
        title  = '{};;{}'.format(title, ytitle)
        stats  = settings.Content['stats']
        xbins  = settings.Content['xbins']
        
        self.content = R.TH1D(name, title, xbins, 1, xbins+1)
        self.content.GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
        self.content.GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
        #bin_labels = settings.Content['bin_labels']
        #for _ in range(1, xbins+1):
        #   self.content.GetXaxis().SetBinLabel(_, bin_labels[_])


        # First Impact
        #---------------------------------------------------------------------
        name   = settings.Impact['name']
        title  = settings.Impact['title']
        ytitle = settings.Impact['ytitle']
        title  = '{};;{}'.format(title, ytitle)
        stats  = settings.Impact['stats']
        xbins  = settings.Impact['xbins']
        ybins  = settings.Impact['ybins']
        xmin   = settings.Impact['xmin']
        ymin   = settings.Impact['ymin']
        xmax   = settings.Impact['xmax']
        ymax   = settings.Impact['ymax']

        self.impact = R.TH2D(name, title, xbins, xmin, xmax, ybins, ymin, ymax)
        self.impact.GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
        self.impact.GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
        self.__log_binning__(self.impact.GetYaxis())
        #bin_labels = settings.Impact['bin_labels']
        #for _ in range(1, xbins+1):
        #   self.impact.GetXaxis().SetBinLabel(_, bin_labels[_])

    ##########################################################################

    def __print_std_err__(self):
        print(' [ERROR] Sorry, data file was not successfully loaded.')
        print("         You'll need to try again from scratch.")


