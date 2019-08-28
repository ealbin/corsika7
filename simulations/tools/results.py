
import os
import numpy as np
import ROOT as R

from . import nuclei
from . import particle
from . import settings

# TODO: replace this with load
__data__ = R.TFile('./proton-1e14/DAT000001.root')


##############################################################################

def relative_diff_1D(hist1, hist2):
    diff = hist1.Clone()
    diff.Add(hist2, hist1, -1.)
    for _ in range(1, 1 + diff.GetNbinsX()):
        h1 = hist1.GetBinContent(_)
        if (h1 == 0.):
            diff.SetBinContent(_, 0.)
        else:
            d  = diff.GetBinContent(_)
            diff.SetBinContent(_, d / h1)
    return diff


##############################################################################

def relative_diff_2D(hist1, hist2):
    diff = hist1.Clone()
    diff.Add(hist2, hist1, -1.)
    for i in range(1, 1 + diff.GetNbinsX()):
        for j in range(1, 1 + diff.GetNbinsY()):
            h1 = hist1.GetBinContent(i, j)
            if (h1 == 0.):
                diff.SetBinContent(i, j, 0.)
            else:
                d  = diff.GetBinContent(i, j)
                diff.SetBinContent(i, j, d / h1)
    return diff


##############################################################################

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
                continue

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

        try:
            self.legend
        except:
            pass
        else:
            del self.legend
        
        name   = 'density_canvas'
        width  = settings.Density['width']
        height = settings.Density['height']
        self.density_canvas = Results.__1D_canvas__(
                name, title=title, width=width, height=height)
        self.density_canvas.SetLogx()
        self.density_canvas.SetLogy()

        for _ in settings.Catagories:
            self.density[_].SetMinimum(np.power(10., settings.Density['ymin']))
            self.density[_].SetMaximum(np.power(10., settings.Density['ymax']))
            self.density[_].Draw('hist pe same')

        self.legend = Results.__1D_legend__()
        self.legend.AddEntry(self.density['photon'],   '#gamma',      'pe')
        self.legend.AddEntry(self.density['proton'],   'p, #bar{p}',  'pe')
        self.legend.AddEntry(self.density['electron'], 'e^{#pm}',     'pe')
        self.legend.AddEntry(self.density['neutron'],  'n, #bar{n}',  'pe')
        self.legend.AddEntry(self.density['muon'],     '#mu^{#pm}',   'pe')
        self.legend.AddEntry(self.density['ocharged'], 'other^{#pm}', 'pe')
        self.legend.AddEntry(self.density['nuclei'],   'nuclei',      'pe')
        self.legend.AddEntry(self.density['oneutral'], 'other^{0}',   'pe')
        self.legend.Draw()
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
                continue

            # "x" is West-to-East (i.e. -y)
            # "y" is South-to-North (i.e. x)
            self.slice[catagory].Fill(-y, x)

        try:
            self.slice_canvases
        except:
            pass
        else:
            del self.slice_canvases

        try:
            self.legend
        except:
            pass
        else:
            del self.legend
        
        R.gStyle.SetPalette(settings.Palette)
        name   = 'slice_canvas'
        width  = settings.Slice['width']
        height = settings.Slice['height']
        self.slice_canvases = {}
        self.legend = []
        for _ in settings.Catagories:
            
            # skip empty plots
            if (self.slice[_].GetEntries() == 0.):
                continue

            canvas = Results.__2D_canvas__(
                    name + '_' + _, title=title, width=width, height=height)
            canvas.SetLogz()
            self.slice_canvases[_] = canvas

            self.slice[_].SetMinimum(np.power(10., settings.Slice['zmin']))
            self.slice[_].SetMaximum(np.power(10., settings.Slice['zmax']))
            self.slice[_].Draw('colz')
            
            legend = Results.__2D_legend__()
            legend.AddEntry(self.slice[_], settings.Symbols[_], '')
            legend.Draw()
            self.legend.append(legend)
            canvas.Update()

    ##########################################################################

    def plot_spectrum(self, obs_level=-1):
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
            self.spectrum[_].SetTitle(title)
            self.spectrum[_].Reset('ICEM')

        for _ in range(self.sim.particle__):
            
            if (obs_level != -1):
                if (self.sim.particle__ObservationLevel[_] != obs_level):
                    continue

            pid = self.sim.particle__ParticleID[_]
            px  = self.sim.particle__Px[_] # GeV
            py  = self.sim.particle__Py[_] # GeV
            pz  = self.sim.particle__Pz[_] # GeV

            if (str(pid) in particle.ID):
                catagory = particle.ID[str(pid)]['catagory']
                mass     = particle.ID[str(pid)]['mass'] # GeV
            elif (str(pid)[-2:] in nuclei.ID):
                catagory = nuclei.ID[str(pid)[-2:]]['catagory']
                mass     = nuclei.ID[str(pid)[-2:]]['mass'] # GeV
            else:
                print(' [ERROR] Unknown particle ID: {}'.format(pid))
                continue

            energy = np.sqrt(px*px + py*py + pz*pz + mass*mass) * 1e9 # eV
            self.spectrum[catagory].Fill(energy)

        xaxis = self.spectrum[settings.Catagories[0]].GetXaxis()
        for _ in range(xaxis.GetNbins()):
            i     = _ + 1
            width = xaxis.GetBinWidth(i)

            for cat in settings.Catagories:
                raw = self.spectrum[cat].GetBinContent(_)
                self.spectrum[cat].SetBinContent(_, raw / width)

                raw = self.spectrum[cat].GetBinError(_)
                self.spectrum[cat].SetBinError(_, raw / width)

        try:
            self.spectrum_canvas
        except:
            pass
        else:
            del self.spectrum_canvas

        try:
            self.legend
        except:
            pass
        else:
            del self.legend
        
        name   = 'spectrum_canvas'
        width  = settings.Density['width']
        height = settings.Density['height']
        self.spectrum_canvas = Results.__1D_canvas__(
                name, title=title, width=width, height=height)
        self.spectrum_canvas.SetLogx()
        self.spectrum_canvas.SetLogy()

        for _ in settings.Catagories:
            self.spectrum[_].SetMinimum(np.power(10., settings.Spectrum['ymin']))
            self.spectrum[_].SetMaximum(np.power(10., settings.Spectrum['ymax']))
            self.spectrum[_].Draw('hist pe same')

        self.legend = Results.__1D_legend__()
        self.legend.AddEntry(self.density['photon'],   '#gamma',      'pe')
        self.legend.AddEntry(self.density['proton'],   'p, #bar{p}',  'pe')
        self.legend.AddEntry(self.density['electron'], 'e^{#pm}',     'pe')
        self.legend.AddEntry(self.density['neutron'],  'n, #bar{n}',  'pe')
        self.legend.AddEntry(self.density['muon'],     '#mu^{#pm}',   'pe')
        self.legend.AddEntry(self.density['ocharged'], 'other^{#pm}', 'pe')
        self.legend.AddEntry(self.density['nuclei'],   'nuclei',      'pe')
        self.legend.AddEntry(self.density['oneutral'], 'other^{0}',   'pe')
        self.legend.Draw()
        self.spectrum_canvas.Update()

    ##########################################################################

    def plot_content(self, obs_level=-1):
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

        self.content.SetTitle(title)
        self.content.Reset('ICEM')

        for _ in range(self.sim.particle__):
            
            if (obs_level != -1):
                if (self.sim.particle__ObservationLevel[_] != obs_level):
                    continue

            pid = self.sim.particle__ParticleID[_]

            if (str(pid) in particle.ID):
                catagory = particle.ID[str(pid)]['catagory']
            elif (str(pid)[-2:] in nuclei.ID):
                catagory = nuclei.ID[str(pid)[-2:]]['catagory']
            else:
                print(' [ERROR] Unknown particle ID: {}'.format(pid))
                continue

            center = self.content.GetBinCenter(settings.Content[catagory])
            self.content.Fill(center)

        try:
            self.content_canvas
        except:
            pass
        else:
            del self.content_canvas

        name   = 'content_canvas'
        width  = settings.Content['width']
        height = settings.Content['height']
        self.content_canvas = Results.__1D_canvas__(
                name, title=title, width=width, height=height)
        self.content_canvas.SetLogy()

        self.content.SetMinimum(np.power(10., settings.Content['ymin']))
        self.content.SetMaximum(np.power(10., settings.Content['ymax']))
        self.content.Draw('hist pe same')
        
        self.content_canvas.Update()

    ##########################################################################

    def plot_impact(self):
        if (not self.all_good):
            self.__print_std_err__()
            return

        title = settings.Impact['title'] 
        self.impact.SetTitle(title)
        self.impact.Reset('ICEM')

        x = self.sim.shower_FirstTarget
        y = self.sim.shower_FirstHeight / 1e5 # [km]
        self.impact.Fill(x, y)

        try:
            self.impact_canvas
        except:
            pass
        else:
            del self.impact_canvas

        R.gStyle.SetPalette(settings.Palette)
        name   = 'impact_canvas'
        width  = settings.Impact['width']
        height = settings.Impact['height']
        self.impact_canvas = Results.__2D_canvas__(
                name, title=title, width=width, height=height)
        self.impact_canvas.SetLogz()

        self.impact.SetMinimum(np.power(10., settings.Impact['zmin']))
        self.impact.SetMaximum(np.power(10., settings.Impact['zmax']))
        self.impact.Draw('colz')
        
        self.impact_canvas.Update()

    ##########################################################################

    def plot_efficiency(self):
        if (not self.all_good):
            self.__print_std_err__()
            return

        self.efficiency.Reset('ICEM')
        max_energy = self.sim.shower_Energy # GeV

        nlevels = self.run.run_ObservationLevel.size()
        for _ in range(nlevels):
            xbin = _ + 1
            obs_level = _ + 1
            if (_ == 9):
                obs_level = 0

            energy_GeV = 0.
            for i in range(self.sim.particle__):
            
                if (self.sim.particle__ObservationLevel[i] != obs_level):
                    continue

                pid = self.sim.particle__ParticleID[i]
                px  = self.sim.particle__Px[i] # GeV
                py  = self.sim.particle__Py[i] # GeV
                pz  = self.sim.particle__Pz[i] # GeV

                if (str(pid) in particle.ID):
                    mass     = particle.ID[str(pid)]['mass'] # GeV
                elif (str(pid)[-2:] in nuclei.ID):
                    mass     = nuclei.ID[str(pid)[-2:]]['mass'] # GeV
                else:
                    print(' [ERROR] Unknown particle ID: {}'.format(pid))
                    continue

                energy_GeV += np.sqrt(px*px + py*py + pz*pz + mass*mass) # GeV

            efficiency = energy_GeV / max_energy
            self.efficiency.SetBinContent(xbin, efficiency)

        try:
            self.efficiency_canvas
        except:
            pass
        else:
            del self.efficiency_canvas

        name   = 'efficiency_canvas'
        width  = settings.Efficiency['width']
        height = settings.Efficiency['height']
        self.efficiency_canvas = Results.__1D_canvas__(
                name, title='Energy Efficiency', width=width, height=height)

#        self.efficiency.SetMinimum(settings.Efficiency['ymin'])
#        self.efficiency.SetMaximum(settings.Efficiency['ymax'])
        self.efficiency.Draw('hist pl')
        
        self.efficiency_canvas.Update()


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
            self.spectrum[_].SetMarkerColor(settings.Colors[_])
            self.spectrum[_].SetMarkerStyle(R.kFullCircle)
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
        self.content.SetStats(stats)
        self.content.SetMarkerColor(R.kBlack)
        self.content.SetMarkerStyle(R.kFullCircle)
        self.content.SetLineColor(R.kBlack)
        self.content.GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
        self.content.GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
        
        axis = self.content.GetXaxis()
        for _ in settings.Catagories:
            axis.SetBinLabel(settings.Content[_], settings.Symbols[_])


        # First Impact
        #---------------------------------------------------------------------
        name   = settings.Impact['name']
        title  = settings.Impact['title']
        xtitle = settings.Impact['xtitle']
        ytitle = settings.Impact['ytitle']
        title  = '{};{};{}'.format(title, xtitle, ytitle)
        stats  = settings.Impact['stats']
        xbins  = settings.Impact['xbins']
        ybins  = settings.Impact['ybins']
        ymin   = settings.Impact['ymin']
        ymax   = settings.Impact['ymax']

        self.impact = R.TH2D(name, title, xbins, 0, xbins, ybins, ymin, ymax)
        self.impact.SetStats(stats)
        self.impact.SetMarkerColor(R.kBlack)
        self.impact.SetMarkerStyle(R.kFullCircle)
        self.impact.SetLineColor(R.kBlack)
        self.impact.GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
        self.impact.GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
        self.__log_binning__(self.impact.GetYaxis())
        
        axis = self.impact.GetXaxis()
        for _ in ['Random', 'Nitrogen', 'Oxygen', 'Argon']:
            axis.SetBinLabel(settings.Impact[_] + 1, _)


        # Energy Efficiency
        #---------------------------------------------------------------------
        name   = settings.Efficiency['name']
        title  = settings.Efficiency['title']
        xtitle = settings.Efficiency['xtitle']
        ytitle = settings.Efficiency['ytitle']
        title  = '{};{};{}'.format(title, xtitle, ytitle)
        stats  = settings.Efficiency['stats']
        xbins  = settings.Efficiency['xbins']
        ymin   = settings.Efficiency['ymin']
        ymax   = settings.Efficiency['ymax']

        self.efficiency = R.TH1D(name, title, xbins, 1, xbins+1)
        self.efficiency.SetStats(stats)
        self.efficiency.SetMarkerColor(R.kBlack)
        self.efficiency.SetMarkerStyle(R.kFullCircle)
        self.efficiency.SetLineColor(R.kBlack)
        self.efficiency.GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
        self.efficiency.GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
        
        axis    = self.efficiency.GetXaxis()
        nlevels = self.run.run_ObservationLevel.size()
        for _ in range(nlevels):
            xbin  = _ + 1
            index = _ + 1
            if (_ == 9):
                index = 0
            altitude = self.run.run_ObservationLevel.at(_) / 1e5
            altstr   = '{:.2f}'.format(altitude).rstrip('0').rstrip('.')
            axis.SetBinLabel(xbin, altstr)


    ##########################################################################

    def __print_std_err__(self):
        print(' [ERROR] Sorry, data file was not successfully loaded.')
        print("         You'll need to try again from scratch.")


    ##########################################################################

    def __1D_canvas__(name, title='', width=800, height=800):
        left   = settings.Margins['left']
        right  = settings.Margins['right']
        bottom = settings.Margins['bottom']
        top    = settings.Margins['top']
        canvas = R.TCanvas(name, title, width, height)
        canvas.SetFillStyle(4000)
        canvas.SetMargin(left, right, bottom, top)
        return canvas


    ##########################################################################

    def __2D_canvas__(name, title='', width=800, height=800):
        left   = settings.Margins['left']
        right  = settings.Margins['right2D']
        bottom = settings.Margins['bottom']
        top    = settings.Margins['top']
        canvas = R.TCanvas(name, title, width, height)
        canvas.SetFillStyle(4000)
        canvas.SetMargin(left, right, bottom, top)
        return canvas


    ##########################################################################

    def __1D_legend__():
        x1 = .6
        x2 = .94
        y1 = .7
        y2 = .89
        legend = R.TLegend(x1, y1, x2, y2)
        legend.SetFillStyle(4000)
        legend.SetBorderSize(0)
        legend.SetNColumns(2)
        return legend

    ##########################################################################

    def __2D_legend__():
        x1 = .6
        x2 = .8
        y1 = .8
        y2 = .89
        legend = R.TLegend(x1, y1, x2, y2)
        legend.SetFillStyle(4000)
        legend.SetBorderSize(0)
        return legend



