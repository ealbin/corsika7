#!/usr/bin/env python

import ROOT as R

import glob
import os

__catagories__ = ['photon',
                  'electron',
                  'muon',
                  'proton',
                  'neutron',
                  'ocharged',
                  'oneutral',
                  'nuclei']
__symbols__ = {}
__symbols__['photon']   = '#gamma'
__symbols__['electron'] = 'e^{#pm}'
__symbols__['muon']     = '#mu^{#pm}'
__symbols__['proton']   = 'p, #bar{p}'
__symbols__['neutron']  = 'n, #bar{n}'
__symbols__['ocharged'] = 'other^{#pm}'
__symbols__['oneutral'] = 'other^{0}'
__symbols__['nuclei']   = 'nuclei'

__levels__ = ['all', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

def make1Dcanvas(name, title='', width=800, height=800):
    left   = 0.12
    right  = 0.05
    bottom = 0.12
    top    = 0.10
    canvas = R.TCanvas(name, title, width, height)
    canvas.SetFillStyle(4000)
    canvas.SetMargin(left, right, bottom, top)
    return canvas


def make1Dlegend():
    x1 = .6
    x2 = .94
    y1 = .7
    y2 = .89
    legend = R.TLegend(x1, y1, x2, y2)
    legend.SetFillStyle(4000)
    legend.SetBorderSize(0)
    return legend


def thin():
    subfolder = 'thin_plots'

    dirs = {}
    for f in glob.glob('**/*.root', recursive=True):
        dirname = os.path.dirname(f)
        if (dirname not in dirs):
            dirs[dirname] = []
        dirs[dirname].append(os.path.basename(f))

    cdensity    = make1Dcanvas('density')
    cspectrum   = make1Dcanvas('spectrum')
    ccontent    = make1Dcanvas('content')
    cefficiency = make1Dcanvas('efficiency')

    cdensity.SetLogx()
    cspectrum.SetLogx()

    cdensity.SetLogy()
    cspectrum.SetLogy()
    ccontent.SetLogy()
    
    for dirname in dirs:
        savedir = os.path.join(dirname, subfolder)
        if (not os.path.isdir(savedir)):
            os.mkdir(savedir)

        tfiles = [[] for _ in range(4)]
        for basename in dirs[dirname]:
            tfile = R.TFile(os.path.join(dirname, basename))

            thinning = tfile.GetName().split('_')[-1].split('.')[0]
            if (thinning == '1E-6'):
                index  = 0
                marker = R.kFullDotMedium #R.kFullCircle
            elif (thinning == '1E-4'):
                index  = 1
                marker = R.kFullCross
            elif (thinning == '1E-2'):
                index  = 2
                marker = R.kOpenCrossX
            elif (thinning == '1'):
                index  = 3
                marker = R.kOpenCircle
                thinning == '1E-0'
            
            tfiles[index] = [tfile, marker, thinning]

   
        # EFFICIENCY
        cefficiency.Clear()
        lefficiency = make1Dlegend()
        efficiency_count = 0
        for tfile in tfiles:
            if (tfile == []):
                continue

            tfile[0].cd()

            cefficiency.cd()
            h = tfile[0].Get('_efficiency')
            h.SetMarkerStyle(tfile[1])
            lefficiency.AddEntry(h, 'thin {}'.format(tfile[2]), 'pl')
            h.SetMinimum(0.)
            h.SetMaximum(1.05)
            if (efficiency_count == 0):
                h.Draw('hist pl')
            else:
                h.Draw('hist pl same')
            efficiency_count += 1

        if (efficiency_count > 0):
            cefficiency.cd()
            lefficiency.Draw()
            cefficiency.Update()
            cefficiency.SaveAs('{}/efficiency.png')
        
        
        for level in __levels__:

            # CONTENT
            ccontent.Clear()
            lcontent = make1Dlegend()
            content_count = 0
            for tfile in tfiles:
                if (tfile == []):
                    continue

                tfile[0].cd()

                ccontent.cd()
                h = tfile[0].Get('_content_{}'.format(level))
                h.SetMarkerStyle(tfile[1])
                lcontent.AddEntry(h, 'thin {}'.format(tfile[2]), 'p')
                if (content_count == 0):
                    min_  = None
                    max_  = None
                    xaxis = h.GetXaxis()
                    for bin_ in range(1, xaxis.GetNbins() + 1):
                        v = h.GetBinContent(bin_)
                        if (v > 0):
                            if (min_ is None or max_ is None):
                                min_ = v
                                max_ = v
                            else:
                                if (v < min_):
                                    min_ = v
                                if (v > max_):
                                    max_ = v
                    if (min_ is not None and max_ is not None):
                        h.SetMaximum(max_ * 100.)
                        h.SetMinimum(0.1)
                        h.Draw('hist p')
                else:
                    h.Draw('hist p same')
                content_count += 1
            
            if (content_count > 0):
                ccontent.cd()
                lcontent.Draw()
                ccontent.Update()
                ccontent.SaveAs('{}/content_{}.png'.format(savedir, level))


            for catagory in __catagories__:
                cdensity.Clear()
                cspectrum.Clear()
                
                ldensity  = make1Dlegend()
                lspectrum = make1Dlegend()

                symbol = __symbols__[catagory]

                density_count  = 0
                spectrum_count = 0

                for tfile in tfiles:
                    if (tfile == []):
                        continue

                    tfile[0].cd()
         
                    # DENSITY
                    cdensity.cd()
                    h = tfile[0].Get('_density_{}_{}'.format(catagory, level))
                    h.SetMarkerStyle(tfile[1])
                    ldensity.AddEntry(h, '{} thin {}'.format(symbol, tfile[2]), 'p')
                    if (density_count == 0):
                        min_  = None
                        max_  = None
                        xaxis = h.GetXaxis()
                        for bin_ in range(1, xaxis.GetNbins() + 1):
                            v = h.GetBinContent(bin_)
                            if (v > 0):
                                if (min_ is None or max_ is None):
                                    min_ = v
                                    max_ = v
                                else:
                                    if (v < min_):
                                        min_ = v
                                    if (v > max_):
                                        max_ = v
                        if (min_ is None or max_ is None):
                            continue
                        h.SetMaximum(max_ * 100.)
                        h.SetMinimum(min_ / 100.)
                        h.Draw('hist p')
                    else:
                        h.Draw('hist p same')
                    density_count += 1

                    # SPECTRUM
                    cspectrum.cd()
                    h = tfile[0].Get('_spectrum_{}_{}'.format(catagory, level))
                    h.SetMarkerStyle(tfile[1])
                    lspectrum.AddEntry(h, '{} thin {}'.format(symbol, tfile[2]), 'p')
                    if (spectrum_count == 0):
                        min_  = None
                        max_  = None
                        xaxis = h.GetXaxis()
                        for bin_ in range(1, xaxis.GetNbins() + 1):
                            v = h.GetBinContent(bin_)
                            if (v > 0):
                                if (min_ is None or max_ is None):
                                    min_ = v
                                    max_ = v
                                else:
                                    if (v < min_):
                                        min_ = v
                                    if (v > max_):
                                        max_ = v
                        h.SetMaximum(max_ * 100.)
                        h.SetMinimum(min_ / 100.)
                        h.Draw('hist p')
                    else:
                        h.Draw('hist p same')
                    spectrum_count += 1

                if (density_count > 0):
                    cdensity.cd()
                    ldensity.Draw()
                    cdensity.Update()
                    cdensity.SaveAs('{}/density_{}_{}.png'.format(savedir, catagory, level))

                if (spectrum_count > 0):
                    cspectrum.cd()
                    lspectrum.Draw()
                    cspectrum.Update()
                    cspectrum.SaveAs('{}/spectrum_{}_{}.png'.format(savedir, catagory, level))


