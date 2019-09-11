#!/usr/bin/env python

import os
import ROOT as R

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

def make1Dcanvas(name, title='', width=600, height=600):
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

def has_key(tfile, key):
    n_keys = tfile.GetNkeys()
    keys   = tfile.GetListOfKeys()
    for _ in range(n_keys):
        key_name = keys.At(_).GetName()
        if (key_name == key):
            return True
    return False


def min_max(hist):
    min_  = None
    max_  = None
    xaxis = hist.GetXaxis()
    for bin_ in range(1, xaxis.GetNbins() + 1):
        v = hist.GetBinContent(bin_)
        if (v > 0):
            if (min_ is None or max_ is None):
                min_ = v
                max_ = v
            else:
                if (v < min_):
                    min_ = v
                if (v > max_):
                    max_ = v
    return (min_, max_)


def glob_files(root, startswith='*'):
    filelist = []
    root.cd()
    path = root.GetPath()
    for key in root.GetListOfKeys():
        if (key.IsFolder()):
            filelist = filelist + glob_files(root.GetDirectory(key.GetName()))
        elif (startswith == '*'):
            filelist.append(path + '/' + key.GetName())
        elif (key.GetName().startswith(startswith)):
            filelist.append(path + '/' + key.GetName())
    return filelist


def thin(infile, outfolder='thin_plots', batchmode=True):

    if (batchmode == True):
        R.gROOT.SetBatch(R.kTRUE)

    thin8_linestyle = 2
    thin8_linewidth = 1

    other_linestyle = 1
    other_linewidth = 1

    cdensity    = make1Dcanvas('density')
    cspectrum   = make1Dcanvas('spectrum')
    ccontent    = make1Dcanvas('content')
    cefficiency = make1Dcanvas('efficiency')

    cdensity.SetLogx()
    cspectrum.SetLogx()

    cdensity.SetLogy()
    cspectrum.SetLogy()
    ccontent.SetLogy()
    
    f = R.TFile(infile)
    filelist = glob_files(f, startswith='THIN_')
    filemap = {}
    for _ in filelist:
        dirname  = os.path.dirname(_)
        basename = os.path.basename(_)
        if (dirname not in filemap):
            filemap[dirname] = []
        filemap[dirname].append(basename)

    for dirname in filemap.keys():

        subfolder = dirname.split(':')[-1]
        if (subfolder[0] == '/'):
            subfolder = subfolder[1:]
        savedir = os.path.join(outfolder, subfolder)
        if (not os.path.isdir(savedir)):
            os.makedirs(savedir, exist_ok=True)
        
        def get_histlist(kind, level=None, catagory=None):
            histlist = [[] for _ in range(5)]
            for histkey in filemap[dirname]:
                tokens = histkey.split('_')

                if (len(tokens) < 4):
                    continue

                if (tokens[3] != kind):
                    continue
                
                if (tokens[0] != 'THIN'):
                    continue

                if (tokens[2] != '1'):
                    continue

                if (level is None and catagory is None):
                    pass
                elif (level is not None and catagory is None):
                    if (tokens[4] != level):
                        continue
                elif (level is not None and catagory is not None):
                    if (tokens[4] != catagory or tokens[5] != level):
                        continue
                else:
                    if (tokens[4] != catagory):
                        continue

                thinning = tokens[1]
                if (thinning == '1E-8'):
                    index  = 0
                    marker = None
                elif (thinning == '1E-6'):
                    index  = 1
                    marker = R.kFullDotMedium
                elif (thinning == '1E-4'):
                    index  = 2
                    marker = R.kFullCross
                elif (thinning == '1E-2'):
                    index  = 3
                    marker = R.kOpenCrossX
                elif (thinning == '1E-0'):
                    index  = 4
                    marker = R.kOpenCircle
                else:
                    continue

                histlist[index] = [os.path.join(dirname, histkey), marker, thinning]
            return histlist

        # EFFICIENCY
        hfiles = get_histlist('efficiency')
        cefficiency.Clear()
        lefficiency = make1Dlegend()
        efficiency_count = 0
        for hfile in hfiles:
            if (hfile == []):
                continue
            h = f.Get(hfile[0])
            
            cefficiency.cd()
            if (hfile[1] is not None):
                h.SetMarkerStyle(hfile[1])
                h.SetLineStyle(other_linestyle)
                h.SetLineWidth(other_linewidth)
                plot_opt = 'p'
            else:
                h.SetMarkerStyle(R.kFullDotSmall)
                h.SetLineStyle(thin8_linestyle)
                h.SetLineWidth(thin8_linewidth)
                plot_opt = 'l'
            lefficiency.AddEntry(h, 'thin {}'.format(hfile[2]), plot_opt)
            h.SetMinimum(0.)
            h.SetMaximum(1.05)
            if (efficiency_count == 0):
                h.Draw('hist ' + plot_opt)
            else:
                h.Draw('hist ' + plot_opt + ' same')
            efficiency_count += 1

        if (efficiency_count > 0):
            cefficiency.cd()
            lefficiency.Draw()
            cefficiency.Update()
            cefficiency.SaveAs('{}/efficiency.png'.format(savedir))
        
        
        for level in __levels__:

            # CONTENT
            hfiles = get_histlist('content', level=level)
            ccontent.Clear()
            lcontent = make1Dlegend()
            content_count = 0
            for hfile in hfiles:
                if (hfile == []):
                    continue
                h = f.Get(hfile[0])

                ccontent.cd()
                if (hfile[1] is not None):
                    h.SetMarkerStyle(hfile[1])
                    h.SetLineStyle(other_linestyle)
                    h.SetLineWidth(other_linewidth)
                    plot_opt = 'p'
                else:
                    h.SetMarkerStyle(R.kFullDotSmall)
                    h.SetLineStyle(thin8_linestyle)
                    h.SetLineWidth(thin8_linewidth)
                    plot_opt = 'l'
                lcontent.AddEntry(h, 'thin {}'.format(hfile[2]), plot_opt)
                if (content_count == 0):
                    h.SetMinimum(0.1)
                    min_, max_ = min_max(h)
                    if (max_ is not None):
                        h.SetMaximum(max_ * 100.)
                        h.Draw('hist ' + plot_opt)
                        content_count = 1
                else:
                    h.Draw('hist ' + plot_opt + ' same')
                    content_count += 1
            
            if (content_count > 0):
                ccontent.cd()
                lcontent.Draw()
                ccontent.Update()
                ccontent.SaveAs('{}/content_{}.png'.format(savedir, level))


            for catagory in __catagories__:

                symbol = __symbols__[catagory]
                
                # DENSITY
                hfiles = get_histlist('density', level=level, catagory=catagory)
                for _ in hfiles:
                    print(_)
                cdensity.Clear()
                ldensity = make1Dlegend()
                density_count  = 0
                for hfile in hfiles:
                    if (hfile == []):
                        continue
                    h = f.Get(hfile[0])
                    
                    cdensity.cd()
                    if (hfile[1] is not None):
                        h.SetMarkerStyle(hfile[1])
                        h.SetLineStyle(other_linestyle)
                        h.SetLineWidth(other_linewidth)
                        plot_opt = 'p'
                    else:
                        h.SetMarkerStyle(R.kFullDotSmall)
                        h.SetLineStyle(thin8_linestyle)
                        h.SetLineWidth(thin8_linewidth)
                        plot_opt = 'l'
                    ldensity.AddEntry(h, '{} thin {}'.format(symbol, hfile[2]), plot_opt)
                    if (density_count == 0):
                        min_, max_ = min_max(h)
                        if (min_ is not None and max_ is not None):
                            h.SetMaximum(max_ * 100.)
                            h.SetMinimum(min_ / 100.)
                            h.Draw('hist ' + plot_opt)
                            density_count = 1
                    else:
                        h.Draw('hist ' + plot_opt + ' same')
                        density_count += 1

                if (density_count > 0):
                    cdensity.cd()
                    ldensity.Draw()
                    cdensity.Update()
                    cdensity.SaveAs('{}/density_{}_{}.png'.format(savedir, catagory, level))
                
                
                # SPECTRUM
                hfiles = get_histlist('spectrum', level=level, catagory=catagory)
                cspectrum.Clear()
                lspectrum = make1Dlegend()
                spectrum_count = 0
                for hfile in hfiles: 
                    if (hfile == []):
                        continue
                    h = f.Get(hfile[0])

                    cspectrum.cd()
                    if (hfile[1] is not None):
                        h.SetMarkerStyle(hfile[1])
                        h.SetLineStyle(other_linestyle)
                        h.SetLineWidth(other_linewidth)
                        plot_opt = 'p'
                    else:
                        h.SetMarkerStyle(R.kFullDotSmall)
                        h.SetLineStyle(thin8_linestyle)
                        h.SetLineWidth(thin8_linewidth)
                        plot_opt = 'l'
                    lspectrum.AddEntry(h, '{} thin {}'.format(symbol, hfile[2]), plot_opt)
                    if (spectrum_count == 0):
                        min_, max_ = min_max(h)
                        if (min_ is not None and max_ is not None):
                            h.SetMaximum(max_ * 100.)
                            h.SetMinimum(min_ / 100.)
                            h.Draw('hist ' + plot_opt)
                            spectrum_count = 1
                    else:
                        h.Draw('hist ' + plot_opt + ' same')
                        spectrum_count += 1

                if (spectrum_count > 0):
                    cspectrum.cd()
                    lspectrum.Draw()
                    cspectrum.Update()
                    cspectrum.SaveAs('{}/spectrum_{}_{}.png'.format(savedir, catagory, level))
