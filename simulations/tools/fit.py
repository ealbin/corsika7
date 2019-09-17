
import os
import numpy as np
import pandas as pd
from scipy.special import gamma
import ROOT as R
from . import plot

import matplotlib.pyplot as plt

def A_func(A):
    return np.log(A)

def E_func(E):
    return np.log10(E)

def H_func(H):
    return H

class NKG:
    name     = 'NKG'
    fitfunc  = '10^[0] * x^(-[1]) * (1. + x/[2])^(-[3])'
    equation = '10^{a_{0}} x^{-a_{1}} (1 + x / a_{2})^{-a_{3}}'

    a0_val =  -5.
    a0_low = -10.
    a0_hi  =  10.

    a1_val =   1.
    a1_low =  -1.
    a1_hi  =   5.

    a2_val =   1.
    a2_low =  -2.
    a2_hi  =  20.

    a3_val =   5.
    a3_low =  -2.
    a3_hi  = 100.


class Albin:
    name     = 'Albin'
    fitfunc  = 'x^(-[0]) * exp(-[1] - [2] * x^[3])'
    equation = 'x^{-a_{0}} exp(-a_{1} - a_{2} x^{a_{3}})' 

    a0_val =   0.5
    a0_low = -10.
    a0_hi  =  20.

    a1_val =   5.
    a1_low = -30.
    a1_hi  =  30.

    a2_val =   5.
    a2_low = -10.
    a2_hi  =  35.

    a3_val =  0.5
    a3_low = -5.
    a3_hi  =  5.

__XMIN__ = 1e-6
__XMAX__ = 1e3

class Density:

    def __init__(self):
        pass

    def plot_density(self, indir, mass_number, energy, altitude, channel='muon'):
        nkg = pd.read_csv(os.path.join(indir, '{}_{}_model.txt'.format(NKG.name, channel)), delim_whitespace=True)
        albin = pd.read_csv(os.path.join(indir, '{}_{}_model.txt'.format(Albin.name, channel)), delim_whitespace=True)

        fnkg = R.TF1('NKG',   NKG.fitfunc,   __XMIN__, __XMAX__)
        falb = R.TF1('Albin', Albin.fitfunc, __XMIN__, __XMAX__)

        A = A_func(mass_number)
        E = E_func(energy)
        H = H_func(altitude)

        X = np.asarray([1., A, E, H, A*A, A*E, A*H, E*E, E*H, H*H])
        for i in range(4):
            fnkg.SetParameter(i, np.dot(X, nkg.iloc[i,1:]))
            falb.SetParameter(i, np.dot(X, albin.iloc[i,1:]))

        self.canvas = plot.make1Dcanvas('canvas')
        self.canvas.SetLogx()
        self.canvas.SetLogy()

        legend = plot.make1Dlegend()
        fnkg.SetLineColor(R.kRed)
        falb.SetLineColor(R.kBlue)
        legend.AddEntry(fnkg, 'NKG', 'l')
        legend.AddEntry(falb, 'Albin', 'l')

        #fnkg.Draw()
        #falb.Draw('same')
        falb.Draw()
        legend.Draw()
        self.canvas.Update()


def likelihood_ratio(infile, indir, model0=NKG, model1=Albin, channel='muon'):
    model0_coeff = pd.read_csv(os.path.join(indir, '{}_{}_model.txt'.format(model0.name, channel)), delim_whitespace=True)
    model1_coeff = pd.read_csv(os.path.join(indir, '{}_{}_model.txt'.format(model1.name, channel)), delim_whitespace=True)
    
    alt_km = [0., 100., 50., 20., 10., 5., 2., 1.4, 1.0, .5]
    MassNumber = {'He':4, 'O':16, 'Fe':56}

    thin      = '1E-6'
    N         = '100'
    primaries = ['He', 'O', 'Fe']
    energies  = ['1e14', '1e15', '1e16', '1e17', '1e18', '1e19', '1e20', '1e21']
    altitudes = ['3', '4', '5', '6', '7', '8', '9', '0']
    models    = ['DPMJET', 'QGSJET', 'QGSII', 'SIBYLL', 'VENUS']
    
    f = R.TFile(infile)

#    def poisson(l, k):
#        return np.power(l, k) * np.exp(-l) / gamma(k + 1.)

    with open('{}/{}_llratio.txt'.format(indir, channel), 'w') as outfile:
        for primary in primaries:
            for energy in energies:
                for altitude in altitudes:

                    model0_sum = 0.
                    model1_sum = 0.

                    name = 'THIN_{}_{}_density_{}_{}'.format(thin, N, channel, altitude)
                    dpmjet = f.Get('{}/{}/{}/{}'.format(energy, 'DPMJET', primary, name))
                    qgsjet = f.Get('{}/{}/{}/{}'.format(energy, 'QGSJET', primary, name))
                    qgsii  = f.Get('{}/{}/{}/{}'.format(energy, 'QGSII',  primary, name))
                    sibyll = f.Get('{}/{}/{}/{}'.format(energy, 'SIBYLL', primary, name))
                    venus  = f.Get('{}/{}/{}/{}'.format(energy, 'VENUS',  primary, name))

                    f0 = R.TF1('model0', model0.fitfunc, __XMIN__, __XMAX__)
                    f1 = R.TF1('model1', model1.fitfunc, __XMIN__, __XMAX__)

                    A = A_func(MassNumber[primary])
                    E = E_func(float(energy))
                    H = H_func(alt_km[int(altitude)])
            
                    X = np.asarray([1., A, E, H, A*A, A*E, A*H, E*E, E*H, H*H])
                    for i in range(4):
                        f0.SetParameter(i, np.dot(X, model0_coeff.iloc[i,1:]))
                        f1.SetParameter(i, np.dot(X, model1_coeff.iloc[i,1:]))

                    for h in [dpmjet, qgsjet, qgsii, sibyll, venus]:
                        if (h != None):
                            xaxis = h.GetXaxis()
                            for bin in range(xaxis.GetNbins() + 1):
                                center = xaxis.GetBinCenterLog(bin)
                                binval = h.GetBinContent(bin)
                                
                                f0_val = f0.Eval(center)
                                f1_val = f1.Eval(center)

                                if (binval <= 0 or f0_val <= 0 or f1_val <= 0):
                                    continue

                                model0_sum += binval * np.log(f0_val) - f0_val - np.log(gamma(binval + 1.))
                                model1_sum += binval * np.log(f1_val) - f1_val - np.log(gamma(binval + 1.))

                                #model0_sum += np.log( poisson(f0.Eval(center), binval) )
                                #model1_sum += np.log( poisson(f1.Eval(center), binval) )

                    outfile.write('{:<20.10e}\n'.format(2*(model1_sum - model0_sum)))


def density_coeff(indir, model=Albin, channel='muon'):
    df = pd.read_csv(os.path.join(indir, '{}_coeff_data_{}.txt'.format(model.name, channel)), delim_whitespace=True)

    df.insert(0, 'C', 1)
    A = df.iloc[:,:10].to_numpy()
    B = np.linalg.pinv(A)

    with open(os.path.join(indir, '{}_{}_model.txt'.format(model.name, channel)), 'w') as outfile:
        outfile.write( ('{:15s}'*11).format('a', 'C', 'A', 'E', 'H', 'AA', 'AE', 'AH', 'EE', 'EH', 'HH') + '\n')
    
        for i in range(4):
            a = df.loc[:, 'a{}'.format(i)].to_numpy()
            string = '{:<15d}'.format(i)
            for coeff in np.dot(B, a):
                string += '{:<15.6e}'.format(coeff)
            outfile.write(string + '\n')


def quick_check(indir, model=Albin, channel='muon'):
    df = pd.read_csv(os.path.join(indir, '{}_coeff_data_{}.txt'.format(model.name, channel)), delim_whitespace=True)

    # no const
    A9 = df.iloc[:,:9].to_numpy()
    B9 = np.linalg.pinv(A9)
    
    # with const
    df.insert(0, 'C', 1)
    A10 = df.iloc[:,:10].to_numpy()
    B10 = np.linalg.pinv(A10)

    a   = []
    c9  = []
    c10 = []
    t9  = []
    t10 = []
    for i in range(4):
        a.append(df.loc[:, 'a{}'.format(i)].to_numpy())
        c9.append(np.dot(B9, a[i]))
        c10.append(np.dot(B10, a[i]))

        t9.append(np.dot(A9, c9[i]))
        t10.append(np.dot(A10, c10[i]))

    a   = np.asarray(a)
    c9  = np.asarray(c9)
    c10 = np.asarray(c10)
    t9  = np.asarray(t9)
    t10 = np.asarray(t10)

    plt.figure(figsize=[15, 15])
    for i in range(4):
        low = np.min([a[i], t9[i], t10[i]])
        hi  = np.max([a[i], t9[i], t10[i]])
        plt.subplot(2, 2, i+1)
        plt.scatter(a[i], t9[i],  color='b', label='no const')
        plt.scatter(a[i], t10[i], color='orange', label='with const')
        plt.plot([low,hi], [low,hi])
        plt.title('a{}'.format(i))
        plt.legend()

        d9  =  t9[i] - a[i]
        d10 = t10[i] - a[i]
        print('a{}, no const mean: {:7.3f} vs const mean: {:7.3f}'.format(i, d9.mean(), d10.mean()))
        print('a{}, no const std:  {:7.3f} vs const std:  {:7.3f}'.format(i, d9.std(), d10.std()))
        print()

    plt.gcf().suptitle('{} {}'.format(model.name, channel))


def density_check(infile, indir, model=Albin, savedir='fit_check', batchmode=True):

    if (batchmode == True):
        R.gROOT.SetBatch(R.kTRUE)

    if (not os.path.isdir(savedir)):
        os.makedirs(savedir, exist_ok=True)
    
    alt_km = [0., 100., 50., 20., 10., 5., 2., 1.4, 1.0, .5]
    MassNumber = {'He':4, 'O':16, 'Fe':56}

    thin      = '1E-6'
    N         = '100'
    primaries = ['He', 'O', 'Fe']
    energies  = ['1e14', '1e15', '1e16', '1e17', '1e18', '1e19', '1e20', '1e21']
    altitudes = ['3', '4', '5', '6', '7', '8', '9', '0']
    models    = ['DPMJET', 'QGSJET', 'QGSII', 'SIBYLL', 'VENUS']
    
    f = R.TFile(infile)
    canvas = plot.make1Dcanvas('check_canvas')
    canvas.SetLogy()
    canvas.SetLogx()
    
    for channel in ['muon', 'photon']:
        df = pd.read_csv(os.path.join(indir, '{}_coeff_data_{}.txt'.format(model.name, channel)), delim_whitespace=True)
        df.insert(0, 'C', 1.)
        A = df.iloc[:,:10].to_numpy()
        B = np.linalg.pinv(A)
        c = []
        for i in range(4):
            a = df.loc[:,'a{}'.format(i)].to_numpy()
            c.append(np.dot(B, a))

        for primary in primaries:
            for energy in energies:
                for altitude in altitudes:

                    name = 'THIN_{}_{}_density_{}_{}'.format(thin, N, channel, altitude)
                    dpmjet = f.Get('{}/{}/{}/{}'.format(energy, 'DPMJET', primary, name))
                    qgsjet = f.Get('{}/{}/{}/{}'.format(energy, 'QGSJET', primary, name))
                    qgsii  = f.Get('{}/{}/{}/{}'.format(energy, 'QGSII',  primary, name))
                    sibyll = f.Get('{}/{}/{}/{}'.format(energy, 'SIBYLL', primary, name))
                    venus  = f.Get('{}/{}/{}/{}'.format(energy, 'VENUS',  primary, name))

                    hists = []
                    wrong_lim = []
                    for h in [dpmjet, qgsjet, qgsii, sibyll, venus]:
                        if (h != None):
                            if (h.GetXaxis().GetXmin() > __XMIN__):
                                wrong_lim.append(h)
                            else:
                                hists.append(h)

                    if (len(wrong_lim) > len(hists)):
                        hists = wrong_lim

                    if (len(hists) > 0):
                        hsum = hists[0]
                        for h in hists[1:]:
                            hsum.Add(h)
                        hsum.Scale(1./float(len(hists)))

                        min_, max_ = plot.min_max(hsum)
                        if (min_ is not None and max_ is not None):
                            hsum.SetMaximum(max_ * 100.)
                            hsum.SetMinimum(min_ / 100.)
                        hsum.SetMarkerColor(R.kBlack)
                        hsum.SetMarkerStyle(R.kFullCircle)

                        xaxis = hsum.GetXaxis()
                        for stop_bin in range(xaxis.GetNbins(), 0, -1):
                            val = hsum.GetBinContent(stop_bin)
                            if (val > 0.):
                                break
                        limit = xaxis.GetBinCenter(stop_bin)

                        fitfunc = model.fitfunc
                        f1 = R.TF1('fitfunc', fitfunc, __XMIN__, limit)

                        A = A_func(MassNumber[primary])
                        E = E_func(float(energy))
                        H = H_func(alt_km[int(altitude)])
                
                        V = np.asarray([1., A, E, H, A*A, A*E, A*H, E*E, E*H, H*H])
                        for i in range(4):
                            f1.SetParameter(i, np.dot(V, c[i]))

                        f1.SetLineColor(R.kBlue)
                        f1.SetLineStyle(1)
                        f1.SetLineWidth(1)

                        legend = plot.make1Dlegend()
                        legend.AddEntry(hsum, '{}: n={}, thin={}, model ave'.format(plot.__symbols__[channel], N, thin), 'p')
                        legend.AddEntry(f1, model.equation, 'l')

                        canvas.cd()
                        canvas.Clear()
                        hsum.GetXaxis().SetRangeUser(__XMIN__, 1e3)
                        hsum.Draw('hist p')
                        f1.Draw('same')
                        legend.Draw()
                        canvas.Update()
                        canvas.SaveAs('{}/check_{}_{}_{}_{}_{}.png'.format(savedir, model.name, channel, primary, energy, altitude))


def supervised_density(infile, savedir='fit_density', model=Albin):

    if (not os.path.isdir(savedir)):
        os.makedirs(savedir, exist_ok=True)

    thin      = '1E-6'
    N         = '100'
    channels  = ['muon', 'photon']
    primaries = ['He', 'O', 'Fe']
    energies  = ['1e14', '1e15', '1e16', '1e17', '1e18', '1e19', '1e20', '1e21']
    altitudes = ['3', '4', '5', '6', '7', '8', '9', '0'] 
    models    = ['DPMJET', 'QGSJET', 'QGSII', 'SIBYLL', 'VENUS']

    alt_km = [0., 100., 50., 20., 10., 5., 2., 1.4, 1.0, .5]
    MassNumber = {'He':4, 'O':16, 'Fe':56}

    f = R.TFile(infile)
    canvas = plot.make1Dcanvas('fitcanvas')
    canvas.SetLogy()
    canvas.SetLogx()
   
    def write_header(file):
        file.write( ('{:15s}'*13).format('A', 'E', 'H', 'AA', 'AE', 'AH', 'EE', 'EH', 'HH', 'a0', 'a1', 'a2', 'a3') + '\n')

    def write_data(file, primary, energy, altitude, a0, a1, a2, a3, flush=True):
        A = A_func(MassNumber[primary])
        E = E_func(float(energy))
        H = H_func(alt_km[int(altitude)])        
        file.write(('{:<15.6e}'*13).format(A, E, H, A*A, A*E, A*H, E*E, E*H, H*H, a0, a1, a2, a3) + '\n')
        if (flush):
            file.flush()

    muon   = open('{}/{}_coeff_data_muon.txt'.format(savedir, model.name),   'w')
    photon = open('{}/{}_coeff_data_photon.txt'.format(savedir, model.name), 'w')

    write_header(muon)
    write_header(photon)

    total_fits = len(channels) * len(primaries) * len(energies) * len(altitudes)
    fit_count  = 0
    for channel in channels:
        for primary in primaries:
            for energy in energies:
                for altitude in altitudes:
                    fit_count += 1

                    name = 'THIN_{}_{}_density_{}_{}'.format(thin, N, channel, altitude)
                    dpmjet = f.Get('{}/{}/{}/{}'.format(energy, 'DPMJET', primary, name))
                    qgsjet = f.Get('{}/{}/{}/{}'.format(energy, 'QGSJET', primary, name))
                    qgsii  = f.Get('{}/{}/{}/{}'.format(energy, 'QGSII',  primary, name))
                    sibyll = f.Get('{}/{}/{}/{}'.format(energy, 'SIBYLL', primary, name))
                    venus  = f.Get('{}/{}/{}/{}'.format(energy, 'VENUS',  primary, name))

                    hists = []
                    wrong_min = []
                    for h in [dpmjet, qgsjet, qgsii, sibyll, venus]:
                        if (h != None):
                            if (h.GetXaxis().GetXmin() > __XMIN__):
                                wrong_min.append(h)
                            else:
                                hists.append(h)

                    if (len(wrong_min) > len(hists)):
                        with open('{}/{}_info.txt'.format(savedir, model.name), 'w') as g:
                            info = ' [INFO]: switching to wrong min for {}'.format(name)
                            g.write(info + '\n')
                            print(info)
                        hists = wrong_min

                    if (len(hists) > 0):
                        hsum = hists[0]
                        for h in hists[1:]:
                            hsum.Add(h)
                        hsum.Scale(1./float(len(hists)))

                        min_, max_ = plot.min_max(hsum)
                        if (min_ is not None and max_ is not None):
                            hsum.SetMaximum(max_ * 100.)
                            hsum.SetMinimum(min_ / 100.)
                        hsum.SetMarkerColor(R.kBlack)
                        hsum.SetMarkerStyle(R.kFullCircle)

                        xaxis = hsum.GetXaxis()
                        for stop_bin in range(xaxis.GetNbins(), 0, -1):
                            val = hsum.GetBinContent(stop_bin)
                            if (val > 0.):
                                break
                        limit = xaxis.GetBinCenter(stop_bin)

                        fitfunc = model.fitfunc
                        f1 = R.TF1('fitfunc', fitfunc, __XMIN__, limit)

                        a0_val = model.a0_val
                        a0_low = model.a0_low
                        a0_hi  = model.a0_hi

                        a1_val = model.a1_val
                        a1_low = model.a1_low
                        a1_hi  = model.a1_hi

                        a2_val = model.a2_val
                        a2_low = model.a2_low
                        a2_hi  = model.a2_hi

                        a3_val = model.a3_val
                        a3_low = model.a3_low
                        a3_hi  = model.a3_hi

                        f1.SetParameters(a0_val, a1_val, a2_val, a3_val)
                        f1.SetParLimits(0, a0_low, a0_hi)
                        f1.SetParLimits(1, a1_low, a1_hi)
                        f1.SetParLimits(2, a2_low, a2_hi)
                        f1.SetParLimits(3, a3_low, a3_hi)

                        f1.SetLineColor(R.kBlue)
                        f1.SetLineStyle(1)
                        f1.SetLineWidth(1)

                        legend = plot.make1Dlegend()
                        legend.AddEntry(hsum, '{}: n={}, thin={}, model ave'.format(plot.__symbols__[channel], N, thin), 'p')
                        legend.AddEntry(f1, model.equation, 'l')

                        def rolldice(f1):
                            def rolldie(mu, low, hi):
                                sig = (hi - low) / 10.
                                while (True):
                                    val = np.random.normal(mu, sig)
                                    if (low < val and val < hi):
                                        return val

                            val0 = rolldie(a0_val, a0_low, a0_hi)
                            val1 = rolldie(a1_val, a1_low, a1_hi)
                            val2 = rolldie(a2_val, a2_low, a2_hi)
                            val3 = rolldie(a3_val, a3_low, a3_hi)
                            f1.SetParameters(val0, val1, val2, val3)

                        def query(f1):
                            while (True):
                                print()
                                print('fitting {} of {} for: {}_{}_{}_{}_{}'.format(fit_count, total_fits, model.name, channel, primary, energy, altitude))
                                action = input('>>> (r)etry, (s)ave, (a)bandon: ')
                                if (action == 'r'):
                                    return 'r'
                                if (action == 's'):
                                    canvas.SaveAs('{}/{}_{}_{}_{}_{}.png'.format(savedir, model.name, channel, primary, energy, altitude))
                                    fit0 = f1.GetParameter(0)
                                    fit1 = f1.GetParameter(1)
                                    fit2 = f1.GetParameter(2)
                                    fit3 = f1.GetParameter(3)
                                    if (channel == 'muon'):
                                        write_data(muon, primary, energy, altitude, fit0, fit1, fit2, fit3) 
                                    if (channel == 'photon'):
                                        write_data(photon, primary, energy, altitude, fit0, fit1, fit2, fit3)
                                    return 's'
                                if (action == 'a'):
                                    return 'a'
                                print(' input invalid, try again ..')

                        trys = 0
                        while (True):
                            fitstatus = int(hsum.Fit('fitfunc', 'RE Q'))
                            if (fitstatus != 0 or f1.GetChisquare() > 1000 and trys < 100):
                                trys += 1
                                rolldice(f1)
                                continue

                            canvas.cd()
                            canvas.Clear()
                            hsum.GetXaxis().SetRangeUser(__XMIN__,1e3)
                            hsum.Draw('hist p')
                            f1.Draw('same')
                            legend.Draw()
                            canvas.Update()
                            
                            action = query(f1)
                            if (action == 'r'):
                                trys = 0
                                rolldice(f1)
                                continue
                            if (action == 's' or action == 'a'):
                                break
    muon.close()
    photon.close()


