
import os
import numpy as np
import pandas as pd
import ROOT as R
from . import plot

import matplotlib.pyplot as plt

def A_func(A):
    return np.log(A)

def E_func(E):
    return np.log10(E)

def H_func(H):
    return H


def quick_check(indir, channel='muon'):
    df = pd.read_csv(os.path.join(indir, 'coeff_data_{}.txt'.format(channel)), delim_whitespace=True)

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

    plt.gcf().suptitle(channel) 


def density_check(infile, indir, savedir='fit_check', batchmode=True):

    if (batchmode == True):
        R.gROOT.SetBatch(R.kTRUE)

    if (not os.path.isdir(savedir)):
        os.makedirs(savedir, exist_ok=True)
    
    alt_km = [0., 100., 50., 20., 10., 5., 2., 1.4, 1.0, .5]
    MassNumber = {'He':4, 'O':16, 'Fe':56}

    thin      = '1E-6'
    N         = '100'
    primaries = ['He', 'O', 'Fe']
    energies  = ['1e15', '1e16', '1e17', '1e18', '1e19', '1e20', '1e21']
    altitudes = ['3', '4', '5', '6', '7', '8', '9', '0']
    models    = ['DPMJET', 'QGSJET', 'QGSII', 'SIBYLL', 'VENUS']
    
    f = R.TFile(infile)
    canvas = plot.make1Dcanvas('check_canvas')
    canvas.SetLogy()
    canvas.SetLogx()
    
    for channel in ['muon', 'photon']:
        df = pd.read_csv(os.path.join(indir, 'coeff_data_{}.txt'.format(channel)), delim_whitespace=True)
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
                    for h in [dpmjet, qgsjet, qgsii, sibyll, venus]:
                        if (h != None):
                            hists.append(h)

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

                        fitfunc = 'x^(-[0]) * exp(-[1] - [2] * x^[3])'
                        f1 = R.TF1('fitfunc', fitfunc, 1e-3, limit)

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
                        legend.AddEntry(f1, 'x^{-a_{0}} exp(-a_{1} - a_{2} x^{a_{3}})', 'l')

                        canvas.cd()
                        canvas.Clear()
                        hsum.Draw('hist p')
                        f1.Draw('same')
                        legend.Draw()
                        canvas.Update()
                        canvas.SaveAs('{}/check_{}_{}_{}_{}.png'.format(savedir, channel, primary, energy, altitude))


def supervised_density(infile, savedir='fit_density'):

    if (not os.path.isdir(savedir)):
        os.makedirs(savedir, exist_ok=True)
    
    thin      = '1E-6'
    N         = '100'
    channels  = ['muon', 'photon']
    primaries = ['He', 'O', 'Fe']
    energies  = ['1e15', '1e16', '1e17', '1e18', '1e19', '1e20', '1e21']
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

    muon   = open('{}/coeff_data_muon.txt'.format(savedir),   'w')
    photon = open('{}/coeff_data_photon.txt'.format(savedir), 'w')

    write_header(muon)
    write_header(photon)

    for channel in channels:
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
                    for h in [dpmjet, qgsjet, qgsii, sibyll, venus]:
                        if (h != None):
                            hists.append(h)

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

                        fitfunc = 'x^(-[0]) * exp(-[1] - [2] * x^[3])'
                        f1 = R.TF1('fitfunc', fitfunc, 1e-3, limit)

                        a0_val = .5
                        a0_low = -10.
                        a0_hi  = 20.

                        a1_val = 5.
                        a1_low = -30.
                        a1_hi  = 30.

                        a2_val = 5.
                        a2_low = -10.
                        a2_hi  = 35.

                        a3_val = .5
                        a3_low = -5.
                        a3_hi  = 5.

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
                        legend.AddEntry(f1, 'x^{-a_{0}} exp(-a_{1} - a_{2} x^{a_{3}})', 'l')

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
                                action = input('(r)etry, (s)ave, (a)bandon: ')
                                if (action == 'r'):
                                    return 'r'
                                if (action == 's'):
                                    canvas.SaveAs('{}/{}_{}_{}_{}.png'.format(savedir, channel, primary, energy, altitude))
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
                            if (fitstatus != 0 and trys < 10):
                                trys += 1
                                rolldice(f1)
                                continue

                            canvas.cd()
                            canvas.Clear()
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


#################################################
# deleteable
###########


def density(infile, savedir='fit_density', batchmode=True, chi2lim=2000):

    if (batchmode == True):
        R.gROOT.SetBatch(R.kTRUE)

    if (not os.path.isdir(savedir)):
        os.makedirs(savedir, exist_ok=True)
    
    thin      = '1E-6'
    N         = '100'
    channels  = ['muon', 'photon']
    primaries = ['He', 'O', 'Fe']
    energies  = ['1e15', '1e16', '1e17', '1e18', '1e19', '1e20', '1e21']
    altitudes = ['3', '4', '5', '6', '7', '8', '9', '0'] 
    models    = ['DPMJET', 'QGSJET', 'QGSII', 'SIBYLL', 'VENUS']

    alt_km = [0., 100., 50., 20., 10., 5., 2., 1.4, 1.0, .5]
    MassNumber = {'He':4, 'O':16, 'Fe':56}

    f = R.TFile(infile)
    canvas = plot.make1Dcanvas('fitcanvas')
    canvas.SetLogy()
    canvas.SetLogx()
   
    def A_func(A):
        return np.log(A)

    def E_func(E):
        return np.log10(E) / 5.

    def H_func(H):
        return H / 5.

    def coefficients(iteration):
        with open(os.path.join(savedir, 'coefficients_{}.txt'.format(iteration)), 'w') as fcoeff:
            fcoeff.write('{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}\n'
                  .format('Channel', 'a', 'Const', 'A', 'E', 'H', 'AA', 'AE', 'AH', 'EE', 'EH', 'HH'))

            for i in range(4):

                inf = os.path.join(savedir, 'a{}_{}.txt'.format(i, iteration)) 
                df  = pd.read_csv(inf, delim_whitespace=True)
                df  = df[ (df.Fitstatus == 0) & (df.Chi2 < chi2lim) ]

                for channel in channels:
                
                    outfile = os.path.join(savedir, 'a{}_{}_matrix_{}.txt'.format(i, channel, iteration))
                    with open(outfile, 'w') as fmatrix:
                        fmatrix.write('{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}\n'
                               .format('Constant', 'LogMassNum_A', 'LogEnergy_E', 'Altitude_H', 'AA', 'AE', 'AH', 'EE', 'EH', 'HH', 'Value'))

                        for index, row in df[df.Channel == channel].iterrows():
                            A = A_func(MassNumber[ row['Primary'] ])
                            E = E_func(row['Energy'])
                            H = H_func(alt_km[ row['Altitude'] ])
                            V = row['Value']

                            fmatrix.write('{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}\n'
                                   .format(1., A, E, H, A*A, A*E, A*H, E*E, E*H, H*H, V))
                    
                    of = pd.read_csv(outfile, delim_whitespace=True)
                    A = of.loc[:, of.columns != 'Value'].to_numpy()
                    B = np.linalg.pinv(A)
                    
                    if (not np.allclose(A, np.dot(A, np.dot(B, A)))):
                        print('Inverse failed for: {}'.format(outfile))
                    
                    V = of.Value.to_numpy()
                    C = np.dot(B, V)
                    fcoeff.write('{:15s}{:15s}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}\n'
                          .format(channel, 'a{}'.format(i), C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9]))


    for iteration in range(7):

        print('\n\n>>>> Starting iteration: {}'.format(iteration))
        n_converged = 0
        n_failed    = 0
        n_saved     = 0
        n_ignored   = 0
        
        if (iteration > 0):
            coefficients(iteration - 1) 

        a0 = open('{}/a0_{}.txt'.format(savedir, iteration), 'w')
        a1 = open('{}/a1_{}.txt'.format(savedir, iteration), 'w')
        a2 = open('{}/a2_{}.txt'.format(savedir, iteration), 'w')
        a3 = open('{}/a3_{}.txt'.format(savedir, iteration), 'w')

        def write_header(file):
            file.write('{:10s}{:10s}{:10s}{:10s}{:15s}{:15s}{:15s}{:10s}\n'.format('Channel', 'Primary', 'Energy', 'Altitude', 'Value', 'Error', 'Chi2', 'Fitstatus'))

        def write_data(file, channel, primary, energy, altitude, value, error, chi2, fitstatus):
            file.write('{:10s}{:10s}{:10s}{:10s}{:<15.4e}{:<15.4e}{:<15.4e}{:<10d}\n'.format(channel, primary, energy, altitude, value, error, chi2, fitstatus))

        write_header(a0)
        write_header(a1)
        write_header(a2)
        write_header(a3)

        for channel in channels:
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
                        for h in [dpmjet, qgsjet, qgsii, sibyll, venus]:
                            if (h != None):
                                hists.append(h)

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

                            fitfunc = 'x^(-[0]) * exp([1] - [2] * x^[3])'
                            f1 = R.TF1('fitfunc', fitfunc, 1e-3, limit)
                            
                            if (iteration == 0):
                                a0_val = .5
                                a0_low = 0.
                                a0_hi  = 5.

                                a1_val = 5.
                                a1_low = 0.
                                a1_hi  = 10.

                                a2_val = 5.
                                a2_low = 0.
                                a2_hi  = 10.

                                a3_val = .5
                                a3_low = 0.
                                a3_hi  = 5.

                                f1.SetParameters(a0_val, a1_val, a2_val, a3_val)
                                f1.SetParLimits(0, a0_low, a0_hi)
                                f1.SetParLimits(1, a1_low, a1_hi)
                                f1.SetParLimits(2, a2_low, a2_hi)
                                f1.SetParLimits(3, a3_low, a3_hi)

                            else:
                                df = pd.read_csv(os.path.join(savedir, 'coefficients_{}.txt'.format(iteration - 1)), delim_whitespace=True) 
                                for i in range(4):
                                    A = A_func(MassNumber[primary])
                                    E = E_func(float(energy))
                                    H = H_func(alt_km[int(altitude)])
                                    
                                    _ = df[ (df.Channel == channel) & (df.a == 'a{}'.format(i)) ].iloc[0]                   
                                    p =  _.Const + _.A*A + _.E*E + _.H*H + _.AA*(A*A) + _.AE*(A*E) + _.AH*(A*H) + _.EE*(E*E) + _.EH*(E*H) + _.HH*(H*H) 

                                    if (p < 0.):
                                        if (i == 0 or i == 3):
                                            p = .5
                                        else:
                                            p = 5.

                                    if (iteration == 1):
                                        if (i != 1):
                                            f1.SetParameter(i, p)
                                            f1.SetParLimits(i, .5*p, 1.5*p) 
                                        else:
                                            f1.SetParameter(1, 5.)
                                            f1.SetParLimits(1, 0., 10.)
                                    elif (iteration == 2):
                                        f1.SetParameter(i, p)
                                        if (i != 1):
                                            f1.SetParLimits(i, .75*p, 1.25*p)
                                        else:
                                            f1.SetParLimits(1, .5*p, 1.5*p)
                                    else: # (iteration == 3):
                                        f1.SetParameter(i, p)
                                        f1.SetParLimits(i, .8*p, 1.2*p)

                            
                            f1.SetLineColor(R.kBlue)
                            f1.SetLineStyle(1)
                            f1.SetLineWidth(1)

                            legend = plot.make1Dlegend()
                            legend.AddEntry(hsum, '{}: n={}, thin={}, model ave'.format(plot.__symbols__[channel], N, thin), 'p')
                            legend.AddEntry(f1, 'x^{-a_{0}} exp(a_{1} - a_{2} x^{a_{3}})', 'l')

                            fitstatus = int(hsum.Fit('fitfunc', 'RE Q'))

                            if (fitstatus == 0):
                                n_converged += 1
                            else:
                                n_failed += 1

                            chi2 = f1.GetChisquare()
                            fit0 = f1.GetParameter(0)
                            fit1 = f1.GetParameter(1)
                            fit2 = f1.GetParameter(2)
                            fit3 = f1.GetParameter(3)
                            err0 = f1.GetParError(0)
                            err1 = f1.GetParError(1)
                            err2 = f1.GetParError(2)
                            err3 = f1.GetParError(3)

                            if (fit0 < 0. or fit1 < 0.or fit2 < 0. or fit3 < 0.):
                                n_ignored += 1
                            else:
                                n_saved += 1
                                write_data(a0, channel, primary, energy, altitude, fit0, err0, chi2, fitstatus)
                                write_data(a1, channel, primary, energy, altitude, fit1, err1, chi2, fitstatus)
                                write_data(a2, channel, primary, energy, altitude, fit2, err2, chi2, fitstatus)
                                write_data(a3, channel, primary, energy, altitude, fit3, err3, chi2, fitstatus)

        a0.close()
        a1.close()
        a2.close()
        a3.close()

        print('   N converged: {:<10d} / failed: {:<10d} = {:5.2f}% success'.format(n_converged, n_failed, 100. * n_converged / float(n_converged + n_failed)))
        print('   N saved: {:<10d} / ignored: {:<10d} = {:5.2f}% success'.format(n_saved, n_ignored, 100. * n_saved / float(n_saved + n_ignored)))

    coefficients(iteration)
   
    for channel in channels:
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
                    for h in [dpmjet, qgsjet, qgsii, sibyll, venus]:
                        if (h != None):
                            hists.append(h)

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

                        fitfunc = 'x^(-[0]) * exp([1] - [2] * x^[3])'
                        f1 = R.TF1('fitfunc', fitfunc, 1e-3, limit)
                        
                        df = pd.read_csv(os.path.join(savedir, 'coefficients_{}.txt'.format(iteration)), delim_whitespace=True) 
                        for i in range(4):
                            A = A_func(MassNumber[primary])
                            E = E_func(float(energy))
                            H = H_func(alt_km[int(altitude)])
                            
                            _ = df[ (df.Channel == channel) & (df.a == 'a{}'.format(i)) ].iloc[0]                   
                            p =  _.Const + _.A*A + _.E*E + _.H*H + _.AA*(A*A) + _.AE*(A*E) + _.AH*(A*H) + _.EE*(E*E) + _.EH*(E*H) + _.HH*(H*H) 
                            f1.SetParameter(i, p)
                        
                        f1.SetLineColor(R.kBlue)
                        f1.SetLineStyle(1)
                        f1.SetLineWidth(1)

                        legend = plot.make1Dlegend()
                        legend.AddEntry(hsum, '{}: n={}, thin={}, model ave'.format(plot.__symbols__[channel], N, thin), 'p')
                        legend.AddEntry(f1, 'x^{-a_{0}} exp(a_{1} - a_{2} x^{a_{3}})', 'l')

                        canvas.cd()
                        canvas.Clear()
                        hsum.Draw('hist p')
                        f1.Draw('same')
                        legend.Draw()
                        canvas.Update()
                        canvas.SaveAs('{}/{}_{}_{}_{}.png'.format(savedir, channel, primary, energy, altitude))

    print('finished!')



def density_coefficients(indir, chi2lim=2000.):

    alt_km = [0., 100., 50., 20., 10., 5., 2., 1.4, 1.0, .5]
    MassNumber = {'He':4, 'O':16, 'Fe':56}

    with open(os.path.join(indir, 'coefficients.txt'), 'w') as fcoeff:
        fcoeff.write('{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}\n'
              .format('Channel', 'a', 'Const', 'A', 'E', 'H', 'AA', 'AE', 'AH', 'EE', 'EH', 'HH'))

        for i in range(4):

            infile = os.path.join(indir, 'a{}.txt'.format(i)) 
            df = pd.read_csv(infile, delim_whitespace=True)
            df = df[ (df.fitstatus == 0) & (df.chi2 < chi2lim) ]

            for channel in ['photon', 'muon']:
            
                outfile = os.path.join(indir, 'a{}_{}_matrix.txt'.format(i, channel))
                with open(outfile, 'w') as f:
                    f.write('{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}\n'
                     .format('Constant', 'MassNum_A', 'LogEnergy_E', 'Altitude_H', 'A2', 'AE', 'AH', 'E2', 'EH', 'H2', 'Value'))

                    for index, row in df[df.channel == channel].iterrows():
                        A = MassNumber[ row['primary'] ]
                        E = np.log10(row['energy'])
                        H = alt_km[ row['altitude'] ]
                        V = row['value']

                        f.write('{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}\n'
                         .format(1., A, E, H, A*A, A*E, A*H, E*E, E*H, H*H, V))
                
                of = pd.read_csv(outfile, delim_whitespace=True)
                A = of.loc[:, of.columns != 'Value'].to_numpy()
                B = np.linalg.pinv(A)
                
                if (not np.allclose(A, np.dot(A, np.dot(B, A)))):
                    print('Inverse failed for: {}'.format(outfile))
                
                V = of.Value.to_numpy()
                C = np.dot(B, V)
                fcoeff.write('{:15s}{:15s}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}\n'
                      .format(channel, 'a{}'.format(i), C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9]))



def density_check_old(infile, coeff_file, channel, savedir='fit_check', batchmode=True):

    if (batchmode == True):
        R.gROOT.SetBatch(R.kTRUE)

    if (not os.path.isdir(savedir)):
        os.makedirs(savedir, exist_ok=True)
    
    df = pd.read_csv(coeff_file, delim_whitespace=True)
    df.insert(0, 'C', 1.)

    coeff = pd.DataFrame(columns=df.columns[:10])
    A = df.iloc[:,:10].to_numpy()
    B = np.linalg.pinv(A)
    for i in range(4):
        V = df.loc[:,'a{}'.format(i)].to_numpy()
        coeff.loc[-1] = np.dot(B, V)

    def a(n, A, E, H):
        A = np.log(A)
        E = np.log10(E)
        H = H

    alt_km = [0., 100., 50., 20., 10., 5., 2., 1.4, 1.0, .5]
    MassNumber = {'He':4, 'O':16, 'Fe':56}

    thin      = '1E-6'
    N         = '100'
    primaries = ['He', 'O', 'Fe']
    energies  = ['1e15', '1e16', '1e17', '1e18', '1e19', '1e20', '1e21']
    altitudes = ['3', '4', '5', '6', '7', '8', '9', '0']
    models    = ['DPMJET', 'QGSJET', 'QGSII', 'SIBYLL', 'VENUS']
    
    f = R.TFile(infile)
    canvas = plot.make1Dcanvas('check_canvas')
    canvas.SetLogy()
    canvas.SetLogx()

    chi2 = []
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
                for h in [dpmjet, qgsjet, qgsii, sibyll, venus]:
                    if (h != None):
                        hists.append(h)

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

                    fitfunc = 'x^(-[0]) * exp([1] - [2] * x^[3])'
                    f1 = R.TF1('fitfunc', fitfunc, 1e-3, limit)

                    for i in range(4):
                        A = A_func(MassNumber[primary])
                        E = E_func(float(energy))
                        H = H_func(alt_km[int(altitude)])        
                        
                        _ = df[ (df.Channel == channel) & (df.a == 'a{}'.format(i)) ].iloc[0]                   
                        p =  _.Const + _.A*A + _.E*E + _.H*H + _.AA*(A*A) + _.AE*(A*E) + _.AH*(A*H) + _.EE*(E*E) + _.EH*(E*H) + _.HH*(H*H) 

                        f1.SetParameter(i, p)
                        if (i > 0):
                            f1.SetParLimits(i, .9*p, 1.1*p) 

                    f1.SetLineColor(R.kBlue)
                    f1.SetLineStyle(1)
                    f1.SetLineWidth(1)

                    legend = plot.make1Dlegend()
                    legend.AddEntry(hsum, '{}: n={}, thin={}, model ave'.format(plot.__symbols__[channel], N, thin), 'p')
                    legend.AddEntry(f1, 'a_{0} x^{-a_{1}} exp(-a_{2} x^{a_{3}})', 'l')

                    canvas.cd()
                    canvas.Clear()
                    hsum.Draw('hist p')
                    if (int(hsum.Fit('fitfunc', 'R')) == 0):
                        n_converged += 1
                        chi2.append(f1.GetChisquare())
                    else:
                        n_failed += 1
                    f1.Draw('same')
                    legend.Draw()
                    canvas.Update()
                    canvas.SaveAs('{}/check_{}_{}_{}_{}.png'.format(savedir, channel, primary, energy, altitude))
    print('\n\nN converged: {:d}, N failed: {:d}, {:.2f}% success'.format(n_converged, n_failed, 100. * n_converged / float(n_converged + n_failed)))
    plt.hist(chi2, bins=20)


def vs_altitude(indir, channel='photon', primary='O', energy=1e18):

    alt_km = [0., 100., 50., 20., 10., 5., 2., 1.4, 1.0, .5]

    c = []
    g = []
    for i in range(4):
        c.append(plot.make1Dcanvas('a{}'.format(i)))

        infile = os.path.join(indir, 'a{}.txt'.format(i))
        df = pd.read_csv(infile, delim_whitespace=True)

        df = df[ df.fitstatus == 0.0 ].reset_index()
        df.drop('fitstatus', axis=1)

        x  = []
        y  = []
        ex = []
        ey = []
        for altitude in df.altitude.unique():
            select = df[ (df.altitude == altitude) 
                       & (df.channel  == channel)
                       & (df.energy   == energy)
                       & (df.primary  == primary)]
            for v,e in zip(select.value, select.error):
                x.append(alt_km[altitude])
                y.append(v)
                ex.append(0.)
                ey.append(e)
            
        x = np.asarray(x)
        y = np.asarray(y)
        ex = np.asarray(ex)
        ey = np.asarray(ey)

        title = os.path.basename('a_{}'.format(i))
        g.append(R.TGraphErrors(x.size, x, y, ex, ey))
        g[-1].SetTitle(title)
        g[-1].GetXaxis().SetTitle('Altitude [km]')
        g[-1].SetMarkerStyle(R.kFullDotSmall)
        c[-1].cd()
        g[-1].Draw('ap')
        c[-1].Update()

    return c, g
