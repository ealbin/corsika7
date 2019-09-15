
import os
import numpy as np
import pandas as pd
import ROOT as R
from . import plot

import matplotlib.pyplot as plt


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
        return np.log10(E) / 4. - 1.

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


    plist = []
    for iteration in range(5):

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

                                a2_val = 10.
                                a2_low = 0.
                                a2_hi  = 20.

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

                                    if (i == 0):
                                        if (p < 0. or p > 5.):
                                            p = .5
                                        else:
                                            pass
                                        f1.SetParLimits(i, 0., 5.)
                                        
                                    if (i == 3):
                                        if (p < 0. or p > 5.):
                                            p = .5
                                        f1.SetParLimits(i, 0., 5.)

                                    if (i == 1):
                                        plist.append(p)
                                        if (p < 0. or p > 10.):
                                            n_ignored += 1
                                            p = 5.
                                        else:
                                            n_saved += 1
                                        f1.SetParLimits(i, 0., 10.)

                                    if (i == 2):
                                        if (p < 0. or p > 20.):
                                            p = 10.
                                        else:
                                            pass
                                        f1.SetParLimits(i, 0., 20.)

                                    f1.SetParameter(i, p)

                                    if (iteration == 2):
                                        f1.SetParLimits(i, .5*p, 2.*p)
                                    elif (iteration == 3):
                                        f1.SetParLimits(i, .6*p, 1.4*p)
                                    elif (iteration == 4):
                                        f1.SetParLimits(i, .7*p, 1.3*p)
                                    else:
                                        f1.SetParLimits(i, .8*p, 1.2*p)

                                    #if (iteration == 1):
                                    #    if (i != 1):
                                    #        f1.SetParameter(i, p)
                                    #        f1.SetParLimits(i, .5*p, 1.5*p) 
                                    #    else:
                                    #        f1.SetParameter(1, 5.)
                                    #        f1.SetParLimits(1, 0., 10.)
                                    #elif (iteration == 2):
                                    #    f1.SetParameter(i, p)
                                    #    if (i != 1):
                                    #        f1.SetParLimits(i, .75*p, 1.25*p)
                                    #    else:
                                    #        f1.SetParLimits(1, .5*p, 1.5*p)
                                    #else: # (iteration == 3):
                                    #    f1.SetParameter(i, p)
                                    #    f1.SetParLimits(i, .8*p, 1.2*p)

                            
                            f1.SetLineColor(R.kBlue)
                            f1.SetLineStyle(1)
                            f1.SetLineWidth(1)

                            legend = plot.make1Dlegend()
                            legend.AddEntry(hsum, '{}: n={}, thin={}, model ave'.format(plot.__symbols__[channel], N, thin), 'p')
                            legend.AddEntry(f1, 'x^{-a_{0}} exp(a_{1} - a_{2} x^{a_{3}})', 'l')

                            fitstatus = int(hsum.Fit('fitfunc', 'R Q'))

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

                            if (fit0 < 0. or fit1 < 0. or fit2 < 0. or fit3 < 0.
                                or 5. < fit0 or 10. < fit1 or 10. < fit2 or 5. < fit3):
                                pass
                            else:
                                write_data(a0, channel, primary, energy, altitude, fit0, err0, chi2, fitstatus)
                                write_data(a1, channel, primary, energy, altitude, fit1, err1, chi2, fitstatus)
                                write_data(a2, channel, primary, energy, altitude, fit2, err2, chi2, fitstatus)
                                write_data(a3, channel, primary, energy, altitude, fit3, err3, chi2, fitstatus)

        a0.close()
        a1.close()
        a2.close()
        a3.close()

        print('   N converged: {:<10d} / failed: {:<10d} = {:5.2f}% success'.format(n_converged, n_failed, 100. * n_converged / float(n_converged + n_failed)))
        if (iteration > 0):
            print('   N saved: {:<10d} / ignored: {:<10d} = {:5.2f}% success'.format(n_saved, n_ignored, 100. * n_saved / float(n_saved + n_ignored)))

    coefficients(iteration)
    plt.hist(plist, bins=30)
    plt.yscale('log')
    return
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



def density_check(infile, coeff_file, savedir='fit_density', batchmode=True):

    if (batchmode == True):
        R.gROOT.SetBatch(R.kTRUE)

    if (not os.path.isdir(savedir)):
        os.makedirs(savedir, exist_ok=True)
    
    df = pd.read_csv(coeff_file, delim_whitespace=True)
    channels = ['photon'] # df.Channel.unique()
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

    n_converged = 0
    n_failed = 0
    chi2 = []

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

                        fitfunc = '[0] * x^(-[1]) * exp(-[2] * x^[3])'
                        f1 = R.TF1('fitfunc', fitfunc, 1e-3, limit)

                        for i in range(4):
                            A = MassNumber[primary]
                            E = np.log10(float(energy))
                            H = alt_km[int(altitude)]
                            
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
