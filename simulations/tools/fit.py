
import os
import ROOT as R
from . import plot


def density(infile, savedir='fit_density', batchmode=True):

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

    a0 = open('{}/a0.txt'.format(savedir), 'w')
    a1 = open('{}/a1.txt'.format(savedir), 'w')
    a2 = open('{}/a2.txt'.format(savedir), 'w')
    a3 = open('{}/a3.txt'.format(savedir), 'w')

    def write_header(file):
        file.write('#channel\tprimary\tenergy\taltitude\tvalue\terror\tchi^2\n#\n')

    def write_data(file, channel, primary, energy, altitude, value, error, chi2):
        file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(channel, primary, energy, altitude, value, error, chi2))

    write_header(a0)
    write_header(a1)
    write_header(a2)
    write_header(a3)

    f = R.TFile(infile)
    canvas = plot.make1Dcanvas('fitcanvas')
    canvas.SetLogy()
    canvas.SetLogx()

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
                        f1.SetParameters(1e3, .5, 10., .4)
                        f1.SetParLimits(1, .1, 1.5)
                        f1.SetParLimits(2, 1., 20.)
                        f1.SetParLimits(3, 0., 1.)
                        f1.SetLineColor(R.kBlue)
                        f1.SetLineStyle(1)
                        f1.SetLineWidth(1)

                        legend = plot.make1Dlegend()
                        legend.AddEntry(hsum, '{}: n={}, thin={}, model ave'.format(plot.__symbols__[channel], N, thin), 'p')
                        legend.AddEntry(f1, 'a_{0} x^{-a_{1}} exp(-a_{2} x^{a_{3}})', 'l')

                        canvas.cd()
                        canvas.Clear()
                        hsum.Draw('hist p')
                        result = hsum.Fit('fitfunc', 'RS')
                        f1.Draw('same')
                        legend.Draw()
                        canvas.Update()
                        canvas.SaveAs('{}/{}_{}_{}_{}.png'.format(savedir, channel, primary, energy, altitude))

                        chi2 = result.Chi2()
                        fit0 = result.Parameter(0)
                        fit1 = result.Parameter(1)
                        fit2 = result.Parameter(2)
                        fit3 = result.Parameter(3)
                        err0 = result.ParError(0)
                        err1 = result.ParError(1)
                        err2 = result.ParError(2)
                        err3 = result.ParError(3)

                        write_data(a0, channel, primary, energy, altitude, fit0, err0, chi2)
                        write_data(a1, channel, primary, energy, altitude, fit1, err1, chi2)
                        write_data(a2, channel, primary, energy, altitude, fit2, err2, chi2)
                        write_data(a3, channel, primary, energy, altitude, fit3, err3, chi2)

    a0.close()
    a1.close()
    a2.close()
    a3.close()


