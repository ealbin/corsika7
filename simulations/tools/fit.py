
import os
import numpy as np
import pandas as pd
from scipy.special import gamma
import ROOT as R
from . import plot
from . import settings

import matplotlib as mpl
import matplotlib.pyplot as plt


__XMIN__ = 1e-3
__XMAX__ = 1e4

def generate(maxA=None, maxE=None, maxH=None):
    models={'linear':[0, 1], 
            'even':[0,2,4], 'odd':[1,3,5], 
            'pol2':[0,1,2], 'pol3':[0,1,2,3], 'pol4':[0,1,2,3,4], 'pol5':[0,1,2,3,4,5]}
    
    for key in models.keys():
        for model in ['Alt', 'NKG']:
            for channel in ['muon', 'photon']:
                Coefficients.make('fit_density/{}_coeff_data_{}.txt'.format(model, channel), 'fit_models/{}_{}_{}_model.txt'.format(model, channel, key), orders=models[key], maxA=maxA, maxE=maxE, maxH=maxH)  

class Coefficients:

    orders = []
    orders.append(['1'])
    orders.append(['A', 'E', 'H'])
    orders.append(['AA', 'AE', 'AH', 'EE', 'EH', 'HH'])
    orders.append(['AAA', 'AAE', 'AAH', 'AEE', 'AEH', 'AHH', 'EEE', 'EEH', 'EHH', 'HHH'])
    orders.append(['AAAA', 'AAAE', 'AAAH', 'AAEE', 'AAEH', 'AAHH', 'AEEE', 'AEEH', 'AEHH', 'AHHH', 'EEEE', 'EEEH', 'EEHH', 'EHHH', 'HHHH'])
    orders.append(['AAAAA', 'AAAAE', 'AAAAH', 'AAAEE', 'AAAEH', 'AAAHH', 'AAEEE', 'AAEEH', 'AAEHH', 'AAHHH', 'AEEEE', 'AEEEH', 'AEEHH', 'AEHHH', 'AHHHH', 'EEEEE', 'EEEEH', 'EEEHH', 'EEHHH', 'EHHHH', 'HHHHH'])

    def make(fit_data, outfile, orders=[1,3], maxA=None, maxE=None, maxH=None):
        df = pd.read_csv(fit_data, delim_whitespace=True)
        orders.sort()

        A = Transform.A(df.A).to_numpy()
        E = Transform.E(df.E).to_numpy()
        H = Transform.H(df.H).to_numpy()
        a0 = df.a0.to_numpy()
        a1 = df.a1.to_numpy()
        a2 = df.a2.to_numpy()
        a3 = df.a3.to_numpy()

        model = ''
        for order in orders:
            for term in Coefficients.orders[order]:
                model += term + ' '
        model.rstrip()
        keys = model.split()
        remove = []
        dic = {}
        for key in keys:
            if (maxA is not None and key.count('A') > maxA):
                remove.append(key)
                continue
            if (maxE is not None and key.count('E') > maxE):
                remove.append(key)
                continue
            if (maxH is not None and key.count('H') > maxH):
                remove.append(key)
                continue

            val = []
            for _ in range(len(df)):
                expr = key
                expr = expr.replace('A', '{}*'.format(A[_]))
                expr = expr.replace('E', '{}*'.format(E[_]))
                expr = expr.replace('H', '{}*'.format(H[_]))
                expr = expr.rstrip('*')
                val.append(eval(expr))
            dic[key] = val
        dic['a0'] = a0
        dic['a1'] = a1
        dic['a2'] = a2
        dic['a3'] = a3

        for _ in remove:
            keys.remove(_)

        raw = pd.DataFrame.from_dict(dic)
        _A = raw.iloc[:,:-4].to_numpy()
        _B = np.linalg.pinv(_A)
        _C = []
        _C.append(np.dot(_B, a0))
        _C.append(np.dot(_B, a1))
        _C.append(np.dot(_B, a2))
        _C.append(np.dot(_B, a3))

        with open(outfile, 'w') as f:
            head = '{:15s}'.format('a')
            for key in keys:
                head += '{:15s}'.format(key)
            f.write(head + '\n')

            for i in range(len(_C)):
                row = '{:<15d}'.format(i)
                for term in _C[i]:
                    row += '{:<15.6e}'.format(term)
                f.write(row + '\n')


    def get(index, model, A, E, H):
        mf = pd.read_csv(model, delim_whitespace=True)
        terms = mf.columns[1:]
        
        mA = Transform.A(A)
        mE = Transform.E(E)
        mH = Transform.H(H)
        mC = mf.iloc[index, 1:].to_numpy()

        mX = []
        for term in terms:
            term = term.replace('A', '{}*'.format(mA))
            term = term.replace('E', '{}*'.format(mE))
            term = term.replace('H', '{}*'.format(mH))
            term = term.rstrip('*')
            mX.append(eval(term))
        mX = np.asarray(mX)

        return np.dot(mX, mC)


    def plot(data, model, yscale=2.):
        df = pd.read_csv(data,  delim_whitespace=True)
        mf = pd.read_csv(model, delim_whitespace=True)
        
        terms = mf.columns[1:]

        e = np.logspace(9, 39, 100)
        A = df.A.unique()
        H = df.H.unique()
        A.sort()
        H.sort()
        H = H[::-1]

        colors = ['m','b','r']

        for h in H:
            dff = df[ df.H == h ]
            plt.figure(figsize=[15,15])
            plt.gcf().suptitle('{:.1f} km, {}'.format(h, model))
            for i in range(4):
                plt.subplot(2, 2, i + 1)
                plt.title('a{}'.format(i))
                plt.xscale('log')
                plt.xlim([e[0], e[-1]])
                a_min = None
                a_max = None
                for j, _ in enumerate(A):
                    E = dff[ dff.A == _ ].E
                    a = dff[ dff.A == _ ].loc[:,'a{}'.format(i)]
                    plt.scatter(E, a, color=colors[j], marker='x', label='A = {}'.format(int(_)))

                    a__min = a.min()
                    a__max = a.max()
                    if (a_min is None or a__min < a_min):
                        a_min = a__min
                    if (a_max is None or a__max > a_max):
                        a_max = a__max

                    mA = Transform.A(_)
                    mE = Transform.E(e)
                    mH = Transform.H(h)
                    mC = mf.iloc[i,1:].to_numpy()
                    ma = []
                    for me in mE:
                        mX = []
                        for term in terms:
                            term = term.replace('A', '{}*'.format(mA))
                            term = term.replace('E', '{}*'.format(me))
                            term = term.replace('H', '{}*'.format(mH))
                            term = term.rstrip('*')
                            mX.append(eval(term))
                        ma.append(np.dot(mX, mC))
                    plt.plot(e, ma, color=colors[j])

                if (a_min < 0.):
                    a_min = a_min * yscale
                else:
                    a_min = a_min / yscale
                if (a_max < 0.):
                    a_max = a_max / yscale
                else:
                    a_max = a_max * yscale
                plt.ylim([a_min, a_max])
                plt.legend()


class Transform:
    def A(A):
        return np.log(A + 1.)
    def Ainv(a):
        return np.exp(a) - 1.

    def E(E):
        return np.log10(E / 1e18)
    def Einv(e):
        return np.power(10., e) * 1e21

    def H(H):
        return np.log(H + 20.)
    def Hinv(h):
        return np.exp(h) - 20.


class NKG:
    name      = 'NKG'
    fitexpr   = 'exp([0]) * (1000.*x)^(-[1]) * (1. + (1000.*x)/exp([2]))^(-[3])'
    printexpr = 'exp(a_{0}) x^{-a_{1}} (1 + x / exp(a_{2}))^{-a_{3}}'
    a_set     = [  1.,  1.,  1.,  1.]
    a_min     = [-30., -5., -5.,  0.]
    a_max     = [ 30., 10., 10., 10.]

class Alt:
    name      = 'Alt'
    fitexpr   = 'exp([0]) * (1000.*x)^(-[1]) * exp(-((1000.*x)^[3]) / exp([2]))'
    printexpr = 'exp(a_{0}) x^{-a_{1}} exp(-x^{a_{3}} / exp(a_{2}))' 
    a_set     = [  1.,  1.,  1.,  1.]
    a_min     = [-30., -5., -5.,  0.]
    a_max     = [ 30., 10., 10., 10.]  


class Density:


    def __init__(self, indir, nkg_models, alt_models):
        self.indir = indir

        self.canvas = plot.make1Dcanvas('canvas')
        self.canvas.SetLogx()
        self.canvas.SetLogy()

        self.h_detection = R.TH1D('detection', ';Distance from Shower Core [km];Detected Particles', 100, 0., 1.)
        self.h_detection.SetStats(False)
        self.h_detection.GetXaxis().SetTitleOffset(settings.TitleOffset['x'])
        self.h_detection.GetYaxis().SetTitleOffset(settings.TitleOffset['y'])
        self.h_detection.SetMarkerStyle(R.kFullCircle)
        self.h_detection.SetMarkerColor(R.kBlack)
        self.h_detection.SetLineColor(R.kBlack)
        self.h_detection.SetLineWidth(1)
        self.h_detection.SetLineStyle(1)
       
        self.nkg_models = nkg_models
        self.alt_models = alt_models
        self.nkg = R.TF1('NKG_muon', NKG.fitexpr, __XMIN__, __XMAX__)
        self.alt = R.TF1('Alt_muon', Alt.fitexpr, __XMIN__, __XMAX__)

    def shower(self, mass_number, energy_eV, altitude_km):
        self.A = mass_number
        self.E = energy_eV
        self.H = altitude_km

    def plot_channel(self, channel):
        for i in range(4):
            self.nkg.SetParameter(i, Coefficients.get(i, self.nkg_models[channel], self.A, self.E, self.H))
            self.alt.SetParameter(i, Coefficients.get(i, self.alt_models[channel], self.A, self.E, self.H))
        
        legend = plot.make1Dlegend()
        self.nkg.SetLineColor(R.kRed)
        self.alt.SetLineColor(R.kBlue)
        legend.AddEntry(self.nkg, 'NKG', 'l')
        legend.AddEntry(self.alt, 'Alt', 'l')

        self.canvas.Clear()
        self.nkg.Draw()
        self.alt.Draw('same')
        legend.Draw()
        self.canvas.Update()


    def plot_detection(self, phone_density_km2, phone_area_cm2=0.2, photon_eff=1e-4, muon_eff=.5, ymin=1e-10, ymax=1e10):
        nbins = self.h_detection.GetXaxis().GetNbins()
        bins = np.logspace(np.log10(__XMIN__), np.log10(__XMAX__), nbins)
        self.h_detection.GetXaxis().Set(nbins - 1, bins)
        for i, bin in enumerate(bins):
            if (i == nbins - 1):
                break
            
            center_km    = np.sqrt(bin * bins[i+1])
            bin_area_km2 = np.pi * (bins[i+1]*bins[i+1] - bin*bin)
            sensor_area_cm2 = phone_area_cm2 * (phone_density_km2 * bin_area_km2)

            for j in range(4):
                self.nkg.SetParameter(j, Coefficients.get(j, self.nkg_models['muon'], self.A, self.E, self.H))
                self.alt.SetParameter(j, Coefficients.get(j, self.alt_models['muon'], self.A, self.E, self.H))
            nkg_muons = self.nkg.Eval(center_km)
            alt_muons = self.alt.Eval(center_km)

            for j in range(4):
                self.nkg.SetParameter(j, Coefficients.get(j, self.nkg_models['photon'], self.A, self.E, self.H))
                self.alt.SetParameter(j, Coefficients.get(j, self.alt_models['photon'], self.A, self.E, self.H))
            nkg_photons = self.nkg.Eval(center_km)
            alt_photons = self.alt.Eval(center_km)

            nkg = (nkg_muons * muon_eff) + (nkg_photons * photon_eff)
            alt = (alt_muons * muon_eff) + (alt_photons * photon_eff)
            self.h_detection.SetBinContent(i + 1, sensor_area_cm2 * np.sqrt(nkg * alt) )

        self.canvas.Clear()
        self.h_detection.Draw('hist pl')
        self.h_detection.SetMinimum(ymin)
        self.h_detection.SetMaximum(ymax)
        self.canvas.Update()   
        
        fsum = 0.
        msum = 0.
        xaxis = self.h_detection.GetXaxis()
        for bin in range(nbins + 1):
            val = self.h_detection.GetBinContent(bin + 1)
            msum += val
            if (val > 1.):
                fsum += np.floor(val)
        print('Floor detections = {:.2e}'.format(fsum))
        print('Max detections = {:.2e}'.format(msum))



def likelihood_ratio(infile, indir, model0=NKG, model1=Alt, channel='muon'):
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

                    f0 = R.TF1('model0', model0.fitexpr, __XMIN__, __XMAX__)
                    f1 = R.TF1('model1', model1.fitexpr, __XMIN__, __XMAX__)

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


def density_coeff(indir, model=Alt, channel='muon'):
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


def quick_check(indir, model=Alt, channel='muon'):
    df = pd.read_csv(os.path.join(indir, '{}_coeff_data_{}.txt'.format(model.name, channel)), delim_whitespace=True)

    uniq_A = df.A.unique()
    uniq_E = df.E.unique()
    uniq_H = df.H.unique()

    color_arr = ['m', 'b', 'r', 'g', 'c', 'k', 'y']
    colors = {}
    for i, _ in enumerate(uniq_A):
        colors[str(_)] = color_arr[i]

    size_low = 50
    size_hi  = 100
    size_step = (size_hi - size_low) / len(uniq_E)
    sizes = {}
    for i, _ in enumerate(uniq_E):
        sizes[str(_)] = size_low + (i * size_step)

    angle_low = 0.
    angle_hi  = 180.
    angle_step = (angle_hi - angle_low) / len(uniq_H)
    angles = {}
    for i, _ in enumerate(uniq_H):
        angles[str(_)] = -i * angle_step

    m = mpl.markers.MarkerStyle(marker='2')

    df.insert(0, 'C', 1)
    A = df.iloc[:,:10].to_numpy()
    B = np.linalg.pinv(A)

    a = []
    c = []
    t = []
    for i in range(4):
        a.append(df.loc[:, 'a{}'.format(i)].to_numpy())
        c.append(np.dot(B, a[i]))
        t.append(np.dot(A, c[i]))

    plt.figure(figsize=[15, 15])
    for i in range(4):
        low = np.min([a[i], t[i]])
        hi  = np.max([a[i], t[i]])
        plt.subplot(2, 2, i+1)
        for j, (aij, tij) in enumerate(zip(a[i], t[i])):
            m._transform = m.get_transform().rotate_deg(angles[str(df.loc[j, 'H'])])
            plt.scatter(aij, tij, 
                    color=colors[str(df.loc[j, 'A'])],
                    marker=m,
                    s=sizes[str(df.loc[j, 'E'])])
                    #facecolors='none')

        plt.plot([low,hi], [low,hi])
        plt.title('a{}'.format(i))

    plt.gcf().suptitle('{} {}'.format(model.name, channel))


def density_check(infile, nkg_models, alt_models, savedir='fit_check', batchmode=True):

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
    markercolors = [R.kBlack, R.kBlue, R.kCyan, R.kGreen, R.kViolet]


    f = R.TFile(infile)
    canvas = plot.make1Dcanvas('check_canvas')
    canvas.SetLogy()
    canvas.SetLogx()
    
    fnkg = R.TF1(NKG.name, NKG.fitexpr, __XMIN__, __XMAX__)
    falt = R.TF1(Alt.name, Alt.fitexpr, __XMIN__, __XMAX__)
   
    fnkg.SetLineColor(R.kOrange - 3)
    fnkg.SetLineWidth(2)
    fnkg.SetLineStyle(1)
    falt.SetLineColor(R.kGreen + 2)
    falt.SetLineWidth(2)
    falt.SetLineStyle(1)

    for channel in ['muon', 'photon']:
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
                    names = []
                    min_ = None
                    max_ = None
                    for h, c, m in zip([dpmjet, qgsjet, qgsii, sibyll, venus], markercolors, models):
                        if (h != None):
                            h.SetLineColor(c)
                            h.SetMarkerColor(c)
                            h.SetMarkerStyle(R.kOpenCircle)
                            mn_, mx_ = plot.min_max(h)
                            if (min_ is None or mn_ < min_):
                                min_ = mn_
                            if (max_ is None or mx_ > max_):
                                max_ = mx_
                            hists.append(h)
                            names.append(m)

                    for i, h in enumerate(hists):
                        if (h.GetXaxis().GetXmin() > __XMIN__):
                            hists.append(hists.pop(i))
                            names.append(names.pop(i))

                    if (len(hists) > 0):
                        if (min_ is not None and max_ is not None):
                            hists[0].SetMaximum(max_ * 100.)
                            hists[0].SetMinimum(min_ / 100.)
                            hists[0].GetXaxis().SetRangeUser(__XMIN__, __XMAX__)

                        A = MassNumber[primary]
                        E = float(energy)
                        H = alt_km[int(altitude)]
                
                        for i in range(4):
                            fnkg.SetParameter(i, Coefficients.get(i, nkg_models[channel], A, E, H))
                            falt.SetParameter(i, Coefficients.get(i, alt_models[channel], A, E, H))

                        legend = plot.make1Dlegend()
                        legend.AddEntry(fnkg, NKG.name, 'l')
                        legend.AddEntry(falt, Alt.name, 'l')
                        for h, n in zip(hists, names):
                            legend.AddEntry(h, n, 'p')

                        canvas.cd()
                        canvas.Clear()
                        hists[0].GetXaxis().SetRangeUser(__XMIN__, __XMAX__)
                        hists[0].Draw('hist pe')
                        for h in hists[1:]:
                            h.Draw('hist pe same')
                        fnkg.Draw('same')
                        falt.Draw('same')
                        legend.Draw()
                        canvas.Update()
                        canvas.SaveAs('{}/check_{}_{}_{}_{}.png'.format(savedir, channel, primary, energy, altitude))


def supervised_density(infile, automate=False, edit_fit=None, savedir='fit_density', model=Alt, channels=['muon', 'photon']):

    if (automate and edit_fit is not None):
        print(' [ERROR]: cannot automate a fit edit')
        return

    if (not os.path.isdir(savedir)):
        os.makedirs(savedir, exist_ok=True)

    thin      = '1E-6'
    N         = '100'
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
        file.write( ('{:15s}'*8).format('A', 'E', 'H', 'a0', 'a1', 'a2', 'a3', 'chi2') + '\n')

    def write_data(file, primary, energy, altitude, a0, a1, a2, a3, chi2, flush=True):
        A = MassNumber[primary]
        E = float(energy)
        H = alt_km[int(altitude)]
        file.write(('{:<15.6e}'*8).format(A, E, H, a0, a1, a2, a3, chi2) + '\n')
        if (flush):
            file.flush()

    def edit_data(file, primary, energy, altitude, a0, a1, a2, a3, chi2):
        A = MassNumber[primary]
        E = float(energy)
        H = alt_km[int(altitude)]
        
        lines = file.readlines()
        lines[edit_fit + 1] = ('{:<15.6e}'*8).format(A, E, H, a0, a1, a2, a3, chi2) + '\n'
        file.close()
        file = open(file.name, 'w')
        file.writelines(lines)

    def delete_data(file):
        lines = file.readlines()
        lines.pop(edit_fit + 1)
        file.close()
        file = open(file.name, 'w')
        file.writelines(lines)

    total_fits = len(channels) * len(primaries) * len(energies) * len(altitudes)
    fit_count  = 0
    for channel in channels:
        mode = 'w'
        if (edit_fit is not None):
            mode = 'r'
        coeff_data = open('{}/{}_coeff_data_{}.txt'.format(savedir, model.name, channel), mode)
        if (mode == 'w'):
            write_header(coeff_data)

        for primary in primaries:
            for energy in energies:
                for altitude in altitudes:
                    fit_count += 1

                    if (edit_fit is not None and fit_count != edit_fit + 1):
                        continue

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
                            if (h.GetXaxis().GetXmin() >= 1e-3):
                                wrong_min.append(h)
                            else:
                                hists.append(h)

                    if (len(wrong_min) > len(hists)):
                        with open('{}/{}_info.txt'.format(savedir, model.name), 'w') as g:
                            info = ' [INFO]: switching to wrong min for {}'.format(name)
                            info += '\n         {} vs {}'.format(len(wrong_min), len(hists))
                            g.write(info + '\n')
                            print(info)
                        hists = wrong_min

                    if (len(hists) > 0):
                        hsum = hists[0]
                        count = 1
                        hNum = hists[0].GetXaxis().GetNbins()
                        hMin = hists[0].GetXaxis().GetXmin()
                        hMax = hists[0].GetXaxis().GetXmax()
                        for h in hists[1:]:
                            x = h.GetXaxis()
                            if (x.GetNbins() == hNum and x.GetXmin() == hMin and x.GetXmax() == hMax):
                                hsum.Add(h)
                                count += 1
                        hsum.Scale(1./float(count))

                        min_, max_ = plot.min_max(hsum)
                        if (min_ is not None and max_ is not None):
                            hsum.SetMaximum(max_ * 100.)
                            hsum.SetMinimum(min_ / 100.)
                        hsum.SetMarkerColor(R.kBlack)
                        hsum.SetMarkerStyle(R.kFullCircle)
                        hsum.SetLineColor(R.kRed)

                        xaxis = hsum.GetXaxis()
                        for stop_bin in range(xaxis.GetNbins(), 0, -1):
                            val = hsum.GetBinContent(stop_bin)
                            if (val > 0.):
                                break
                        limit = xaxis.GetBinCenter(stop_bin)

                        f1 = R.TF1('fitfunc', model.fitexpr, 1e-3, limit)
                        for i in range(4):
                            f1.SetParameter(i, model.a_set[i])
                            f1.SetParLimits(i, model.a_min[i], model.a_max[i])

                        f1.SetLineColor(R.kBlue)
                        f1.SetLineStyle(1)
                        f1.SetLineWidth(2)

                        legend = plot.make1Dlegend()
                        legend.AddEntry(hsum, '{}: n={}, thin={}, ave'.format(plot.__symbols__[channel], N, thin), 'p')
                        legend.AddEntry(f1, model.printexpr, 'l')

                        def rolldice(f1):
                            def rolldie(mu, low, hi):
                                sig = (hi - low) / 10.
                                while (True):
                                    val = np.random.normal(mu, sig)
                                    if (low < val and val < hi):
                                        return val

                            for i in range(4):
                                f1.SetParameter(i, rolldie(model.a_set[i], model.a_min[i], model.a_max[i]))

                        def query(f1):
                            while (True):
                                print()
                                print('fitting {} of {} for: {}_{}_{}_{}_{}'.format(fit_count, total_fits, model.name, channel, primary, energy, altitude))
                                print('current chi2 = {:.1f}'.format(f1.GetChisquare()))
                                if (automate):
                                    action = 's'
                                else:
                                    if (edit_fit is None):
                                        action = input('>>> (r)etry, (s)ave, (a)bandon: ')
                                    else:
                                        action = input('>>> (r)etry, (s)ave, (a)bandon, (d)elete: ')
                                if (action == 'r'):
                                    return 'r'
                                if (action == 's'):
                                    canvas.SaveAs('{}/{:03d}_{}_{}_{}_{}_{}.png'.format(savedir, fit_count-1, model.name, channel, primary, energy, altitude))
                                    fit0 = f1.GetParameter(0)
                                    fit1 = f1.GetParameter(1)
                                    fit2 = f1.GetParameter(2)
                                    fit3 = f1.GetParameter(3)
                                    chi2 = f1.GetChisquare()
                                    if (edit_fit is None):
                                        write_data(coeff_data, primary, energy, altitude, fit0, fit1, fit2, fit3, chi2) 
                                    else:
                                        edit_data(coeff_data, primary, energy, altitude, fit0, fit1, fit2, fit3, chi2)
                                    return 's'
                                if (action == 'a'):
                                    return 'a'
                                if (edit_fit is not None and action == 'd'):
                                    return 'd'
                                print(' input invalid, try again ..')

                        trys = 0
                        best_chi2 = None
                        while (True):
                            fitstatus = int(hsum.Fit('fitfunc', 'RE Q'))
                            if (best_chi2 is None):
                                best_chi2 = f1.GetChisquare()
                            elif (best_chi2 > f1.GetChisquare()):
                                best_chi2 = f1.GetChisquare()

                            if (fitstatus != 0 and trys < 100):
                                trys += 1
                                rolldice(f1)
                                continue
                            elif (f1.GetChisquare() > 1e5 and trys < 200):
                                trys += 1
                                rolldice(f1)
                                continue
                            elif (trys < 10):
                                trys += 1
                                rolldice(f1)
                                continue
                            elif (f1.GetChisquare() != best_chi2):
                                trys += 1
                                rolldice(f1)
                                continue

                            fit = f1.GetParameter(0)
                            fit1 = f1.GetParameter(1)
                            fit2 = f1.GetParameter(2)
                            fit3 = f1.GetParameter(3)

                            print()
                            for i in range(4):
                                fit = f1.GetParameter(i)
                                if (fit == model.a_min[i] or fit == model.a_max[i]):
                                    print(' ***** a{} at limit: {}'.format(i, fit))

                            canvas.cd()
                            canvas.Clear()
                            hsum.GetXaxis().SetRangeUser(__XMIN__,__XMAX__)
                            hsum.Draw('hist pe')
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
                            if (action == 'd'):
                                delete_data(coeff_data)
                                break
        coeff_data.close()


