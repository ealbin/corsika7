
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.special import gamma
from scipy.sparse import dok_matrix
from scipy.interpolate import CubicSpline

import ROOT as R

from . import plot
from . import settings


class Transform:
    def A(A):
        return np.log(A + 1.)
    def Ainv(a):
        return np.exp(a) - 1.

    def E(E):
        return np.log10(E / 1e18)
    def Einv(e):
        return np.power(10., e) * 1e18

    def H(H):
        return 100. * (1. - np.log10(10. - H/10.)) 
    def Hinv(h):
        return 10.*(10. - np.power(10., 1. - h/100.)) 


class Coefficients:
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

    indir = 'fit_density/'
    alt_models = {'photon':'fit_models/Alt_photon_pol3_model.txt',
                  'muon':'fit_models/Alt_muon_pol3_model.txt'}

    eff_photon = 0.0001
    eff_muon   = 0.5

    __xmin__ = 1e-24 # km
    __xmax__ = 1e4   # km


    def __init__(self):
        self.alt = R.TF1('Alt_muon', Alt.fitexpr, Density.__xmin__, Density.__xmax__)


    def set_shower(self, mass_number, energy_eV, altitude_km):
        self.A = mass_number
        self.E = energy_eV
        self.H = altitude_km

        x = np.logspace(np.log10(Density.__xmin__), np.log10(Density.__xmax__), 10000)
        
        for i in range(4):
            self.alt.SetParameter(i, Coefficients.get(i, Density.alt_models['photon'], self.A, self.E, self.H))
        y1 = np.zeros(x.size)
        y2 = np.zeros(x.size)
        for i in range(x.size):
            y1[i] = self.alt.Eval(x[i])
            y2[i] = x[i] * y1[i]
        y2 = 2.*np.pi * 1e10 * y2
        self.photon  = CubicSpline(x, y1)
        photon2 = CubicSpline(x, y2)
        y3 = np.zeros(x.size)
        for i in range(x.size):
            y3[i] = photon2.integrate(0., x[i])
        self.N_photons = y3[-1]
        y3 = y3 / y3[-1]
        cut = -1
        while (y3[cut] == 1.):
            cut -= 1
        cut += 2
        self.photon_cdf = CubicSpline(y3[:cut], x[:cut])


        for i in range(4):
            self.alt.SetParameter(i, Coefficients.get(i, Density.alt_models['muon'], self.A, self.E, self.H))
        y1 = np.zeros(x.size)
        y2 = np.zeros(x.size)
        for i in range(x.size):
            y1[i] = self.alt.Eval(x[i])
            y2[i] = x[i] * y1[i]
        y2 = 2.*np.pi * 1e10 * y2
        self.muon  = CubicSpline(x, y1)
        muon2 = CubicSpline(x, y2) 
        y3 = np.zeros(x.size)
        for i in range(x.size):
            y3[i] = muon2.integrate(0., x[i])
        self.N_muons = y3[-1]
        y3 = y3 / y3[-1]
        cut = -1
        while (y3[cut] == 1.):
            cut -= 1
        cut += 2
        self.muon_cdf = CubicSpline(y3[:cut], x[:cut])


class Detection:

    swidth  = 6.25  / 100. / 1000. # km
    sheight = 11.25 / 100. / 1000. # km 

    Aratio = 0.15 / (6.25 * 11.25)
    eff_muon   = 0.5
    eff_photon = 0.0001

    R = 70. # km


    def __init__(self):
        self.density = Density()
        self.xoffset, self.yoffset = self.xy2index(Detection.R, Detection.R)

    def xy2index(self, x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        ix = np.sign(x) * np.floor( np.abs(x + Detection.swidth/2.)  / Detection.swidth )
        iy = np.sign(y) * np.floor( np.abs(y + Detection.sheight/2.) / Detection.sheight )
        return np.asarray(ix, dtype=int), np.asarray(iy, dtype=int)

    def index2key(self, ixiy):
        return ixiy[0] + self.xoffset, ixiy[1] + self.yoffset

    def key2index(self, key):
        return key[0] - self.xoffset, key[1] - self.yoffset

    def index2xy(self, ixiy):
        ix = ixiy[0]
        iy = ixiy[1]
        x = ix * Detection.swidth
        y = iy * Detection.sheight
        return np.asarray(x, dtype=float), np.asarray(y, dtype=float)

    def set_shower(self, mass_number, energy_eV, altitude_km):
        self.density.set_shower(mass_number, energy_eV, altitude_km)

    def run_shower(self, step=10000000):
        #self.phone_hits  = dok_matrix( (2*self.xoffset, 2*self.yoffset), dtype=int )
        #self.sensor_hits = dok_matrix( (2*self.xoffset, 2*self.yoffset), dtype=int )

        self.photon_hits = dok_matrix( (2*self.xoffset, 2*self.yoffset), dtype=int )
        self.muon_hits   = dok_matrix( (2*self.xoffset, 2*self.yoffset), dtype=int )

        N_photons = self.density.N_photons
        N_muons = self.density.N_muons

        def get_keys(cdf, size):
            r = cdf(np.random.random(size))
            theta = 2.*np.pi * np.random.random(size)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            return self.index2key( self.xy2index(x, y) )

        def save_sparse(keys, arr):
            ix, iy = keys
            hit = np.random.choice([True, False], ix.size, p=[self.p_hit, 1. - self.p_hit])
            for x, y in zip(ix[hit], iy[hit]):
                arr[x, y] += 1

#            for x, y, h in zip(ix, iy, hit):
#                self.phone_hits[x, y] += 1
#                if h is True:
#                    self.sensor_hits[x, y] += 1

        self.p_hit = Detection.Aratio * Detection.eff_photon
        while (N_photons > step):
            save_sparse( get_keys(self.density.photon_cdf, step), self.photon_hits )
            N_photons -= step
        save_sparse( get_keys(self.density.photon_cdf, int(np.ceil(N_photons))), self.photon_hits )

        self.p_hit = Detection.Aratio * Detection.eff_muon
        while (N_muons > step):
            save_sparse( get_keys(self.density.muon_cdf, step), self.muon_hits )
            N_muons -= step
        save_sparse( get_keys(self.density.muon_cdf, int(np.ceil(N_muons))), self.muon_hits )

    swidth  = 6.25  / 100. / 1000. # km
    sheight = 11.25 / 100. / 1000. # km 

    Aratio = 0.15 / (6.25 * 11.25)
    eff_muon   = 0.5
    eff_photon = 0.0001

    R = 70. # km
    def save(self, path='./'):
        header = ['# mass number = {:.0f}'.format(self.density.A),
                  '# energy eV   = {:e}'.format(self.density.E),
                  '# altitude km = {:.1f}'.format(self.density.H),
                  '# R km        = {:.0f}'.format(Detection.R),
                  '#',
                  '# swidth cm   = {:.2f}'.format(Detection.swidth * 100. * 1000.),
                  '# sheight cm  = {:.2f}'.format(Detection.sheight * 100. * 1000.),
                  '# Aratio      = {:.2e}'.format(Detection.Aratio),
                  '# eff photon  = {:.1e}'.format(Detection.eff_photon),
                  '# eff muon    = {:.1e}'.format(Detection.eff_muon),
                  '#',
                  '# N photons   = {:.3e}'.format(self.density.N_photons),
                  '# N muons     = {:.3e}'.format(self.density.N_muons),
                  '#',
                  '# xoffset     = {:.0f}'.format(self.xoffset),
                  '# yoffset     = {:.0f}'.format(self.yoffset)] 

        filename = 'hits_A={:.0f}_E={:.0e}_H={:.1f}.txt'.format(self.density.A, self.density.E, self.density.H)
        with open(os.path.join(path, filename), 'w') as f:
            for line in header:
                f.write(line + '\n')
            f.write('\n')
            f.write('{:^12s}{:^12s}{:^12s}{:^12s}{:^12s}{:^12s}\n'.format('particle', 'ix', 'iy', 'x', 'y', 'counts'))
            for k, v in self.photon_hits.items():
                x, y = self.index2xy( self.key2index(k) )
                f.write('{:<12s}{:>12d}{:>12d}{:>12.3e}{:>12.3e}{:>12d}\n'.format('photon', k[0], k[1], x, y, v))
            for k, v in self.muon_hits.items():
                x, y = self.index2xy( self.key2index(k) )
                f.write('{:<12s}{:>12d}{:>12d}{:>12.3e}{:>12.3e}{:>12d}\n'.format('muon', k[0], k[1], x, y, v))


class ShowerContent:

    def __init__(self):
        self.density = Density()

        self.mass_numbers = [0, 1, 4, 16, 56, 238]
        self.energies = [1e18, 1e19, 1e20, 1e21, 1e22]
        self.altitude = 0.

        self.N_photons = []
        self.N_muons = []

        for A in self.mass_numbers:
            for E in self.energies:
                print('A = {}, E = {}'.format(A, E))
                self.density.set_shower(A, E, self.altitude)
                self.N_photons.append({'A':A, 'E':E, 'N':self.density.N_photons})
                self.N_muons.append({'A':A, 'E':E, 'N':self.density.N_muons})


    def plot(self, save=None):
        plt.rc('font', size=16)
        plt.figure(figsize=[9,5])

        x_0 = []
        y_0 = []
        x_238 = []
        y_238 = []
        for p, m in zip(self.N_photons, self.N_muons):
            if (p['A'] == 0):
                x_0.append(p['N'])
                y_0.append(m['N'])
            elif (p['A'] == 238):
                x_238.append(p['N'])
                y_238.append(m['N'])
        plt.plot(x_0, y_0, 'k--')
        plt.plot(x_238, y_238, 'k--')

        done_0 = False
        done_238 = False
        for p, m in zip(self.N_photons, self.N_muons):

            facecolors = 'w'
            label = None 
            if (p['A'] == 0):
                marker = 'o'
                if (done_0 is False):
                    label = r'$\gamma$'
                    done_0 = True
            elif (p['A'] != 238):
                marker = '.'
                facecolors = 'k'
            else:
                marker = 's'
                if (done_238 is False):
                    label = r'$^{238}$U'
                    done_238 = True

            plt.scatter(p['N'], m['N'], marker=marker, facecolors=facecolors, edgecolor='k', s=80, label=label, alpha=1.)

            if (p['A'] == 238):
                plt.text(p['N'] / 1.5, 2*m['N'], r'$10^{{{:.0f}}}$'.format(np.log10(p['E'])))

        xmin = self.N_photons[0]['N'] / 10.
        xmax = self.N_photons[-1]['N'] * 10.
        ymin = self.N_muons[0]['N'] / 10.
        ymax = self.N_muons[-1]['N'] * 10.

        plt.legend(frameon=False)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        plt.xlabel('Number of Photons')
        plt.ylabel('Number of Muons')
        plt.tight_layout()
        plt.show()
        if (save is not None):
            plt.savefig(save)

