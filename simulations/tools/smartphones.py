
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.special import gamma
from scipy.sparse import dok_matrix
from scipy.interpolate import CubicSpline, make_interp_spline

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
        self.photon2 = CubicSpline(x, y2)
        y3 = np.zeros(x.size)
        for i in range(x.size):
            y3[i] = self.photon2.integrate(0., x[i])
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
        self.muon2 = CubicSpline(x, y2) 
        y3 = np.zeros(x.size)
        for i in range(x.size):
            y3[i] = self.muon2.integrate(0., x[i])
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
            fr = cdf(np.random.random(size))
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


def plot_detection(filename, path=None):
    with open(filename) as f:
        A = int( f.readline().split('=')[-1] )
        E = float( f.readline().split('=')[-1] )
        H = float( f.readline().split('=')[-1] )
        for _ in range(4):
            f.readline()
        Aratio = float( f.readline().split('=')[-1] )
        eff_ph = float( f.readline().split('=')[-1] )
        eff_mu = float( f.readline().split('=')[-1] )


    d = Density()
    d.set_shower(A, E, H)
    
    df = pd.read_csv(filename, delim_whitespace=True, comment='#')
    lim = 10
    bins = np.linspace(-lim, lim, 20*lim + 1)

    # bins equate to 100m x 100m in size
    # white means >=100 events in bin,
    # therefore hit density is >= 100 / (100x100), or 1/10 by 1/10
    # aka, at least one event inside a 10m x 10m box

    plt.rc('font', size=16)
    plt.figure(figsize=[14,5], tight_layout=True)
    plt.subplot(1, 2, 2)

    plt.hist2d(df.x, df.y, bins=(bins, bins), norm=mpl.colors.LogNorm(vmin=1e0, vmax=1e2), cmap='Greys_r')
    plt.colorbar()
    plt.xlim(bins[0], bins[-1])
    plt.ylim(bins[0], bins[-1])
    plt.xlabel('Kilometers')
    plt.ylabel('Kilometers')
    plt.xticks(np.arange(bins[0], bins[-1] + 1, 2))
    plt.yticks(np.arange(bins[0], bins[-1] + 1, 2))
    plt.grid()


    plt.subplot(1, 2, 1)
    photons = df[df.particle == 'photon']
    muons   = df[df.particle == 'muon']

    r_photon = np.sqrt(photons.x**2 + photons.y**2)
    r_muon   = np.sqrt(muons.x**2   + muons.y**2)

    plt.hist(r_muon, bins=bins[bins >= 0], histtype='step', linewidth=2, edgecolor='k', weights=np.ones(r_muon.size)/100.)
    plt.hist(r_muon, bins=bins[bins >= 0], histtype='stepfilled', facecolor='0.8', hatch='***', weights=np.ones(r_muon.size)/100., label=r'$\mu^{\pm}$')

    plt.hist(r_photon, bins=bins[bins >= 0], histtype='step', linewidth=2, edgecolor='w', weights=np.ones(r_photon.size)/100.)
    plt.hist(r_photon, bins=bins[bins >= 0], histtype='stepfilled', facecolor='0.3', hatch='---', weights=np.ones(r_photon.size)/100., label=r'$\gamma$')

    plt.hist(np.r_[r_muon, r_photon], bins=bins[bins >= 0], histtype='step', linewidth=2, edgecolor='k', weights=np.ones(r_muon.size + r_photon.size)/100.)

    x = bins[bins > 0]
    plt.plot(x, Aratio * (eff_mu * d.muon2(x) + eff_ph * d.photon2(x)) / 1000., linewidth=2.5, c='k', label='Model')

    plt.yscale('log')
    plt.xlabel('Kilometers from Core')
    plt.ylabel('Counts / 100 m')
    plt.xlim(0, bins[-1])
    plt.ylim(1e-3, 1e3)
    plt.legend(frameon=False)
    plt.show()

    if (path is not None):
        plt.savefig(os.path.join(path, os.path.basename(filename).rstrip('.txt') + '.png'))


def get_sensitivity(filelist, trials=100, path='./'):
    
    threshold = 5
    radius = 10.
    Area = np.pi * radius**2

    d = Detection()

    test_densities = np.logspace(0, 3, 20) #np.logspace(0, 8, 9) # smartphones / km^2

    with open(os.path.join(path, 'sensitivity.txt'), 'w') as out:
        out.write('{:^12s}{:^12s}{:^12s}{:^12s}{:^12s}{:^12s}\n'.format('A', 'E', 'H', 'density', 'count', 'prob'))
        
        for i, filename in enumerate(filelist):
            print('working on {} of {}: {}'.format(i+1, len(filelist), filename), flush=True)

            with open(filename) as f:
                A = int( f.readline().split('=')[-1] )
                E = float( f.readline().split('=')[-1] )
                H = float( f.readline().split('=')[-1] )
                R = float( f.readline().split('=')[-1] )
                f.readline()
                swidth  = float( f.readline().split('=')[-1] )
                sheight = float( f.readline().split('=')[-1] )
                Aratio  = float( f.readline().split('=')[-1] )
                eff_ph  = float( f.readline().split('=')[-1] )
                eff_mu  = float( f.readline().split('=')[-1] )
                f.readline()
                N_photons = float( f.readline().split('=')[-1] )
                N_muons   = float( f.readline().split('=')[-1] )
                f.readline()
                xoffset = int( f.readline().split('=')[-1] )
                yoffset = int( f.readline().split('=')[-1] )

            d.swidth  = swidth
            d.sheight = sheight
            d.Aratio  = Aratio
            d.eff_muon   = eff_mu
            d.eff_photon = eff_ph
            d.R = R
            d.xoffset = xoffset
            d.yoffset = yoffset
            
            df = pd.read_csv(filename, delim_whitespace=True, comment='#')
            for density in test_densities:

                counts = []
                probs  = []
                for _ in range(trials):

                    N_phones = int( np.rint(Area * density) )
                    r = radius * np.random.random(N_phones)
                    theta = 2.*np.pi * np.random.random(N_phones)
                    
                    x = r * np.cos(theta)
                    y = r * np.sin(theta)

                    ix, iy = d.index2key( d.xy2index(x, y) )
                    selection = df[( df.loc[:, 'ix'].isin(ix) ) & ( df.loc[:, 'iy'].isin(iy) )]
                    counts.append( len(selection.drop_duplicates(subset=['ix', 'iy'])) )
                    if (counts[-1] >= 5):
                        probs.append(1)
                    else:
                        probs.append(0)
                
                out.write('{:>12d}{:>12.3e}{:>12.3e}{:>12.3e}{:>12.3e}{:>12.3e}\n'.format(A, E, H, density, np.mean(counts), np.mean(probs)))
                out.flush()


def plot_sensitivity(filename):
    df = pd.read_csv(filename, delim_whitespace=True)

    x = df[(df.A == 0) & (df.E == 1e15)].density
    ix = np.linspace(x.min(), x.max(), 500)

    y0_15c = df[(df.A == 0) & (df.E == 1e15)].loc[:,'count']
    y0_15p = df[(df.A == 0) & (df.E == 1e15)].loc[:,'prob']

    y0_17c = df[(df.A == 0) & (df.E == 1e17)].loc[:,'count']
    y0_17p = df[(df.A == 0) & (df.E == 1e17)].loc[:,'prob']

    y0_19c = df[(df.A == 0) & (df.E == 1e19)].loc[:,'count']
    y0_19p = df[(df.A == 0) & (df.E == 1e19)].loc[:,'prob']
    
    kp = 1

    iy0_15c = make_interp_spline(x, y0_15c)(ix)
    iy0_15p = make_interp_spline(x, y0_15p, k=kp)(ix)
    iy0_17c = make_interp_spline(x, y0_17c)(ix)
    iy0_17p = make_interp_spline(x, y0_17p, k=kp)(ix)
    iy0_19c = make_interp_spline(x, y0_19c)(ix)
    iy0_19p = make_interp_spline(x, y0_19p, k=kp)(ix)

    y238_15c = df[(df.A == 238) & (df.E == 1e15)].loc[:,'count']
    y238_15p = df[(df.A == 238) & (df.E == 1e15)].loc[:,'prob']

    y238_17c = df[(df.A == 238) & (df.E == 1e17)].loc[:,'count']
    y238_17p = df[(df.A == 238) & (df.E == 1e17)].loc[:,'prob']

    y238_19c = df[(df.A == 238) & (df.E == 1e19)].loc[:,'count']
    y238_19p = df[(df.A == 238) & (df.E == 1e19)].loc[:,'prob']

    iy238_15c = make_interp_spline(x, y238_15c)(ix)
    iy238_15p = make_interp_spline(x, y238_15p, k=kp)(ix)
    iy238_17c = make_interp_spline(x, y238_17c)(ix)
    iy238_17p = make_interp_spline(x, y238_17p, k=kp)(ix)
    iy238_19c = make_interp_spline(x, y238_19c)(ix)
    iy238_19p = make_interp_spline(x, y238_19p, k=kp)(ix)

    plt.rc('font', size=16)
    plt.figure(figsize=[14,5], tight_layout=True)

    ax = plt.subplot(1, 2, 1)
    
    ax.add_patch( mpl.patches.Rectangle((170, 1), 10, 1e6, edgecolor=None, facecolor='gray') )
    plt.plot(ix, iy0_15c, 'k-')
    plt.plot(ix, iy0_17c, 'k-')
    plt.plot(ix, iy0_19c, 'k-')
    plt.plot(ix, iy238_15c, 'k--')
    plt.plot(ix, iy238_17c, 'k--')
    plt.plot(ix, iy238_19c, 'k--')
    plt.yscale('log')
    plt.xlim(0, x.max())
    plt.ylim(1, 1e6)
    plt.xlabel('Smartphone Density [#/km$^2$]')
    plt.ylabel('Average Phones Hit')
    plt.text(400, 2e1, r'$10^{15}$')
    plt.text(600, 2e3, r'$10^{17}$')
    plt.text(800, 2e5, r'$10^{19}$')
    plt.text(200, 2e5, 'LA 1%')


    plt.subplot(1, 2, 2)
    plt.plot(ix, iy0_15p, 'k-')
    plt.plot(ix, iy0_17p, 'k-')
    plt.plot(ix, iy0_19p, 'k-')
    plt.plot(ix, iy238_15p, 'k--')
    plt.plot(ix, iy238_17p, 'k--')
    plt.plot(ix, iy238_19p, 'k--')
    plt.xscale('log')
    plt.xlim(1, x.max())
    plt.ylim(0, 1.1)
    plt.xlabel('Smartphone Density [#/km$^2$]')
    plt.ylabel('Probability of >5 Phones Hit')
    plt.text(2, .4, r'$10^{19}$')
    plt.text(20, .6, r'$10^{17}$')
    plt.text(500, .8, r'$10^{15}$')

    plt.show()

