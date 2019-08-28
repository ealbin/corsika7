
import ROOT as R


Catagories = ['photon',
              'electron',
              'muon',
              'proton',
              'neutron',
              'ocharged',
              'oneutral',
              'nuclei']


Symbols = {}
Symbols['photon']   = '#gamma'
Symbols['electron'] = 'e^{#pm}'
Symbols['muon']     = '#mu^{#pm}'
Symbols['proton']   = 'p, #bar{p}'
Symbols['neutron']  = 'n, #bar{n}'
Symbols['ocharged'] = 'other^{#pm}'
Symbols['oneutral'] = 'other^{0}'
Symbols['nuclei']   = 'nuclei'


Colors = {}
Colors['photon']   = R.kRed
Colors['electron'] = R.kMagenta
Colors['muon']     = R.kGreen
Colors['proton']   = R.kBlue
Colors['neutron']  = R.kCyan
Colors['ocharged'] = R.kOrange
Colors['oneutral'] = R.kGray
Colors['nuclei']   = R.kBlack


Palette = 55
# 53 = R.kDarkBodyRadiator
# 54 = R.kBlueYellow
# 55 = R.kRainbow
# 56 = R.kInvertedDarkBodyRadiator
# 57 = R.kBird
#113 = R.kCividis


Margins = {}
Margins['left']   = 0.12
Margins['right']  = 0.05
Margins['top']    = 0.1
Margins['bottom'] = 0.12


TitleOffset = {}
TitleOffset['x'] = 1.5
TitleOffset['y'] = 1.7
TitleOffset['z'] = 1.5


Density = {}
Density['name']   = 'density_'
Density['title']  = ''
Density['xtitle'] = 'Distance from Shower Core [km]'
Density['ytitle'] = 'Number Density [1 / cm^{2}]'
Density['stats']  = False
Density['xbins']  = 30
Density['xmin']   =  -3. #0.   # [km]
Density['ymin']   = -12.       # [1 / cm^2]
Density['xmax']   =   3.  #10. # [km]  
Density['ymax']   =   1.       # [1 / cm^2]
Density['width']  = 800
Density['height'] = 800


Slice = {}
Slice['name']   = 'slice_'
Slice['title']  = ''
Slice['xtitle'] = 'West-to-East [km]'
Slice['ytitle'] = 'South-to-North [km]'
Slice['ztitle'] = 'Counts'
Slice['stats']  = False
Slice['xbins']  = 100
Slice['ybins']  = 100
Slice['xmin']   = -10. # [km]
Slice['ymin']   = -10. # [km]
Slice['xmax']   =  10. # [km]
Slice['ymax']   =  10. # [km]
Slice['zmin']   = -.1
Slice['zmax']   =  6.
Slice['width']  = 800
Slice['height'] = 800


Spectrum = {}
Spectrum['name']   = 'spectrum_'
Spectrum['title']  = ''
Spectrum['xtitle'] = 'Energy [eV]'
Spectrum['ytitle'] = 'dN/dE [Counts / eV]'
Spectrum['stats']  = False
Spectrum['xbins']  = 100
Spectrum['xmin']   =  3. # 10^(xmin) [eV]
Spectrum['ymin']   = -12.
Spectrum['xmax']   = 21. # 10^(xmax) [eV]
Spectrum['ymax']   = 1.
Spectrum['width']  = 800
Spectrum['height'] = 800


Content = {}
Content['name']   = 'content'
Content['title']  = ''
Content['ytitle'] = 'Counts'
Content['stats']  = False
Content['xbins']  = 8 #37
Content['ymin']   = -.1
Content['ymax']   = 9
Content['width']  = 800
Content['height'] = 800
Content['photon']   = 1
Content['electron'] = 2
Content['muon']     = 3
Content['proton']   = 4
Content['neutron']  = 5
Content['nuclei']   = 6
Content['ocharged'] = 7
Content['oneutral'] = 8


Impact = {}
Impact['name']   = 'impact'
Impact['title']  = ''
Impact['xtitle'] = 'First Target' 
Impact['ytitle'] = 'Altitude [km]'
Impact['stats']  = False
Impact['xbins']  = 4
Impact['ybins']  = 1000
Impact['ymin']   = -3. # 10^(ymin) [km]
Impact['zmin']   = -.1
Impact['ymax']   = 2.3 # 10^(ymax) [km]
Impact['zmax']   = 6.
Impact['width']  = 800
Impact['height'] = 800
Impact['Random']   = 0
Impact['Nitrogen'] = 1
Impact['Oxygen']   = 2
Impact['Argon']    = 3


Efficiency = {}
Efficiency['name']   = 'efficiency'
Efficiency['title']  = ''
Efficiency['xtitle'] = 'Observation Altitude [km]'
Efficiency['ytitle'] = 'Energy Efficiency of Simulation'
Efficiency['stats']  = False
Efficiency['xbins']  = 10
Efficiency['ymin']   = .8
Efficiency['ymax']   = 1.05
Efficiency['width']  = 800
Efficiency['height'] = 800

