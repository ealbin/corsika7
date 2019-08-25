
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
Spectrum['xbins']  = 1000
Spectrum['xmin']   =  3. # 10^(xmin) [eV]
Spectrum['xmax']   = 21. # 10^(xmax) [eV]
Spectrum['width']  = 800
Spectrum['height'] = 800


Content = {}
Content['name']   = 'content'
Content['title']  = ''
Content['ytitle'] = 'Counts'
Content['stats']  = False
Content['xbins']  = 37
Content['width']  = 800
Content['height'] = 800


Impact = {}
Impact['name'] = 'impact'
Impact['title'] = ''
Impact['ytitle'] = 'Altitude [km]'
Impact['stats']  = False
Impact['xbins']  = 3
Impact['ybins']  = 1000
Impact['xmin']   = 0
Impact['ymin']   = -3. # 10^(ymin) [km]
Impact['xmax']   = 3
Impact['ymax']   = 2.3 # 10^(ymax) [km]
Impact['width']  = 800
Impact['height'] = 800

