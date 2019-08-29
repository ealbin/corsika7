
def mmass2GeV(mmass):
    N_A = 6.0221409e23   # Avogadro's number
    c   = 299_792_458.   # speed of light [m/s]
    e   = 1.60217662e-19 # unit of charge [C]

    kg  = mmass / N_A / 1000. # [kg/atom]
    J   = kg * c * c          # E = m c^2 [J]
    eV  = J / e               # [eV]
    return eV / 1e9           # [GeV]

# Nuclei IDs are A x 100 + Z
# But this list only considers the "Z" part, e.g. the ones and tens place
ID = {}
ID['01'] = {'catagory':'nuclei', 'symbol':'H',  'mass':mmass2GeV(1.0079)}
ID['02'] = {'catagory':'nuclei', 'symbol':'He', 'mass':mmass2GeV(4.0026)}
ID['03'] = {'catagory':'nuclei', 'symbol':'Li', 'mass':mmass2GeV(6.941) }
ID['04'] = {'catagory':'nuclei', 'symbol':'Be', 'mass':mmass2GeV(9.0122)}
ID['05'] = {'catagory':'nuclei', 'symbol':'B',  'mass':mmass2GeV(10.811)}
ID['06'] = {'catagory':'nuclei', 'symbol':'C',  'mass':mmass2GeV(12.011)}
ID['07'] = {'catagory':'nuclei', 'symbol':'N',  'mass':mmass2GeV(14.067)}
ID['08'] = {'catagory':'nuclei', 'symbol':'O',  'mass':mmass2GeV(15.999)}
ID['09'] = {'catagory':'nuclei', 'symbol':'F',  'mass':mmass2GeV(18.998)}
ID['10'] = {'catagory':'nuclei', 'symbol':'Ne', 'mass':mmass2GeV(20.180)}
ID['11'] = {'catagory':'nuclei', 'symbol':'Na', 'mass':mmass2GeV(22.990)}
ID['12'] = {'catagory':'nuclei', 'symbol':'Mg', 'mass':mmass2GeV(24.305)}
ID['13'] = {'catagory':'nuclei', 'symbol':'Al', 'mass':mmass2GeV(26.982)}
ID['14'] = {'catagory':'nuclei', 'symbol':'Si', 'mass':mmass2GeV(28.086)}
ID['15'] = {'catagory':'nuclei', 'symbol':'P',  'mass':mmass2GeV(30.974)}
ID['16'] = {'catagory':'nuclei', 'symbol':'S',  'mass':mmass2GeV(32.065)}
ID['17'] = {'catagory':'nuclei', 'symbol':'Cl', 'mass':mmass2GeV(35.453)}
ID['18'] = {'catagory':'nuclei', 'symbol':'Ar', 'mass':mmass2GeV(39.948)}
ID['19'] = {'catagory':'nuclei', 'symbol':'K',  'mass':mmass2GeV(39.098)}
ID['20'] = {'catagory':'nuclei', 'symbol':'Ca', 'mass':mmass2GeV(40.078)}
ID['21'] = {'catagory':'nuclei', 'symbol':'Sc', 'mass':mmass2GeV(44.956)}
ID['22'] = {'catagory':'nuclei', 'symbol':'Ti', 'mass':mmass2GeV(47.867)}
ID['23'] = {'catagory':'nuclei', 'symbol':'V',  'mass':mmass2GeV(50.942)}
ID['24'] = {'catagory':'nuclei', 'symbol':'Cr', 'mass':mmass2GeV(51.996)}
ID['25'] = {'catagory':'nuclei', 'symbol':'Mn', 'mass':mmass2GeV(54.938)}
ID['26'] = {'catagory':'nuclei', 'symbol':'Fe', 'mass':mmass2GeV(55.845)}

