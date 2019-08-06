reset

set title "radio pulses"
set xlabel "time [ns]"
set ylabel "field strength [microV/m]"
set style data lines

set key top right

a = 1.e9		#s to ns
b = 2.99792458e10	#cgs to microVolt/m

#set xrange [*:500]
#set yrange [.05:2000]

set mxtics 10
set mytics 10

plot	"SIM000001_coreas/raw_pole_100m_0deg.dat" using (a*$1):(b*$3) w l title "Ewest" lt 1, \
        "SIM000001_coreas/raw_pole_100m_0deg.dat" using (a*$1):(b*$2) w l title "Enorth" lt 3, \
        "SIM000001_coreas/raw_pole_100m_0deg.dat" using (a*$1):(b*$4) w l title "Eup" lt 4
