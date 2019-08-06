#!/bin/bash
#
# high energy simulations to a single huge particle data file
#      done on kit.scc.campus.north.grid (opus cluster)
# ===========================================================
#      /data/corsdat7/joe/e17m100ra
#      /data/corsdat7/joe/e17m316ra
#      /data/corsdat7/joe/e18m100ra
#      /data/corsdat8/joe/e18m100ra
#
# to get infos from corsika files DAT38*.lst about CPU times:
# -----------------------------------------------------------
#
  for i in $( ls -1 DAT38*.lst ); do ls -1 $i ; grep "PRESENT TIME :" $i | head -3 | tail -2 ; done
#
# under /bin/csh:
# ----------------
# foreach i ( `ls -1 DAT38*.lst` )
#    echo $i
#    grep "PRESENT TIME :" $i | head -3 | tail -2
# end    
