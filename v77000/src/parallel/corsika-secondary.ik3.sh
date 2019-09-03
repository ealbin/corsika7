#!/bin/bash -l
#
########## ik3 script for secondary runs (corsika-secondary.ik3.sh) ##########
#   using /cr/auger02/joe/corsika.trunk/run/ as main path; check it.
#
#$ -S /bin/bash 
#$ -cwd
#$ -j y
#$ -e /cr/auger02/joe/corsika.trunk/run/
#$ -o /cr/auger02/joe/corsika.trunk/run/
#$ -v FILE,NR,CHECK,NUM,EXEC,JOBINFO,PARENT,seed,id1,id2,format,mpiid,RNNUM
#  $ -m be
#  $ -M juergen.oehlschlaeger@kit.edu
#
# set -x
########## corsika_run:
  t=$(date '+%s');
  echo $mpiid $PARENT " 0 " $RNNUM $seed $mpiid $id1 $id2 $t >> ./status_start
  oput=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.lst", $1,$2)}'`
  end=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.end", $1,$2)}'`
  a=$(/usr/bin/time -f "\t%E Real-time \t%S TotalCPUseconds" 2>&1 ./$EXEC < $FILE > $oput);
  a=$(date '+%s')$a;
  echo $a $oput >> $JOBINFO
########## postprocess:
  oput=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.lst", $1,$2)}'`
  end=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.end", $1,$2)}'`
  while [ `grep -c "END OF RUN" $oput` -ne 1 ]
     do
        ./$EXEC < $FILE > $oput
     done
  t=$(date '+%s');    echo $mpiid $PARENT " 0 " $RNNUM $seed $mpiid $id1 $id2 $t >> ./status_finish
  touch $end
  echo "Secondary Run ID" $mpiid | cat >> $CHECK
  echo "Secondary Run ID" $mpiid | cat >> $end
