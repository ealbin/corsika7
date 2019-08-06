#!/bin/bash -l
#
########## ik3 script for primary run (corsika-initial.ik3.sh) ##########
#          using /cr/data02/joe/corsika.trunk/run/ as main path.
#
#$ -S /bin/bash 
#$ -cwd
#$ -j y
#$ -e /cr/data02/joe/corsika.trunk/run/
#$ -o /cr/data02/joe/corsika.trunk/run/
#$ -v FILE,EXEC,JOBINFO,mpiid,RNNUM,CHECK
#
# set -x
########## corsika_run first command block
  oput=`echo $RNNUM $mpiid | awk '{printf("./DAT%06d-000000000-%09d.lst", $1,$2)}'`
  end=`echo $RNNUM $mpiid | awk '{printf("./DAT%06d-000000000-%09d.end", $1,$2)}'`
  t=$(date '+%s');
  echo $mpiid " 0 1 " $RNNUM " 0 " $mpiid " 0 0 " $t > ./status_start
  a=$(/usr/bin/time -f "\t%E Real-time \t%S TotalCPUseconds" 2>&1 ./$EXEC < $FILE > $oput);
  a=$(date '+%s')$a;
  echo $a $oput >> $JOBINFO
########## postprocess second command block 
  t=$(date '+%s');
  echo $mpiid " 0 1 " $RNNUM " 0 " $mpiid " 0 0 " $t > ./status_finish
  echo "Initial Run" | cat > $end
  echo "Initial Run" | cat > $CHECK
