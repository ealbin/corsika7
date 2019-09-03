#!/bin/bash -l
########## loadleveler script for secondary runs (corsika-secondary.ll.sh) #########
#@ job_type = parallel
#@ step_name = corsika_run
#@ output = ./corsika-secondary.$(jobid).$(stepid).out
#@ error = ./corsika-secondary.$(jobid).$(stepid).err
#@ environment = $FILE;$NR;$CHECK;$NUM;$EXEC;$JOBINFO;$PARENT;$seed;$id1;$id2;$format;$mpiid;$RNNUM;
#@ class = medium
## requirements = ( Feature=="penryn" )
#@ node = 1
#@ total_tasks = 1
#@ queue
#@ step_name = postprocess
#@ dependency = (corsika_run == 0)
#@ job_type = serial
#@ class = medium
#@ node_usage = shared
#@ queue
# set -x

case $LOADL_STEP_NAME in

   corsika_run)
      t=$(date '+%s');
      echo $mpiid $PARENT " 0 " $RNNUM $seed $mpiid $id1 $id2 $t >> ./status_start
      oput=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.lst", $1,$2)}'`
      end=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.end", $1,$2)}'`
      a=$(/usr/bin/time -f "\t%E Real-time \t%S TotalCPUseconds" 2>&1 ./$EXEC < $FILE > $oput);
      a=$(date '+%s')$a;
      echo $a $oput >> $JOBINFO
      ;;

   postprocess)
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
      ;;

   *) 
      echo "Nothing to do for $LOADL_STEP_NAME"
      ;;
esac
 
