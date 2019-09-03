#!/bin/bash -l
###########  loadleveler script for primary run (corsika-initial.ll.sh) #########
#@ step_name = corsika_run
#@ job_type = parallel
#@ output = ./corsika-initial.$(jobid).$(stepid).out
#@ error = ./corsika-initial.$(jobid).$(stepid).err
#@ environment = $FILE;$EXEC;$JOBINFO;$mpiid;$RNNUM;$CHECK;
#@ class = small
#@ blocking = unlimited
## requirements = ( Feature=="opteron" )
#@ requirements = ( Feature=="penryn" )
#@ total_tasks = 1
#@ queue
#@ step_name = postprocess
#@ dependency = (corsika_run == 0)
#@ job_type = serial
#@ class = small
#@ queue

case $LOADL_STEP_NAME in

   corsika_run)
      oput=`echo $RNNUM $mpiid | awk '{printf("./DAT%06d-000000000-%09d.lst", $1,$2)}'`
      end=`echo $RNNUM $mpiid | awk '{printf("./DAT%06d-000000000-%09d.end", $1,$2)}'`
      t=$(date '+%s');
      echo $mpiid " 0 1 " $RNNUM " 0 " $mpiid " 0 0 " $t > ./status_start
      a=$(/usr/bin/time -f "\t%E Real-time \t%S TotalCPUseconds" 2>&1 ./$EXEC < $FILE > $oput);
      a=$(date '+%s')$a;
      echo $a $oput >> $JOBINFO
      ;;

   postprocess) 
      oput=`echo $RNNUM $mpiid | awk '{printf("./DAT%06d-000000000-%09d.lst", $1,$2)}'`
      end=`echo $RNNUM $mpiid | awk '{printf("./DAT%06d-000000000-%09d.end", $1,$2)}'`
      t=$(date '+%s');
      echo $mpiid " 0 1 " $RNNUM " 0 " $mpiid " 0 0 " $t > ./status_finish
      echo "Initial Run" | cat > $end
      echo "Initial Run" | cat > $CHECK
      ;;

   *)    
      echo "Nothing to do for $LOADL_STEP_NAME"
      ;;

esac

