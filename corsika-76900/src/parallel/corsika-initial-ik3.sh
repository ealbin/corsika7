#!/bin/bash -l
#
###################################################################################
#       Scripts are developed for parallelized run of CORSIKA in `ik3` at KIT
#       using /cr/data02/joe/corsika.trunk/run/ as main path.
# Usage:
# ./corsika-initial-ik3.sh parallel-002413 corsika74088_stnd_QGSII4_gheisha_pll ik3
###################################################################################
# set -x

########## test on correct number of arguments ##########
if [ $# -lt 2 ] ; then
   echo "Submitting parallel corsika simulation by consecutive scripts:"
   echo "Execute preprocessing by runninmg script preprocess-ik3.sh:"
   echo "./preprocess-ik3.sh parallel-002413 renamedexecutable"
   echo "In csk002413/ start the parallel corsika simulation by typing:"
   echo "./corsika-initial-ik3.sh parallel-002413 renamedexecutable ik3"

########## correct number of arguments (2 or 3) ##########
else

   mpiid=1
   numrunning=0
   c=`grep RUNNR $1 | awk '{print $2}'`
   inputfile=`echo $c $mpiid | awk '{printf("./DAT%06d-000000000-%09d.inp", $1,$2)}'`
   cp $1 $inputfile
   touch ./job-file

   ########## submit the parallel-inputs or new parallel-inputs ##########
   Exec=`echo $2`
   export JOBINFO=./job-file
   export EXEC=`echo $2`
   export FILE=$inputfile
   export TYPE=`echo $3`
   export mpiid
   export RNNUM=$c
   export CHECK=./status_check-file

   case $TYPE in
       "PBS")      
                echo "under construction"
                qsub corsika-initial.pbs.sh
                exit 2;;
       "ik3")      
                numrunning=`expr $numrunning + 1`
                qsub corsika-initial.ik3.sh;;
       "LL")      
                numrunning=`expr $numrunning + 1`
                llsubmit corsika-initial.ll.sh;;
       "BOINC")    
                echo "not available"
                exit 2;;
       "GRID")   
                echo "not tested, but could be done"
                glite-job-submit --vo auger-vo -r remotegrid -o jobid.$c jdlscript > jobsbmt.output 2>&1
                exit 2;;
       *)      :
                numrunning=`expr $numrunning + 1`
                ./corsika-initial.local.sh & 
                ;;
   esac

   ########## generate secondary files and submit them ##########
   primoutfile=`echo $c $mpiid | awk '{printf("./DAT%06d-000000000-%09d.end", $1,$2)}'`
   while [ ! -s $primoutfile ] # wait until initial submission is completed #
      do
         sleep 4
      done
   e=`echo $c | awk '{printf("%06d", $1)}'`
   touch ./status_check-file
   export numrunning

   ########## loop will generate secondaries from others, submit them, continue till end ##########
   call_loop()
   {
   for i in $(find . -name "*.job")
      do
         endcheck=`grep -c "END" $i`
         if [ $endcheck -eq 1 ]
         then
            sed -i '/END/d' $i
            parentid=`echo $i | sed 's/.*-\(.*\).job/\1/'`
            seed=`echo $i | sed 's/.*-\(.*\)-.*/\1/'`
            source corsika-secondary-ik3.sh
            function_start $i $inputfile $seed $parentid $Exec $e
            mv $i './'$i'done'
         fi
      done
   for k in $(find . -name "*.end")
      do
         numrunning=`expr $numrunning - 1`
         rm $k
      done
   }
   ########## test on file names *.job ##########
   numjobfiles=$( find . -name '*.job' | wc -l)
   while [ $numrunning -ne 0 -o $numjobfiles -ne 0 ]
   do
      if [ $numjobfiles -eq 0 ]
      then
         sleep 3
         call_loop
      else   
         call_loop
      fi
      numjobfiles=$( find . -name '*.job' | wc -l)
   done
fi
