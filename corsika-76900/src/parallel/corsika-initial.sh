#!/bin/bash -l

#################################################################################
#                                                                               #  
#    These scripts are developed for parallelized run of CORSIKA code           #
#                                                                               #
# Usage:                                                                        #
#   ./corsika-initial.sh "steering-inputs" "corsika-executable" "jobtype"       #
# Example:                                                                      #
#   ./corsika-initial.sh parallel-002413 corsika72495_stnd_QGSII_gheisha_pll LL #
#                                                                               #
#################################################################################
# set -x

########## test on correct number of arguments ##########
if [ $# -lt 2 ] ;
   then
      echo "script to submit parallel corsika simulations by scripts;"
      echo "an input-steering-file and a corsika-executable must exist;"
      echo "usage: to start the simulation write"
      echo "  ./corsika-initial.sh input-steering-file corsika-executable [job-distributor]"

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
                echo "under contruction"
                qsub corsika-initial.pbs.sh
                exit 2;;
       "LL")      
                numrunning=`expr $numrunning + 1`
                llsubmit corsika-initial.ll.sh;;
       "BOINC")    
                echo "comming soon"
                exit 2;;
       "GRID")   
                echo "could be done"
                glite-job-submit --vo auger-vo -r remotegrid -o jobid.$c jdlscript > jobsbmt.output 2>&1
                exit 2;;
       *)      
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
            source corsika-secondary.sh
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

