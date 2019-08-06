#/bin/bash -l
#
########## create secondary files and submit them to job queue ##########
#       using /cr/data02/joe/corsika.trunk/run/ as main path.
# set -x

function_start()
{
   if (( $# != 6 )) ;
   then
      echo "script to generate secondary input files and change run number"
      echo -e "usage: To run the script write corsika-secondary.sh DAT.JOB_file_name input_file_name seed_of_DAT parent_id executable format"
   else
      export fname=`echo $1`
      export file=`echo $2`
      export seed=`echo $3`
      export parentid=`echo $4`
      export format=`echo $6`
      while read line
         do
            export mpiid=`expr $mpiid + 1`
            inputfile=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.inp", $1,$2)}'`
            echo $line > $inputfile
            cat $file >> $inputfile
            sed -i '/CUTFILE/ s/CUTFILE/\nCUTFILE/'  $inputfile
         # remove duplicate lines from input file containing parallel and cutfile keywords 
         { rm $inputfile && awk 'BEGIN{ nodup["PARALLEL"] nodup["CUTFILE"] } !_[$1]++ || !($1 in nodup)' > $inputfile; } < $inputfile
            oldmpiid=`grep PARALLEL $inputfile | awk '{print $4}'`
            sed -i "/PARALLEL/ s/ $oldmpiid / $mpiid /" $inputfile
            id1=`grep CUTFILE $inputfile | awk '{print $3}'`
            id2=`grep CUTFILE $inputfile | awk '{print $4}'`
            export id1
            export id2
            export EXEC=`echo $Exec`
            export FILE=$inputfile
            export CHECK=./status_check-file
            export JOBINFO=./job-file
            export PARENT=$oldmpiid
      case $TYPE in
            "PBS")    
                     echo "under construction"
                     qsub corsika-initial.pbs.sh
                     exit 2;;
            "ik3")
                     numrunning=`expr $numrunning + 1`
                     qsub corsika-secondary.ik3.sh;;
            "LL")                
                     numrunning=`expr $numrunning + 1`
                     llsubmit corsika-secondary.ll.sh;;
            "BOINC")                 
                     echo "not available"
                     exit 2;;
            "GRID")        
                     echo "not tested, but could be done"
                     glite-job-submit --vo auger-vo -r remotegrid -o jobid.$c jdlscript > jobsbmt.output 2>&1
                     exit 2;;
            *)        
                     export numrunning=`expr $numrunning + 1`
                     ./corsika-secondary.local.sh &
                     jobs -l
                     ;;
      esac
   done < $fname
fi
}
