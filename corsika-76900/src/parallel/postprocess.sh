#!/bin/bash
#
# postprocess.sh:
# ===============
#        new postprocessing for hc3 parallel corsika simulations
#        after december 2012 (because of missing statisticis by
#        the mpi parallelization version). 
# ------------------------------------------------------------------------
#   (p0) MPI job output files Job00nnnn_%j.[err,out] must exist and stay
#        in the upper directory with length 0 and length > 0 resp.; 
#        rename Job_hc3_%j.[err,out] to Job00nnnn_%j.[err,out];
#   (p1) jobhc3-00nnnn and parallel-00nnnn must exist in this path.
#   (p2) time.txt must contain 2 lines like
#       START TIME          STOP TIME       TIME (min)
#   1366262084.316990   1366263487.694263    23.389621
#   (p3) postprocess.sh must exist in this path.
#   (p4) script summ, executable summe and source code summe.f
#        must exist in this path.
#   (p5) executable totaltime and source code totaltime.f  
#        must exist in this path.
# ------------------------------------------------------------------------
# cd csk00nnnn/
#        usage (a): ./postprocess.sh
# if Job00nnnn_%j.[err,out] already copied:
#        usage (b): ./postprocess.sh -j
# ------------------------------------------------------------------------
#        data will be written to Job00nnnn_%j.[err,out] and to time.txt.
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# - - - - - - does file time.txt exist? 
ls -l time.txt* | wc | awk '{printf("%7d\n",$1)}' > jobfilerr.txt
ntimetxt=`cat jobfilerr.txt | awk '{printf("%7d\n",$1)}'`
if [ $ntimetxt -eq 0 ] ; then
 echo "           ... postprocess.sh cannot be executed in this path,"
 echo "               because file time.txt is not available."
 /bin/rm jobfilerr.txt
 exit
fi
# - - - - - - get current run number:
pwd > jobfilerr.txt
sed -i "/csk/ s/csk/   /" jobfilerr.txt 
jrunnr=`cat jobfilerr.txt | awk '{printf("%7d\n",$2)}'`
/bin/rm jobfilerr.txt
# - - - - - - does file Job00_*.err exist?
jobrunnr=`echo $jrunnr | awk '{printf("Job%.6d_*",$1)}'`
# - - - - - - test on already available jobfiles:
movejobfiles=`echo $1`
if [ -n "$movejobfiles" ]; then
  echo "           ... jobfiles $jobrunnr already available in this path."
else
  /bin/mv ../$jobrunnr .
fi
# - - - - - - now jobfiles should be here:
ls -l Job* | wc | awk '{printf("%7d\n",$1)}' > jobfilerr.txt
njobfiles=`cat jobfilerr.txt | awk '{printf("%7d\n",$1)}'`
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $njobfiles -eq 0 ] ; then
 echo "           ... jobfiles Job00*.[err,out] not found in this directory" 
 echo "           ...... and must be copied from the upper level directory."
 /bin/rm jobfilerr.txt
 # end-of case Job00* do not exist.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else   
  # - - - - - - check time.txt if postprocess.sh was already done:
  wc time.txt* | awk '{printf("%7d\n",$1)}' > jobfilerr.txt 
  ntimetxt=`cat jobfilerr.txt | awk '{printf("%7d\n",$1)}'`
  # - - - - - - get current runnr and job id:
  /bin/ls -1 Job00*.err > jobfilerr.txt
  sed -i "/Job/ s/_/ /" jobfilerr.txt
  sed -i "/Job/ s/\./ /" jobfilerr.txt
  sed -i "/Job/ s/err//" jobfilerr.txt
  sed -i "/Job/ s/Job//" jobfilerr.txt
  jrunnr=`cat jobfilerr.txt | awk '{printf("%06d",$1)}'`
  jobnumb=`cat jobfilerr.txt | awk '{printf("%d",$2)}'`
  /bin/rm jobfilerr.txt
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if [ $ntimetxt -gt 2 ] ; then
    echo "           ... postprocess.sh already done for runnr" $jrunnr
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else
    /bin/rm -f *.scratch
    # - - - - - - write number of lst files to Job00*.err:
    jobfilerr=`echo $jrunnr $jobnumb | awk '{printf("Job%06d_%d.err",$1,$2)}'`
    # mv $jobfilerr $jobfilerr.org # keep original in case of errors.
    ls -l DAT*.lst | wc | awk '{printf("%7d\n",$1)}' >> $jobfilerr
    echo "               ... number of files written to $jobfilerr" 
    # - - - - - - reduce protocol lines in Job00*.out: 
    jobfilout=`echo $jrunnr $jobnumb | awk '{printf("Job%06d_%d.out",$1,$2)}'`
    wc $jobfilout > joboutwc.txt
    joblines=`cat joboutwc.txt | awk '{printf("%7d\n",$1)}'`
    head $jobfilout > joboutwc.txt
    /bin/mv joboutwc.txt $jobfilout
    echo "               ... append simulation infos to $jobfilout" 
    # - - - - - - jobhc3 script and steering file to Job00*.out:
    echo "______ jobhc3-$jrunnr ______" >> $jobfilout
    cat jobhc3-$jrunnr >> $jobfilout
    echo "______ parallel-$jrunnr ______" >> $jobfilout
    cat parallel-$jrunnr >> $jobfilout
    # - - - - - - write sum (in GBy) of DAT files to Job00*.out:  
    ./summ DAT >> $jobfilout
    echo "               ... writing total Gigabytes to $jobfilout"
    # - - - - - - time.txt (create 4 more info lines):
    timeinpt=`echo "totaltime.input" | awk '{ printf("%s",$1) }'`
    timefile=`echo "time.txt" | awk '{ printf("%s",$1) }'`
    linestxt=`wc $timefile | awk '{ printf("%d",$1) }'`
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if [ $linestxt -le 2 ] ; then
      # - - - - - - add number of lst files as 3rd line to file time.txt:
      tail -1 $jobfilerr >> $timefile
      # - - - - - - totaltime.input must be empty!
      if [ -e $timeinpt ] ; then
         rm $timeinpt ;
      fi
      touch $timeinpt
      listfiles=`echo "DAT" | awk '{ printf("%s*.lst",$1) }'`
      echo "               ... sum up time infos of all DAT$jrunnr-*.lst files"
      # - - - - - - check line number of start time infos:
      infoline=`echo "111" | awk '{ printf("%d",$1) }'`
      /bin/ls -1 $listfiles | head -1 > jobfilist.txt
      joblist1=`cat jobfilist.txt | awk '{ printf("%s",$1) }'`
      grep -n " START OF RUN " $joblist1 > jobfilist.txt
      listfirst=`grep -n " START OF RUN " $joblist1 | awk '{ printf("%s",$1) }'`
      sed -i "/======/ s/:/ /" jobfilist.txt
      infoline=`cat jobfilist.txt | awk '{ printf("%d",$1) }'`
      let "infoline = infoline + 2"
      /bin/rm jobfilist.txt
      # - - - - - - write wall time infos to `totaltime.input`:
      for i in $( ls -1 $listfiles ) ; do
        /bin/ls -1 $i >> $timeinpt ;
        head -$infoline $i | tail -1 >> $timeinpt ;
        tail -3 $i >> $timeinpt ;
      done
      # - - - - - - totaltime.f must be compiled and linked:
      ifort -C totaltime.f -o totaltime
      ./totaltime < $timeinpt > totaltime.print
      cat $timefile
      rm $timeinpt
      # time.txt looks now like former ones after `postprocess.c`.
      # echo "               ...... if  TOTAL CPU TIME (days) =    0.000000"
      # echo "               ...... check line number of  PRESENT TIME"
      # echo "               ...... after  === START OF RUN ===  in the" 
      # echo "               ...... first file of DAT$jrunnr*.lst"
      # - - - - - - add runnr to name `time.txt`:
      newtname=`echo $timefile $jrunnr | awk '{ printf("%s%06d",$1,$2) }'`
      echo $newtname >> $timefile
      /bin/mv $timefile $newtname
      echo $newtname 
      # - - - - - - end-of regular application. 
    fi
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fi
  # end-of case both files Job00* exist.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit 
