#!/bin/bash
#
# postprocessnew.sh:
# ==================
#        new postprocessing for the hc3 parallel corsika simulations
#        after september 2016 - also valid for parallel `coreas` radio 
#        simulations - , DAT*.long files are automatically summed to
#        DAT*-999999999.long; no DAT*-...-...lst files are written 
#        to csk... subpath - exept in debug mode by `T` as: 
# job_submit -cp -p32 -t290 -m2000 -oJob001485_%j.out -eJob001485_%j.err
#   mpirun mpi_corsika75028_mthi_QGSII4_gheisha_runner parallel-001485 T
# ------------------------------------------------------------------------
#   (p0) one or two MPI job protocol file(s) must exist in the upper 
#        (i.e. submit) directory or already in the subdirectory; 
#        the last lines of Job00nnnn_%j.out must be: 
#        SUCCESSFULLY FINALIZED
#        =============================================================
#        Job 220836 completed at 08.12.2016/11:46:06 (COMPLETED)
#   (p1) jobhc3-00nnnn and parallel-00nnnn must exist in this path too.
#   (p2) time.txt must exist and contain two lines like
#       START TIME          STOP TIME       TIME (min)
#   1466262084.316990   1466263487.694263    23.389621
#   (p3) `postprocessnew.sh` must exist in this path.
#   (p4) script `summ`, executable `summe` and source code `summe.f`
#        must exist in this path.
#   (p5) executable/source code `totaltimenew*` must exist in this path.
# ------------------------------------------------------------------------
#        cd csk00nnnn/
# usage: ./postprocessnew.sh 
# ------------------------------------------------------------------------
#        data will be written to Job00nnnn_%j.[err,out] and time.txt00nnnn
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
/bin/rm -f *.scratch*
computype="ikp"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - get run number from current directory name:
pwd > jobfileerr.tmp
sed -i "/csk/ s/csk/   /" jobfileerr.tmp
jrunnr=`cat jobfileerr.tmp | awk '{printf("%06d\n",$2)}'`
if [ `cut -b 7-9 jobfileerr.tmp` = "fh1" ] ; then           # on $WORK
  computype="fh1"
else
  if [ `cut -b 2-5 jobfileerr.tmp` = "work" ] ; then
    computype="hc3"
  else
    if [ `cut -b 2-5 jobfileerr.tmp` = "home" ] ; then
      computype="hc3"
    else
      if [ `echo $HOSTNAME | cut -b 1-4` = "uc1n" ] ; then
        computype="bwc"
      else
        if [ `cut -b 2-5 jobfileerr.tmp` = "univ" ] ; then
          computype="stu"
        fi
      fi
    fi
  fi
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - test on existence of file time.txt:
/bin/rm jobfileerr.tmp
if [ ! -e time.txt ] ; then
  echo "           ... postprocessing cannot be executed, because \`time.txt\` is missing"
  echo "           ...... or postprocessing already done for runnr" $jrunnr "or last line"
  echo "           ...... Job 163862 completed at 04.05.2015/09:45:51 (COMPLETED) missing" 
  exit
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - test on syntax of time table files:
if [ -e corsika_timetable-000001 ] ; then
  # newer names of time table files:
  chtitle=`head -1 corsika_timetable-000001 | awk '{printf("%s\n",$1)}'` 
  if [ "$chtitle" == "#Columns" ] ; then
    # time table files with shell comment as first line:
    echo "               ... corsika time tables with comment line at the beginning."
    listfiles="corsika_timetable-*"
    /bin/rm -f jobfiletimes.tmp
    touch jobfiletimes.tmp 
    for filename in $( ls -1 $listfiles ) ; do
      jlines=`wc $filename | awk '{ printf("%s",$1) }'`
      let "jlines = $jlines - 1"
      if [ $jlines -gt 0 ] ; then
        tail -$jlines $filename >> jobfiletimes.tmp
      fi
    done
    mv jobfiletimes.tmp corsika_timetable
  else 
    # time table files without shell comments as first line:
    cat corsika_timetable-* > corsika_timetable
  fi
else
  if [ -e corsika_timetable-1 ] ; then
    # old names of time table files:
    cat corsika_timetable-* > corsika_timetable
  fi
fi
jdatnames=`wc "corsika_timetable" | awk '{ printf("%s",$1) }'`
let "jdatnames = $jdatnames / 2" 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - check on already available jobfiles:
if [ $computype = "hc3" ] ; then
  jobrunnr=`echo $jrunnr | awk '{printf("Job%06d_*",$1)}'`
else
  if [ $computype = "fh1" ] ; then
    jobrunnr=`echo $jrunnr | awk '{printf("jobfh1-%06d_*",$1)}'`
  else
    if [ $computype = "bwc" ] ; then
        jobrunnr=`echo $jrunnr | awk '{printf("job_uc1_jobwcl-%06d_*",$1)}'`
    else
      if [ $computype = "stu" ] ; then
        jobrunnr=`echo $jrunnr | awk '{printf("jobpll-%06d.[e,o]*",$1)}'`
      fi
    fi
  fi
fi
# - - - - - - njobfiles must be 2:
njobfiles=`ls -l $jobrunnr 2> jobfileerr.tmp | wc | awk '{printf("%7d\n",$1)}'`
/bin/rm jobfileerr.tmp
if [ $njobfiles -eq 2 ] ; then
  echo "               ... jobfiles $jobrunnr already available in this path."
else
  if [ $njobfiles -eq 1 ] ; then
    if [ $computype = "fh1" ] ; then
      echo "               ... file $jobrunnr already available in this path."
    else
      if [ $computype = "bwc" ] ; then
        echo "               ... file $jobrunnr already available in this path."
      else
        echo "               ... jobfile $jobrunnr still missing in this path."
        exit
      fi
    fi
  else
    /bin/mv ../$jobrunnr .
    ljobfiles=`ls -l $jobrunnr | wc | awk '{printf("%7d\n",$1)}'`
    if [ $ljobfiles -ne 2 ] ; then
      if [ $ljobfiles -eq 1 ] ; then
        if [ $computype = "fh1" ] ; then
          echo "               ... jobfile $jobrunnr moved from the upper directory."
        else
          if [ $computype = "bwc" ] ; then
            echo "               ... jobfile $jobrunnr moved from the upper directory."
          else
            echo "               ... jobfile $jobrunnr still missing in this path."
            exit
          fi
        fi
      else
        echo "               ... jobfile $jobrunnr still missing or several moves done,"
        echo "               ...... keep only newest jobfile (check file list)."
        exit
      fi
    else
      echo "               ... jobfiles $jobrunnr moved from the upper directory."
    fi
  fi
fi
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - get current job-id from error protocol file [Job,job]*.e*:
if [ $computype = "hc3" ] ; then
  /bin/ls -1 Job$jrunnr*.err > jobfileerr.tmp
  sed -i "/Job/ s/_/ /" jobfileerr.tmp
  sed -i "/Job/ s/\./ /" jobfileerr.tmp
  sed -i "/Job/ s/err//" jobfileerr.tmp
  sed -i "/Job/ s/Job//" jobfileerr.tmp
else
  if  [ $computype = "fh1" ] ; then
    jobrunout=`echo $jrunnr | awk '{printf("jobfh1-%06d_*.out",$1)}'`
    /bin/ls -1 $jobrunout > jobfileerr.tmp
    sed -i "/job/ s/_/ /" jobfileerr.tmp
    sed -i "/job/ s/\./ /" jobfileerr.tmp
  else
    if [ $computype = "bwc" ] ; then
      jobrunout=`echo $jrunnr | awk '{printf("job_uc1_jobwcl-%06d_*.out",$1)}'`
      /bin/ls -1 $jobrunout > jobfileerr.tmp
      sed -i "/job/ s/_/-/" jobfileerr.tmp
      sed -i "/job/ s/_/-/" jobfileerr.tmp
      sed -i "/job/ s/_/ /" jobfileerr.tmp
      sed -i "/job/ s/\./ /" jobfileerr.tmp
    else
      if  [ $computype = "stu" ] ; then
        /bin/ls -1 jobpll-$jrunnr.e* > jobfileerr.tmp
        sed -i "/job/ s/-/ /" jobfileerr.tmp
        sed -i "/job/ s/\./ /" jobfileerr.tmp
        sed -i "/job/ s/e/ /" jobfileerr.tmp
        sed -i "/job/ s/jobpll/ /" jobfileerr.tmp
      fi
    fi
  fi
fi
jobidnumber=`cat jobfileerr.tmp | awk '{printf("%d",$2)}'`
/bin/rm jobfileerr.tmp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - check Job00*.out if postprocessing was already done:
jobfileout=`echo $jrunnr $jobidnumber | awk '{printf("Job%06d_%d.out",$1,$2)}'`
if  [ $computype = "stu" ] ; then
  jobstutout=`echo $jrunnr $jobidnumber | awk '{printf("jobpll-%06d.o%s",$1,$2)}'`
  /bin/cp $jobstutout $jobfileout
else
  if [ $computype = "fh1" ] ; then
    jobfh1out=`echo $jrunnr $jobidnumber | awk '{printf("jobfh1-%06d_%d.out",$1,$2)}'`
    /bin/cp $jobfh1out $jobfileout
    echo "               ... jobfh1-..._$jobidnumber successfully finished (fh1)."
  else
    if [ $computype = "bwc" ] ; then
      touch $jobfileout
      echo "               ... job_uc1_..._$jobidnumber successfully completed (bwc)."
    else # other `$computype` specifications (hc3):
      tail -1 $jobfileout > jobfileerr.tmp
      jobcomplete=`cat jobfileerr.tmp | awk '{printf("%s",$3)}'`
      /bin/rm jobfileerr.tmp
      if [ $jobcomplete = "completed" ] ; then
        echo "               ... job $jobidnumber successfully completed (hc3)."
      else
        if [ $jobcomplete = "resources:" ] ; then
          echo "               ... job $jobidnumber successfully completed (stu)."
        else
          echo "           ... job $jobidnumber not completed or"
          echo "           ...... postprocessing seems already done for runnr" $jrunnr 
          echo "           ...... check content of file \`time.txt\`"
          exit
        fi
      fi
    fi
  fi
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - write original and real number of files to Job00*.err:
jobfileerr=`echo $jrunnr $jobidnumber | awk '{printf("Job%06d_%d.err",$1,$2)}'`
du -k | ./summe > summe.out
lstfiles=`echo $jrunnr | awk '{printf("DAT%06d-*.lst",$1)}'`
/bin/ls -1 $lstfiles 2> joblstfiles.err | wc > joblstfiles.txt
/bin/rm joblstfiles.err
jlstnames=`cat joblstfiles.txt | awk '{printf("%7d\n",$1)}'`
/bin/rm joblstfiles.txt
datfiles=`echo $jrunnr | awk '{printf("DAT%06d-?????????-?????????",$1)}'`
/bin/ls -1 $datfiles > joboutempty.txt 2> /dev/null 
jdatempty=`wc joboutempty.txt | awk '{printf("%7d\n",$1)}'` 
if [ -e DAT$jrunnr-000001 ] ; then
   ls -1 DAT$jrunnr-?????? > jobfileerr.tmp
   jpartfiles=`wc "jobfileerr.tmp" | awk '{printf("%7d\n",$1)}'`
   /bin/rm -f jobfileerr.tmp
else
   jpartfiles=`echo $jdatnames | awk '{printf("%7d\n",$1)}'`
fi
echo $jdatnames $jpartfiles | awk '{printf("%7d%7d\n",$1,$2)}' > $jobfileerr
echo "               ... writing number of files to $jobfileerr: " $jdatnames
# - - - - - - check >=SIM010000.*-files on antennas and gamma factor:
gammfact=`echo "0.0"`
gexpfact=`echo "0.0"`
jantenns=`echo "0"`
if [ $jrunnr -ge 10000 ] || [ $jrunnr -lt 100 ] ; then
  jantenns=`grep "AntennaPosition" "SIM$jrunnr.list" | wc | awk '{printf("%d",$1)}'`
  jgamfact=`grep "gamma" "SIM$jrunnr.list" | wc | awk '{printf("%d",$1)}'`
  if [ $jgamfact -gt 0 ] ; then
    head -1 "SIM$jrunnr.list" > jobgamfact.txt
    gammfact=`cat "jobgamfact.txt" | awk '{printf("%s",$8)}'`
    gexpfact=`cat "jobgamfact.txt" | awk '{printf("%s",$9)}'`
    /bin/rm -f jobgamfact.txt
  fi
  echo "               ... writing number of antennas to $jobfileerr: " $jantenns
fi
echo $jantenns $gammfact $gexpfact | awk '{printf("%7d%7.1f%7s\n",$1,$2,$3)}' >> $jobfileerr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - add third line to time.txt (necessary to run totaltimenew):
echo $jdatnames >> time.txt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - delete empty particle data files (if any):
if [ $jdatempty -gt 0 ] ; then
  for filename in $( cat joboutempty.txt ) ; do
    jdatfsize=`/bin/ls -1s $filename | awk '{printf("%d",$1)}'` 
    if [ $jdatfsize -eq 0 ] ; then
      # echo $filename ;
      /bin/rm -f $filename ;
    fi
  done
  echo "               ... deleting empty files $datfiles (if any)"
fi 
/bin/rm -f joboutempty.txt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - delete almost all lst files (keep max 5 in this path, if any):
if [ $jlstnames -gt 0 ] ; then
  /bin/ls -1 DAT$jrunnr-*.lst > erase-lstfiles-all.tmp
  nlst=`wc "erase-lstfiles-all.tmp" | awk '{printf("%d",$1)}'`
  let "ndel = $nlst - 5"
  if [ $ndel -gt 0 ] ; then
    cat erase-lstfiles-all.tmp | head -$ndel > erase-lstfiles-del.tmp
    for lstfile in $( cat "erase-lstfiles-del.tmp" ) ; do
      # echo $lstfile ;  
      /bin/rm -f $lstfile ;
    done
    /bin/rm -f erase-lstfiles-del.tmp
    echo "               ... deleting almost all DAT$jrunnr-?????????-?????????.lst"
  fi
  /bin/rm -f erase-lstfiles-all.tmp
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - reduce lines in Job00*.out, copy script and steering file: 
head -5 $jobfileout > joboutprotc.tmp 
tail -1 $jobfileout >> joboutprotc.tmp
echo "==========================================================================" >> joboutprotc.tmp
/bin/mv joboutprotc.tmp $jobfileout
echo "               ... append simulation infos to $jobfileout" 
if [ $computype = "hc3" ] ; then 
  echo "______ jobhc3-$jrunnr ______" >> $jobfileout
  cat jobhc3-$jrunnr >> $jobfileout
else 
  if [ $computype = "fh1" ] ; then
    echo "______ jobfh1-$jrunnr ______" >> $jobfileout
    cat jobfh1-$jrunnr >> $jobfileout
  fi
  if [ $computype = "bwc" ] ; then
    echo "______ jobwcl-$jrunnr ______" >> $jobfileout
    cat jobwcl-$jrunnr >> $jobfileout   
  fi
  if [ $computype = "stu" ] ; then
    echo "______ jobpll-$jrunnr ______" >> $jobfileout
    cat jobpll-$jrunnr >> $jobfileout
  fi
fi
echo "______ parallel-$jrunnr ______" >> $jobfileout
cat parallel-$jrunnr >> $jobfileout
echo "               ... writing total Gigabytes to $jobfileout"
./summ DAT >> $jobfileout
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - time.txt00nnnn (create 5 more info lines):
newtname=`echo $jrunnr | awk '{printf("time.txt%06d",$1)}'`
/bin/mv time.txt corsika_showertime
./totaltimenew > totaltimenew.out
echo $newtname >> totaltimenew.out 
/bin/mv totaltimenew.out $newtname
echo "               ... writing full time information to $newtname"
cat $newtname
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit 
