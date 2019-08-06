#!/bin/bash
#
# postprocesscray.sh:
# ===================
#        new postprocessing for the copied parallel corsika simulations
#        from CRAY after july 2013 (long files are automatically summed)
# ------------------------------------------------------------------------
#   (p0) one or two MPI job protocol file(s) must exist in the upper
#        (i.e. the hc3 submit) directory or already in the subdirectory;
#        after many lines containing `SUCCESSFULLY FINALIZED` the last
#        two lines of Job00jklm_359836.out must be:
#   =============================================================
#   Job 359836 completed at 31.07.2015/11:46:06 (COMPLETED)
#        These two lines must be appended by hand to the original CRAY
#        protocol file named like CSK_2013_10^20_4800.o385750;
#   (p2) time.txt must exist and contain 2 lines like
#       START TIME          STOP TIME       TIME (min)
#   1466262084.316990   1466263487.694263   423.389621
#   (p3) postprocesscray.sh must exist in this path (copied by hand);
#        first of all utilities scripts, executables and fortran source
#        codes will be copied first to this CRAY simulation path;
# ------------------------------------------------------------------------
#        cd csk00jklm/
# usage: ./postprocesscray.sh 
# ------------------------------------------------------------------------
#        data will be written to Job00jklm_%j.[err,out] and time.txt00jklm
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
/bin/rm -f *.scratch*
computype="ikp"
if [ $computype = "ikp" ] ; then
  /bin/cp "time.txt" "time.orig"
  /bin/cp ../totaltimenew* .
  /bin/cp ../readcsk2asci* .
  /bin/cp ../readmthinprtcls* .
  /bin/cp ../sortaugerhit* .
  /bin/cp ../summ* .
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - get run number from current directory name:
pwd > jobfileerr.tmp
sed -i "/csk/ s/csk/   /" jobfileerr.tmp
jrunnr=`cat jobfileerr.tmp | awk '{printf("%06d\n",$2)}'`
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - test on existence of file time.txt:
/bin/rm jobfileerr.tmp
if [ ! -e time.txt ] ; then
  echo "           ... postprocessing cannot be executed, because \`time.txt\` is missing"
  echo "           ...... or postprocessing already done for runnr" $jrunnr "or last lines"
  echo"            ...... ==============================================================="
  echo "           ...... Job 163862 completed at 30.05.2015/09:45:51 (COMPLETED) missing" 
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
fi
jdatnames=`wc "corsika_timetable" | awk '{ printf("%s",$1) }'`
let "jdatnames = $jdatnames / 2" 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - check on already available jobfiles:
jobrunnr=`echo $jrunnr | awk '{printf("Job%06d_*",$1)}'`
# - - - - - - njobfiles must be 2:
njobfiles=`ls -l $jobrunnr 2> jobfileerr.tmp | wc | awk '{printf("%7d\n",$1)}'`
/bin/rm jobfileerr.tmp
if [ $njobfiles -eq 2 ] ; then
  echo "               ... jobfiles $jobrunnr already available in this path."
else
  # forgot to copy jobfiles then first create them:
  jobfileerr=`echo $jrunnr | awk '{printf("Job%06d_2345670.err",$1)}'`
  jobfileout=`echo $jrunnr | awk '{printf("Job%06d_2345670.out",$1)}'`
  touch $jobfileerr
  /bin/cp "moabsub*" $jobfileout
fi
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - get current job-id from error protocol file [Job,job]*.e*:
if [ $computype = "ikp" ] ; then
  /bin/ls -1 Job$jrunnr*.err > jobfileerr.tmp
  sed -i "/Job/ s/_/ /" jobfileerr.tmp
  sed -i "/Job/ s/\./ /" jobfileerr.tmp
  sed -i "/Job/ s/err//" jobfileerr.tmp
  sed -i "/Job/ s/Job//" jobfileerr.tmp
fi
jobidnumber=`cat jobfileerr.tmp | awk '{printf("%d",$2)}'`
/bin/rm jobfileerr.tmp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - check Job00*.out if postprocessing was already done:
jobfileout=`echo $jrunnr $jobidnumber | awk '{printf("Job%06d_%d.out",$1,$2)}'`
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - write original and real number of files to Job00*.err:
jobfileerr=`echo $jrunnr $jobidnumber | awk '{printf("Job%06d_%d.err",$1,$2)}'`
du -k DAT$jrunnr-* | ./summe > summe.out
lstfiles=`echo $jrunnr | awk '{printf("DAT%06d-*.lst",$1)}'`
/bin/ls -1 $lstfiles > joblstfiles.txt
jlstnames=`wc joblstfiles.txt | awk '{printf("%7d\n",$1)}'`
/bin/rm joblstfiles.txt
datfiles=`echo $jrunnr | awk '{printf("DAT%06d-?????????-?????????",$1)}'`
/bin/ls -1 $datfiles > joboutempty.txt 2> /dev/null 
jdatempty=`wc joboutempty.txt | awk '{printf("%7d\n",$1)}'` 
if [ -e DAT$jrunnr-000001 ] ; then
   ls -1 DAT$jrunnr-?????? > jobfileerr.tmp
   jpartfiles=`wc jobfileerr.tmp | awk '{printf("%7d\n",$1)}'`
   /bin/rm jobfileerr.tmp
else
   jpartfiles=`echo jdatnames | awk '{printf("%7d\n",$1)}'`
fi
echo $jdatnames $jpartfiles | awk '{printf("%7d%7d\n",$1,$2)}' > $jobfileerr
echo "               ... writing number of files to $jobfileerr: " $jdatnames
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
# - - - - - - delete almost all lst files (keep max 5 in this path):
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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - reduce lines in Job00*.out, copy script and steering file: 
head -5 $jobfileout > joboutprotc.tmp 
tail -1 $jobfileout >> joboutprotc.tmp
echo "==========================================================================" >> joboutprotc.tmp
/bin/mv joboutprotc.tmp $jobfileout
echo "               ... append simulation infos to $jobfileout" 
echo "______ jobhc3-$jrunnr ______" >> $jobfileout
cat jobhc3-$jrunnr >> $jobfileout
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
