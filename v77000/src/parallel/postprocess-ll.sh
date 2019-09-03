#!/bin/bash
#
# postprocess-ll.sh:
# ==================
#     for the given run number it creates a shell script file to write a
#     protocol file named `Job00nnnn_scc.out`, which is necessary to get
#     a tabular of parallel simulations (Load Leveler system at the scc
#     opus cluster, iwrcgvor2.fzk.de);
#     parallel run numbers taken up to 009999.
# ------------------------------------------------------------------------
# usage: ./postprocess-ll.sh <runnr>
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
runnr=`echo $1`
if [ -n "$runnr" ]; then
  # - - - - - - length of argument is > 0:
  cskdirect=`echo $runnr | awk '{printf("csk%06d/",$1)}'`
  if [ ! -e $cskdirect ] ; then
      echo "           ... check runnr: " $cskdirect " does not exist."
      exit
  fi
  # - - - - - - create postjob script (run in subdir):
  postjob=`echo $runnr | awk '{printf("postjobscc-%06d.sh",$1)}'`
  echo "#!/bin/bash" > $postjob
  echo "# " >> $postjob
  echo "#" $postjob >> $postjob
  echo "# ====================" >> $postjob
  echo "# copy auxiliary files, create Job..._scc.err and Job..._scc.out:" >> $postjob
  echo "# -------------------------------------------------------------------" >> $postjob
  echo "  /bin/cp ../sum* ." >> $postjob
  echo "  /bin/cp ../totaltime* ." >> $postjob
  parallel=`echo $runnr | awk '{printf("parallel-%06d",$1)}'`
  echo "# - - - - creating job information files:" >> $postjob
  jobawk="awk '{printf(\"%7d\\n\",\$1)}'" 
  joberr=`echo $runnr | awk '{printf("> Job%06d_scc.err",$1)}'`
  echo "  wc job-file |" $jobawk $joberr >> $postjob
  jhead=`echo $runnr | awk '{printf("    head -1 job-file > Job%06d_scc.out",$1)}'`
  echo " " $jhead >> $postjob
  jtail=`echo $runnr | awk '{printf("    tail -1 job-file >> Job%06d_scc.out",$1)}'`
  echo " " $jtail >> $postjob
  jtime=`echo $runnr | awk '{printf("  ./totaltime >> Job%06d_scc.out",$1)}'`
  echo " " $jtime >> $postjob
  jpara=`echo $runnr | awk '{printf("  /bin/cat parallel-%06d >> Job%06d_scc.out",$1,$1)}'`
  echo " " $jpara >> $postjob
  kfiles=`echo $runnr $jobawk | awk '{printf("grep Files Job%06d_scc.out | %s %s%s",$1,$2,$3,$4)}'`
  echo "  nfiles=\`$kfiles\`" >> $postjob
  echo "           ... cd" $cskdirect
  echo "           ......" $postjob
  echo "if [ \$nfiles -gt 0 ]; then" >> $postjob
  echo "# - - - - number of files for /bin/ls command > 0:" >> $postjob
  # - - - - - - sum up < 1000 files:
  echo "  if [ \$nfiles -lt 1000 ]; then" >> $postjob
  echo "    # - - - - - - sum < 1000 files:" >> $postjob
  jsumm=`echo $runnr | awk '{printf("  ./summ DAT%06d- >> Job%06d_scc.out",$1,$1)}'`
  echo "   " $jsumm >> $postjob
  jobcnt=`echo $runnr | awk '{printf("/bin/ls -1 DAT%06d-* | grep t -v | grep n -v | wc",$1)}'`
  echo "   " $jobcnt "|" $jobawk $joberr >> $postjob
  echo "  else" >> $postjob
  # - - - - - - sum up >= 1000 files:
  echo "    # - - - - - - sum >= 1000 files:" >> $postjob
  datnam=`echo $runnr | awk '{printf("DAT%06d",$1)}'`
  echo "    jfiles=0" >> $postjob
  echo "    sumfil=0" >> $postjob
  echo "    for i in 0 1 2 3 4 5 6 7 8; do" >> $postjob
  echo "       ifiles=\`/bin/ls -l $datnam-\$i* | grep t -v | grep n -v | wc -l\`" >> $postjob
  echo "       echo \"               $datnam-\$i* files: \" \$ifiles" >> $postjob 
  echo "       du -k $datnam-\$i* | grep t -v | grep n -v > $datnam-\$i-kbytes" >> $postjob
  echo "       let \"jfiles = jfiles + ifiles\"" >> $postjob  
  echo "    done" >> $postjob
  echo "    echo \"                              sum:\" \$jfiles" >> $postjob
  echo "    cat $datnam-*-kbytes > summe.j" >> $postjob
  jsumm=`echo $runnr | awk '{printf("  ./summe < summe.j >> Job%06d_scc.out",$1)}'`
  echo "   " $jsumm >> $postjob
  echo "  fi" >> $postjob
  echo "else" >> $postjob
  echo "# - - - - check directly fortran code totaltime:" >> $postjob
  echo "  echo \"    incorrect time entries by totaltime process: 0. Files\"" >> $postjob
  echo "fi" >> $postjob
  echo "# " >> $postjob
  chmod +x $postjob
  /bin/mv $postjob $cskdirect
else
  # - - - - - - length of argument = 0 (i.e. missing argument):
  echo "           ... missing run number argument:  0 <= runnr <= 999999"
fi
