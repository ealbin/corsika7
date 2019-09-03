#!/bin/bash
#
# sumprocessing-ll.sh:
# ====================
#     create 4 shell scripts for the given run number,
#     which than have to be in the existing csk00???? subdirectory.        
# ------------------------------------------------------------------------
# compilation for analysis program:
#     f77 -fbounds-check -m32 showanalyscc.f -o showanalyscc
#     gfortran -fbounds-check showanalyscc.f -o showanalyscc
#     ifort -C showanalyscc.f -o showanalyscc
# ------------------------------------------------------------------------
# usage: ./sumprocessing-ll.sh <runnr>
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
runnr=`echo $1`
if [ -n "$runnr" ]; then
  if [ $runnr -lt 0 ]; then
    echo "           ... run number argument out of range:  0 <= runnr <= 999999"
  elif [ $runnr -gt 999999 ]; then
    echo "           ... run number argument out of range:  0 <= runnr <= 999999"
  else
  # - - - - - - length of argument is > 0, but subdir must exist:
  cskdirect=`echo $runnr | awk '{printf("csk%06d/", $1)}'`
  if [ ! -e $cskdirect ] ; then
    echo "           ... subdirectory " $cskdirect " does not exist."
    exit
  fi 
  echo "           ... cd" $cskdirect 
  # - - - - - - create script of NKG averages summation of .lst files:
  sumlistnkg=`echo $runnr | awk '{printf("sumlistnkginfo.sh%06d",$1)}'`
  echo "#!/bin/bash" > $sumlistnkg
  echo "# " >> $sumlistnkg
  echo "# sum up all NKG averages of \`.lst\` files:" >> $sumlistnkg
  echo "# --------------------------------------------------------------------" >> $sumlistnkg
  input=`echo $runnr | awk '{printf("ls -1 DAT%06d-*.lst > sumlistnkginfo.i%06d",$1,$1)}'`
  echo " " "$input" >> $sumlistnkg
  output=`echo $runnr | awk '{printf("./sumlistnkginfo < sumlistnkginfo.i%06d > sumlistnkginfo.out%06d",$1,$1)}'`
  echo " " $output >> $sumlistnkg
  chmod +x $sumlistnkg
  /bin/mv $sumlistnkg $cskdirect
  /bin/cp sumlistnkginfo.f $cskdirect
  /bin/cp sumlistnkginfo $cskdirect 
  echo "           ......" $sumlistnkg 
  # - - - - - - create script of summation of .long files:
  sumfiles=`echo $runnr | awk '{printf("sumlongifiles.sh%06d",$1)}'`
  echo "#!/bin/bash" > $sumfiles
  echo "# " >> $sumfiles
  echo "# sum up all \`.long\` files:" >> $sumfiles
  echo "# -----------------------------------------------------------------" >> $sumfiles
  input=`echo $runnr | awk '{printf("ls -1 DAT%06d-*.long > sumlongifiles.i%06d",$1,$1)}'`
  echo " " "$input" >> $sumfiles
  output=`echo $runnr | awk '{printf("./sumlongifiles < sumlongifiles.i%06d > sumlongifiles.out%06d",$1,$1)}'`
  echo " " $output >> $sumfiles
  chmod +x $sumfiles
  /bin/mv $sumfiles $cskdirect
  /bin/cp sumlongifiles.f $cskdirect
  /bin/cp sumlongifiles $cskdirect
  echo "           ......" $sumfiles 
  # - - - - - - - - create analysis script: 
  alyscript=`echo $runnr | awk '{ printf("showanalyscc.sh%06d",$1) }'`
  echo "#!/bin/bash" > $alyscript
  echo "#" >> $alyscript
  echo "#" $alyscript >> $alyscript
  echo "# =====================" >> $alyscript
  echo "# sum up all particle data to histograms (ascii coded):" >> $alyscript
  echo "# -------------------------------------------------------------------" >> $alyscript
  inpfilist=`echo $runnr | awk '{ printf("showanalyscc.i%06d",$1) }'`
  for i in 0 1 2 3 4 5 6 7 8; do  
    datafiles=`echo $runnr $i | awk '{ printf("DAT%06d-%d*",$1,$2) }'`
    echo "  ls -1" $datafiles "| grep t -v | grep n -v > $inpfilist-$i" >> $alyscript
  done
  echo "  cat $inpfilist-* > $inpfilist" >> $alyscript
  outscript=`echo $runnr | awk '{ printf("showanalyscc.out%06d",$1) }'`
  echo "  ./showanalyscc <" $inpfilist ">" $outscript >> $alyscript
  writfort9=`echo $runnr | awk '{ printf("showanalyscc.fort%06d",$1) }'`
  echo "  mv fort.9" $writfort9 >> $alyscript
  chmod +x $alyscript
  /bin/mv $alyscript $cskdirect
  # - - - - - - - - create llsubmit script for the shower analysis: 
  runscript=`echo $runnr | awk '{ printf("showanalyrun.sh%06d",$1) }'`
  echo "#!/bin/bash" > $runscript
  echo "# @ job_type = serial" >> $runscript
  echo "# @ initialdir = /fzk/cgwork/joe/corsika.conex/run/$cskdirect" >> $runscript
  echo "# @ environment = COPY_ALL;TMPDIR=/tmp" >> $runscript
  echo "# @ restart = no" >> $runscript
  echo "# @ requirements = (Arch == \"x86_64\") && (Feature == \"penryn\")" >> $runscript
  echo "# @ class = medium" >> $runscript
  echo "# @ wall_clock_limit = 140:00:00" >> $runscript
  echo "# @ error = showanalyrun.\$(cluster).err" >> $runscript
  echo "# @ output = showanalyrun.\$(cluster).out" >> $runscript
  echo "# @ queue" >> $runscript
  echo "#" >> $runscript
  echo " " ./$alyscript >> $runscript 
  echo "#" >> $runscript
  chmod +x $runscript
  /bin/mv $runscript $cskdirect
  /bin/cp showanalyscc.f $cskdirect 
  /bin/cp showanalyscc $cskdirect
  echo "           ...... llsubmit" $runscript "(using $alyscript)"
  # - - - - - - - - - - - - - - - - - - - -
  fi
#
# - - - - - - length of argument = 0 (missing argument):
else
  echo "           ... missing run number argument:  0 <= runnr <= 999999"
fi
#
