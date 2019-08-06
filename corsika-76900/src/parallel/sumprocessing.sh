#!/bin/bash
#           sumprocessing.sh:
# ===========================
#           create 3 shell scripts for the given run number
#           in the existing directory or subdirectory csk00????/,
#           to be used on the hc3 applying MPI parallelization.         
# ------------------------------------------------------------------------
# usage:    ./sumprocessing.sh <runnr>
# ------------------------------------------------------------------------
#        f77 -fbounds-check showanalyhc3.f -o showanalyhc3
#        gfortran -fbounds-check showanalyhc3.f -o showanalyhc3
#        ifort -C -check bounds showanalyhc3.f -o showanalyhc3
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
  # - - - - - - create script of NKG averages summation of .lst files:
  summinfo=`echo $runnr | awk '{printf("sumlistnkginfo.sh%06d",$1)}'`
  echo "#!/bin/bash" > $summinfo
  echo "# " >> $summinfo
  echo "# job_submit -p1 -cp -t30 -m2000 $summinfo" >> $summinfo
  echo "# " >> $summinfo
  echo "# sum up all NKG averages of \`.lst\` files:" >> $summinfo
  echo "# --------------------------------------------------------------------" >> $summinfo
  echo "           ... creating  $summinfo" 
  input=`echo $runnr | awk '{printf("ls -1 DAT%06d-*.lst > sumlistnkginfo.i%06d",$1,$1)}'`
  echo " " "$input" >> $summinfo
  output=`echo $runnr | awk '{printf("./sumlistnkginfo < sumlistnkginfo.i%06d > sumlistnkginfo.out%06d",$1,$1)}'`
  echo " " $output >> $summinfo
  chmod +x $summinfo
  # - - - - - - create script of summation of .long files:
  sumfiles=`echo $runnr | awk '{printf("sumlongifiles.sh%06d",$1)}'`
  echo "#!/bin/bash" > $sumfiles
  echo "# " >> $sumfiles
  echo "# job_submit -p1 -cp -t30 -m2000 $sumfiles" >> $sumfiles
  echo "# " >> $sumfiles
  echo "# sum up all \`.long\` files:" >> $sumfiles
  echo "# -----------------------------------------------------------------" >> $sumfiles
  echo "           ... creating  $sumfiles" 
  input=`echo $runnr | awk '{printf("ls -1 DAT%06d-*.long > sumlongifiles.i%06d",$1,$1)}'`
  echo " " "$input" >> $sumfiles
  output=`echo $runnr | awk '{printf("./sumlongifiles < sumlongifiles.i%06d > sumlongifiles.out%06d",$1,$1)}'`
  echo " " $output >> $sumfiles
  chmod +x $sumfiles
  # - - - - - - create histogram summation script:
  showhist=`echo $runnr | awk '{printf("showanalyhc3.sh%06d",$1)}'`
  echo "#!/bin/bash" > $showhist
  echo "# " >> $showhist
  echo "# job_submit -p1 -cp -t660 -m2000 $showhist" >> $showhist
  echo "# " >> $showhist
  echo "# sum up all particle data files to histograms (ascii coded):" >> $showhist
  echo "# ----------------------------------------------------------------" >> $showhist
  echo "           ... creating  $showhist"
  input1=`echo $runnr | awk '{printf("  if [ -e DAT%06d-000001 ] ; then",$1)}'`
  echo " " "$input1" >> $showhist
  input2=`echo $runnr | awk '{printf("ls -1 DAT%06d-?????? > showanalyhc3.i%06d",$1,$1)}'`
  echo " " "$input2" >> $showhist
  echo "  else" >> $showhist
  input=`echo $runnr | awk '{printf("ls -1 DAT%06d-* | grep t -v | grep n -v > showanalyhc3.i%06d",$1,$1)}'`
  echo " " "$input" >> $showhist
  echo "  fi" >> $showhist
  output=`echo $runnr | awk '{printf("./showanalyhc3 < showanalyhc3.i%06d > showanalyhc3.out%06d",$1,$1)}'`
  echo " " $output >> $showhist
  mvfort9=`echo $runnr | awk '{printf("mv fort.9 showanalyhc3.fort%06d",$1)}'`
  echo " " $mvfort9 >> $showhist
  chmod +x $showhist
  fi
#
else
  echo "           ... missing run number argument:  0 <= runnr <= 999999"
fi
#
