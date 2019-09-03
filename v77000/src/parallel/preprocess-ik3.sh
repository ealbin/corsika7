#!/bin/bash
#
# preprocess-ik3.sh:
# ==================
#    prepare a parallel corsika simulation by scripts at the ik3 computing
#    cluster; a corsika steering file `parallel-003210` and the (renamed)
#    corsika executable `corsika75061_stnd_QGSII4_gheisha_pll` must exist;
#    directory name `DIRECT csk003210/` is valid if you start the 
#    simulation in the corsika main path, whereas `DIRECT ./`, `DIRECT ' '`
#    are valid starting the simulation from the subdirectoy csk003210.
#    Check names jobik3-..... of submit scripts.
# -------------------------------------------------------------------------
# usage:
# ./preprocess-ik3.sh parallel-003210 corsika76400_stnd_QGSII4_gheiatm_pll
# then:
# cd csk003210/
# then:
# ./corsika-initial-ik3.sh parallel-003210 corsika76400_stnd_QGSII4_gheiatm_pll ik3
# -------------------------------------------------------------------------
#                                        juergen.oehlschlaeger@kit.edu
# -------------------------------------------------------------------------
if [ $# -eq 2 ]; then
# - - - - - - copy script arguments:
  jobprepair=`echo $0`
  parallinp=`echo $1`
  if [ ! -e $parallinp ]; then
     ls -l $parallinp
     exit 1;
  fi
  corsexec=`echo $2`
  if [ ! -e $corsexec ]; then
     ls -l $corsexec
     exit 1;
  fi
# - - - - - - get run number from file name parallinp (argument #1):
  echo $1 > parallel-nnnnnn.tmp
  sed -i "/parallel/ s/-/ /" parallel-nnnnnn.tmp
  parun=`grep "parallel" parallel-nnnnnn.tmp | awk '{ printf("%d",$2)}'`
  rm parallel-nnnnnn.tmp 
# - - - - - - make sure that cut-files will be written and kept:
  sed -i "/PARALLEL/ s/ F/ T/" $parallinp
# - - - - - - extract current corsika run number and create files:
  runnr=`grep RUNNR $parallinp | awk '{ printf("%d",$2)}'`
  if [ $runnr -ne $parun ]; then
     sed -i "/RUNNR/ s/$runnr/$parun/" $parallinp 
     let 'runnr = parun'
  fi
  jrunnr=`echo $runnr | awk '{ printf("%06d",$1)}'`
  jobscript=`echo $runnr | awk '{ printf("jobik3-%06d",$1)}'`
# - - - - - - create name of corsika data path: 
  cskdirect=`echo $runnr | awk '{ printf("csk%06d/",$1)}'`
  mkdir -p $cskdirect
  echo "           ... creating directory" $cskdirect "and copying utilities;"
  echo "#!/bin/bash" > $jobscript
  echo "#" >> $jobscript
  echo "#" $jobscript >> $jobscript
  echo "# =============" >> $jobscript
  echo "# CORSIKA simulation using parallel scripts at ik3 cluster;" >> $jobscript
  echo "# parallel-$jrunnr and jobik3-$jrunnr must exist." >> $jobscript 
  echo "# ---------------------------------------------------------" >> $jobscript
  echo "#" >> $jobscript
  echo "# (a) corsika.trunk/run/" >> $jobscript
  echo "# preprocess-ik3.sh" $parallinp $corsexec >> $jobscript
  echo "#" >> $jobscript
  echo "# (b) corsika.trunk/run/csk$jrunnr/" >> $jobscript
  echo "./corsika-initial-ik3.sh" $parallinp $corsexec "ik3" >> $jobscript 
  chmod +x $jobscript
  if [ -e "SIM$jrunnr.reas" ] ; then
     cp -p "SIM$jrunnr.reas" $cskdirect
     cp -p "SIM$jrunnr.list" $cskdirect
     cp -p $parallinp "SIM$jrunnr.inp"
     cp -p "SIM$jrunnr.inp" $cskdirect
     cp -p $jobscript "SIM$jrunnr.sh"
     cp -p "SIM$jrunnr.sh" $cskdirect
  fi
# - - - - - - write runnr to current directory name in $parallinp:
  directpll=`grep 'DIRECT' $parallinp | awk '{print $2}'`
  if [ $directpll != "./" ]; then
     sed -i "/DIRECT/ s/1/0/g" $parallinp
     sed -i "/DIRECT/ s/2/0/g" $parallinp
     sed -i "/DIRECT/ s/3/0/g" $parallinp
     sed -i "/DIRECT/ s/4/0/g" $parallinp
     sed -i "/DIRECT/ s/5/0/g" $parallinp
     sed -i "/DIRECT/ s/6/0/g" $parallinp
     sed -i "/DIRECT/ s/7/0/g" $parallinp
     sed -i "/DIRECT/ s/8/0/g" $parallinp
     sed -i "/DIRECT/ s/9/0/g" $parallinp
     sed -i s/csk000000/\./ $parallinp
  fi
# - - - - - - clear and create subdirectory and copy scripts:
  cp -p corsika-*ik3.sh $cskdirect
  cp -p $jobscript $cskdirect
  cp -p $corsexec $cskdirect
  cp -p $parallinp $cskdirect
# - - - - - - switch to subdirectory and optionally create epos path: 
  cd $cskdirect
  rm -rf epos
  mkdir epos/
  cp ../../epos/* epos/ 
  cp ../atmprof*.dat .
  cp ../brems_fin.bin .
  cp ../cohff.bin .
  cp ../dpmCT14LL.pds .
  cp ../dpmjpar.dat .
  cp ../EGSDAT6_.05 .
  cp ../EGSDAT6_1. .
  cp ../EGSDAT6_.15 .
  cp ../EGSDAT6_.25 .
  cp ../EGSDAT6_3. .
  cp ../EGSDAT6_.4 .
  cp ../elasct.bin .
  cp ../fluodt.dat .
  cp ../glaub*.glb .
  cp ../gxsect.bin .
  cp ../mirreff.dat .
  cp ../neuxsc-ind_260.bin .
  cp ../nuclear.bin .
  cp ../nunstab.data .
  cp ../NUCNUCCS .
  echo "           ... copying all corsika model-dependent tabular files;"
  cp ../QGSDAT01 .
  cp ../qgsdat-II-04 .
  cp ../quanteff.dat .
  cp ../SECTNU .
  cp ../sectnu-II-04 .
  cp ../sigmapi.bin .
  cp ../UrQMD-1.3.1-xs.dat .
  cp ../VENUSDAT .
# copy utilities:
  cp ../readcsk2* .
  cp ../readpart* .
  cp ../readmthi* .
  cp ../suml* .
  cp ../summ* .
  cp ../concatcorsika* .
  cp ../totaljobfile* .
  cp ../totaltimeik3* .
  cp ../postprocess-ik3.* .
  cp ../sumrawcoreas* .
  cd ..
# - - - - - - end of working in corsika subdirectory.
  echo "           ... cd $cskdirect"
  echo "           ... execute following line to start simulation:"
  tail -1 $jobscript
# - - - - - - number of arguments < 2:
else
  if [ $# -eq 1 ]; then
     echo "           ... missing second argument: corsika-executable"
     parallinp=`echo $1`
     if [ ! -e $parallinp ]; then
        ls -l $parallinp
        exit 1;
     fi
  else   
     echo "           ... missing two arguments: parallel-0ijklm corsika-executable"
  fi
fi
exit
