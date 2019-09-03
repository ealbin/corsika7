#!/bin/bash
#
# preprocess-ll.sh:
# =================
# prepare a parallel corsika simulation by scripts at the opus cluster 
# under `iwrcgvor2.fzk.de1 with Load Leveler commands (LL or ll);
# a corsika input steering file `parallel-002345` and a (renamed)
# corsika executable `corsika72495_stnd_QGSII_gheisha_pll` must exist;
# directories `DIRECT csk002345/`, `DIRECT ./`, `DIRECT ' '` are valid.
# -------------------------------------------------------------------------
# usage:
#    ./preprocess-ll.sh parallel-002345 corsika72495_stnd_QGSII_gheisha_pll
# start simulation:
#    cd csk002345/
#    ./jobscc-2345 &   
# -------------------------------------------------------------------------
#                                        juergen.oehlschlaeger@kit.edu
# -------------------------------------------------------------------------
  if [ $# -eq 2 ]; then
# - - - - - - copy script arguments:
    jobprepair=`echo $0`
    pllsteer=`echo $1`
    if [ ! -e $pllsteer ]; then
       /bin/ls -l $pllsteer
       exit 1;
    fi
    corsexec=`echo $2`
    if [ ! -e $corsexec ]; then
       /bin/ls -l $corsexec
       exit 1;
    fi
# - - - - - - get run number from file name pllsteer (argument #1):
    echo $1 > parallel-nnnnnn.tmp
    sed -i "/parallel/ s/-/ /" parallel-nnnnnn.tmp
    parun=`grep "parallel" parallel-nnnnnn.tmp | awk '{ printf("%d", $2) }'`
    rm parallel-nnnnnn.tmp 
# - - - - - - extract current corsika run number and create file:
    runnr=`grep RUNNR $pllsteer | awk '{ printf("%d", $2) }'`
    if [ $runnr -ne $parun ]; then
       sed -i "/RUNNR/ s/$runnr/$parun/" $pllsteer 
       let 'runnr = parun'
    fi
# - - - - - - name of shell script to submit the simulation:
    jobscript=`echo $runnr | awk '{ printf("jobscc-%06d", $1) }'`
# - - - - - - create name of corsika data path: 
    cskdirect=`echo $runnr | awk '{ printf("csk%06d/", $1) }'`
    echo "           ... creating directory" $cskdirect "and scripts;"
    echo "#!/bin/bash" > $jobscript
    echo "#" >> $jobscript
    echo "#" $jobscript >> $jobscript
    echo "# =============" >> $jobscript
    echo "# CORSIKA using parallel scripts at scc on iwrcgvor2 (load leveler)" >> $jobscript
    echo "#" >> $jobscript
    echo "  ./corsika-initial.ll.sh" $pllsteer $corsexec "LL" >> $jobscript 
    echo "#" >> $jobscript
    chmod +x $jobscript
# - - - - - - write runnr to current directory name in $pllsteer:
    directpll=`grep 'DIRECT' $pllsteer | awk '{print $2}'`
    if [ $directpll != "./" ]; then
       sed -i "/DIRECT/ s/1/0/g" $pllsteer
       sed -i "/DIRECT/ s/2/0/g" $pllsteer
       sed -i "/DIRECT/ s/3/0/g" $pllsteer
       sed -i "/DIRECT/ s/4/0/g" $pllsteer
       sed -i "/DIRECT/ s/5/0/g" $pllsteer
       sed -i "/DIRECT/ s/6/0/g" $pllsteer
       sed -i "/DIRECT/ s/7/0/g" $pllsteer
       sed -i "/DIRECT/ s/8/0/g" $pllsteer
       sed -i "/DIRECT/ s/9/0/g" $pllsteer
       sed -i s/csk000000/\./ $pllsteer
    fi
# - - - - - - create or clear corsika sub path:
    if [ ! -e $cskdirect ] ; then
       /bin/mkdir $cskdirect
    else
       /bin/rm -f $cskdirect\*
    fi
# - - - - - - copy script files to corsika subdirectory:
    cp -p ../src/parallel/corsika-*.ll.sh $cskdirect
    cp -p $jobscript $cskdirect
    cp -p $corsexec $cskdirect
    cp -p $pllsteer $cskdirect
# - - - - - - switch to subdirectory and create/clear epos path: 
    if [ -e $cskdirect ] ; then
    cd $cskdirect
    if [ ! -e epos/ ] ; then
      /bin/mkdir epos/
    else
      /bin/rm -f epos/*
    fi 
#   /bin/cp ../epos/* epos/ 
    /bin/cp ../atmabs.dat .
    /bin/cp ../EGSDAT6_.05 .
    /bin/cp ../EGSDAT6_1. .
    /bin/cp ../EGSDAT6_.15 .
    /bin/cp ../EGSDAT6_.25 .
    /bin/cp ../EGSDAT6_3. .
    /bin/cp ../EGSDAT6_.4 .
    /bin/cp ../elasct.bin .
    /bin/cp ../GLAUBTAR.DAT .
    /bin/cp ../gxsect.bin .
    /bin/cp ../mirreff.dat .
    echo "           ... copying all corsika model-dependent tabular files;"
    /bin/cp ../nuclear.bin .
    /bin/cp ../NUCLEAR.BIN .
    /bin/cp ../NUCNUCCS .
    /bin/cp ../nunstab.data .
    /bin/cp ../QGSDAT01 .
    /bin/cp ../qgsdat-II-04 .
    /bin/cp ../quanteff.dat .
    /bin/cp ../SECTNU .
    /bin/cp ../sectnu-II-04 .
    /bin/cp ../sigmapi.bin .
    /bin/cp ../tables.dat .
    /bin/cp ../UrQMD-1.3.1-xs.dat .
    /bin/cp ../VENUSDAT .
    cd ..
    fi
# - - - - - - end of working in corsika subdirectory.
    echo "           ...... then switch to $cskdirect and execute ./$jobscript ......"
# - - - - - - number of arguments < 2:
  else
    if [ $# -eq 1 ]; then
       echo "           ... missing second argument: corsika-executable"
       pllsteer=`echo $1`
       if [ ! -e $pllsteer ]; then
          /bin/ls -l $pllsteer
          exit 1;
       fi
    else   
       echo "           ... missing two arguments: parallel-00ijkl corsika-executable"
    fi
  fi
exit
