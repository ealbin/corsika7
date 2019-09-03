#!/bin/bash
#
# postprocess-ik3.sh:
# ===================
#     For the given run number it creates a shell script file to write a
#     protocol file named `Jobjrunnr_ik3.out`, which is necessary to get
#     a tabular of parallel simulations by applying showparallel.sh; 
#     this script should only be executed in the `cskjrunnr`
#     subdirectory, otherwise all necessary corsika tabular files
#     will be deleted in the wrong (working) path.
# ------------------------------------------------------------------------
# source codes sumlongifiles.f, sumlistnkginfo.f, totaljobfile.f, and
# totaltimeik3.f must exist in this path.
# ------------------------------------------------------------------------
# if error case:
#     Check run time messages in file `job-file` or `Jobjrunnr_ik3.out`
#     i.e. `runtime`, `error`, `exited`, therefore fortran program 
#     `totaltimeik3` cannot run.
# ------------------------------------------------------------------------
# usage: ./postprocess-ik3.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
pwd > postprocess.tmp
sed -i "/csk/ s/csk/   /" postprocess.tmp
jrunnr=`cat postprocess.tmp | awk '{printf("%06d\n",$2)}'`
timelist=`echo $jrunnr | awk '{printf("time.txt%06d",$1)}'`
# cskdirect=`echo $jrunnr | awk '{printf("csk%06d",$1)}'`
rm -f postprocess.tmp
if [ ! -e "$timelist" ]; then
# - - - delete and/or copy files in the subdirectory:
rm -f *.scratch
rm -f ../corsika-*ik3.sh.e*
rm -f ../corsika-*ik3.sh.o*
rm -f corsika-*ik3.sh
rm -f corsika*pll
rm -rf epos/
rm -f atmprof*.dat
rm -f fluodt.dat
rm -f glaub*
rm -f dpm*
rm -f *.bin
rm -f *.dat
echo "         ... deleting corsika tabular files ..."
rm -f [E,G,m,N,n,Q,q,U,V]*
rm -f SECTNU 
rm -f sectnu*
rm status_* # now not removed for checks.
rm -f DAT*.jobdone
rm -f DAT*.cut
rm -f DAT*.end
rm -f DAT*.inp
rm -f DAT*.fl*
# - - - check existence of Real-time protocol file `job-file`:
nfiles=`wc job-file | awk '{printf("%d",$1)}'`
if [ $nfiles -le 0 ]; then
  echo "          Warning: Real-time protocol file 'job-file' may be empty."
  exit
fi
# - - - clear important file `job-file` of all IEEE messages:
sed -i -e 's/ IEEE_DENORMAL / /' job-file 
sed -i -e 's/Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG / /' job-file 
cp job-file corsika_timetable
# - - - create `sumlongscript.shjrunnr` and run summation:
longfile=`echo $jrunnr | awk '{printf("DAT%06d-999989999.long",$1)}'`  
if [ -e "fort.48" ] ; then
  mv fort.48 $longfile
fi 
echo "         ... sum up 'DAT*.long' files to 'DAT$jrunnr-999989999.long' ..."  
sumlongscript=`echo $jrunnr | awk '{printf("sumlongifiles.sh%06d",$1)}'`
echo "#!/bin/bash" > $sumlongscript
echo "# " >> $sumlongscript
echo "# sum up all *.long files to file sumlongifiles.sum$jrunnr:" >> $sumlongscript
echo "# ------------------------------------------------------------" >> $sumlongscript
echo "#                           juergen.oehlschlaeger@kit.edu" >> $sumlongscript
echo "# ------------------------------------------------------------" >> $sumlongscript
echo "# gfortran -O0 -fbounds-check sumlongifiles.f -o sumlongifiles" >> $sumlongscript
echo "# " >> $sumlongscript
gfortran -O0 -fbounds-check sumlongifiles.f -o sumlongifiles
echo "  ls -1 DAT$jrunnr-*.long > sumlongifiles.i$jrunnr" >> $sumlongscript 
echo "  ./sumlongifiles < sumlongifiles.i$jrunnr > sumlongifiles.out$jrunnr" >> $sumlongscript
echo "  cp -p sumlongifiles.sum$jrunnr DAT$jrunnr-999989999.long" >> $sumlongscript
chmod +x $sumlongscript
./$sumlongscript
mkdir longfiles
longfile8=`echo $jrunnr | awk '{printf("DAT%06d-[0-8]*.long",$1)}'`
mv $longfile8 longfiles
# - - - create names of output files:
datname=`echo $jrunnr | awk '{printf("DAT%06d",$1)}'`
parallel=`echo $jrunnr | awk '{printf("parallel-%06d",$1)}'`
jobfout=`echo $jrunnr | awk '{printf("Job%06d_ik3.out",$1)}'`
jobferr=`echo $jrunnr | awk '{printf("Job%06d_ik3.err",$1)}'`
jobik=`echo $jrunnr | awk '{printf("jobik3-%06d",$1)}'`
# - - - pattern of names of `DAT*.lst-files:
lstfiles=`echo $jrunnr | awk '{printf("DAT%06d-*.lst",$1)}'`
ls -1 $lstfiles 2> joblstfiles.err | wc > joblstfiles.txt
rm joblstfiles.err
# - - - number of parts of the parallel corsika simulation:
jlstnames=`cat joblstfiles.txt | awk '{printf("%7d\n",$1)}'`
rm joblstfiles.txt
if [ $jlstnames -gt 0 ] ; then
  # - - - - sum up nkg averages from `DAT*.lst`-files:
  ls -1 $datname-*.lst > sumlistnkginfo.i$jrunnr
  echo "         ... sum up NKG averages from 'DAT*.lst'-files ..."
  gfortran -O0 -fbounds-check sumlistnkginfo.f -o sumlistnkginfo
  ./sumlistnkginfo < sumlistnkginfo.i$jrunnr > sumlistnkginfo.out$jrunnr
  # - - - - delete almost all lst files (keep max 3 in this path):
  ls -1 DAT$jrunnr-*.lst > erase-lstfiles-all.tmp
  nlst=`wc "erase-lstfiles-all.tmp" | awk '{printf("%d",$1)}'`
  let "ndel = $nlst - 3"
  if [ $ndel -gt 0 ] ; then
    cat erase-lstfiles-all.tmp | head -$ndel > erase-lstfiles-del.tmp
    for lstname in $( cat "erase-lstfiles-del.tmp" ) ; do
      # echo $lstname ;  
      rm -f $lstname ;
    done
    rm -f erase-lstfiles-del.tmp
    echo "         ... deleting almost all DAT$jrunnr-*.lst ..."
  fi
  rm -f erase-lstfiles-all.tmp
fi
# - - - - calculate wall time from file 'job-file':
echo "         ... calculate total wall time using file 'job-file' ..."
gfortran -O0 -fbounds-check totaljobfile.f -o totaljobfile
./totaljobfile > totaljobfile.out
secnds=`cat totaljobfile.out | awk '{printf("%d",$2)}'`
cp totaljobfile.out corsika_showertime
# - - - - opt. check `SIM??????.*`-files on antennas and gamma factor:
gammfact=`echo "0.0"`
gexpfact=`echo "0.0"`
jantenns=`echo "0"`
simreasinp=`echo $jrunnr | awk '{ printf("SIM%06d.inp",$1) }'`
if [ -e "$simreasinp" ] ; then
  jantenns=`grep "AntennaPosition" "SIM$jrunnr.list" | wc | awk '{printf("%d",$1)}'`
  jgamfact=`grep "gamma" "SIM$jrunnr.list" | wc | awk '{printf("%d",$1)}'`
  if [ $jgamfact -gt 0 ] ; then
    head -1 "SIM$jrunnr.list" > jobgamfact.txt
    gammfact=`cat "jobgamfact.txt" | awk '{printf("%s",$8)}'`
    gexpfact=`cat "jobgamfact.txt" | awk '{printf("%s",$9)}'`
    rm -f jobgamfact.txt
  fi
  echo "         ... writing number of antennas to $jobferr: " $jantenns
fi
# - - - - copy infos to Job??????_ik3.[err,out] output files:
timetext=`echo $jrunnr | awk '{ printf("time.txt%06d",$1) }'`
echo $nfiles | awk '{printf("%7d%7d\n",$1,$1)}' > $jobferr
echo $jantenns $gammfact $gexpfact | awk '{printf("%7d%7.1f%7s\n",$1,$2,$3)}' >> $jobferr
gfortran -O0 -fbounds-check summe.f -o summe
du -k $datname-* | ./summe > summe.out
head -1 job-file > $jobfout
tail -1 job-file >> $jobfout
echo $nfiles | awk '{printf("%10d Files \n",$1)}' >> $jobfout
echo $secnds | awk '{printf("%10d sec \n",$1)}' >> $jobfout 
cat $jobik >> $jobfout
cat $parallel >> $jobfout
cat summe.out >> $jobfout
# - - - - now create info file `time.txt00ijkl`:
gfortran -O0 -fbounds-check totaltimeik3.f -o totaltimeik3
./totaltimeik3 < $jobfout
rm totaljobfile.out
cat $timetext
fi # end-of if [ ! -e "$timelist" ];
exit
