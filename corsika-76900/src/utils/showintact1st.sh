#!/bin/bash
#
# showintact1st.sh:
# =================
#     read one/first particle data file of parallel corsika simulations 
#     and print out tabular of current heights of first interaction 
#     of each air shower simulation, i.e. an additional tabular info to
#     text file `showparallel.bwc-work-jobinfos` created by script
#     `showparallel.sh`;
# ------------------------------------------------------------------------
# check:
#     gfortran -O0 -fbounds-check showintact1st.f -o showintact1st
#     ifort -C -O0 -check bounds showintact1st.f -o showintact1st
# ------------------------------------------------------------------------
# usage: ./showintact1st.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
datnampart=`echo $1`
if [ $# -eq 1 ] ; then
  echo $datnampart > /dev/null
fi
# ------------------------------------------------------------------------
#
computype="ik3"
if [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
  computype="bwc"
else
  if [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
    computype="fh1"
  fi
fi
#
userhomepath=`echo $USER | awk '{printf("/cr/users/%s",$1)}'`
if [ $computype == "bwc" ] ; then
  userhomepath=`echo $USER | awk '{printf("/work/kit/ikp/%s",$1)}'`
else
  if [ $computype == "fh1" ] ; then
    userhomepath=`echo $USER | awk '{printf("/home/kit/ikp/%s",$1)}'`
  fi
fi
#
# ------------------------------------------------------------------------
#
pwdnamkit=`echo $PWD | cut -b 6-13 | awk '{printf("%s",$1)}'`
#
# ------------------------------------------------------------------------
if [ $pwdnamkit == "/kit/ikp" ]; then
  ifort -C -O0 -check bounds showintact1st.f -o showintact1st
  ls -1d csk?????? > showintact1st.directs
  rm -f showintact1st.partfile
  touch showintact1st.partfile
  # - - - loop on all sub paths of parallel simulations:
  for name in $( cat "showintact1st.directs" ) ; do
    echo $name | awk '{printf("%s\n",$1)}' > showintact1st.tmp
    jrunnr=`cut -b 4-9 "showintact1st.tmp" | awk '{printf("%06d",$1)}'` 
    ls -1 $name/Job$jrunnr* 2> showintact1st.err | wc > showintact1st.tmp
    # to test sub paths use only: ls -1 $name/Job$jrunnr* | wc 
    njobfiles=`cat "showintact1st.tmp" | awk '{printf("%d",$1)}'`
    if [ $njobfiles -gt 0 ]; then
      # - - - - check h1stme only for finished parallel simulations:
      datpattern=`echo $name $jrunnr | awk '{printf("%s/DAT%06d-*",$1,$2)}'`
      # - - - - ignore all file names containing l,a,t:
      ls -1 $datpattern | grep -v l | grep -v a | grep -v t | tail -1 > showintact1st.partfile
      ./showintact1st < showintact1st.partfile
    fi
  done
  rm -f showintact1st.directs
  rm -f showintact1st.err
  rm -f showintact1st.h1stme
  rm -f showintact1st.partfile
  rm -f showintact1st.tmp
fi
# ------------------------------------------------------------------------
#
pwdnamik3=`echo $PWD | cut -b 1-4  | awk '{printf("%s",$1)}'`
#
# ------------------------------------------------------------------------
if [ $pwdnamik3 == "/cr/" ]; then
  userhomepath=`echo $USER | awk '{printf("/cr/users/%s",$1)}'`
  gfortran -O0 -fbounds-check "showintact1st.f" -o "$userhomepath/showintact1st"
  touch "$userhomepath/showintact1st.prtcls" 
  if [ "$datnampart" == "part" ] ; then
    ls -1 DAT??????.part >> "$userhomepath/showintact1st.prtcls"
  else
    ls -1 DAT?????? >> "$userhomepath/showintact1st.prtcls"
  fi
  jfiles=`wc "$userhomepath/showintact1st.prtcls" | awk '{printf("%d",$1)}'`
  if [ $jfiles -gt 0 ] ; then  
    # - - - - regular single particle data files:
    rm -f "$userhomepath/showintact1st.partfile"
    touch "$userhomepath/showintact1st.partfile"
    # - - - loop on all data files of sequential simulations:
    for name in $( cat "$userhomepath/showintact1st.prtcls" ) ; do
      echo $name | awk '{printf("%s\n",$1)}' > "$userhomepath/showintact1st.partfile"
      "$userhomepath/showintact1st" < "$userhomepath/showintact1st.partfile"
    done
    rm -f $userhomepath/showintact1st.prtcls
    rm -f $userhomepath/showintact1st.h1stme
    rm -f $userhomepath/showintact1st.partfile
    rm -f $userhomepath/showintact1st.tmp
    rm -f $userhomepath/showintact1st
  else # [ $jfiles -eq 0 ]
    # - - - - subdirectories of type csk?????? of parallel simulations:
    ls -1d csk?????? > "$userhomepath/showintact1st.directs"
    ifiles=`wc "$userhomepath/showintact1st.directs" | awk '{printf("%d",$1)}'`
    if [ $ifiles -gt 0 ] ; then 
      rm -f "$userhomepath/showintact1st.partfile"
      touch "$userhomepath/showintact1st.partfile"
      # - - - loop on all sub paths of parallel simulations:
      for name in $( cat "$userhomepath/showintact1st.directs" ) ; do
        echo $name | awk '{printf("%s\n",$1)}' > "$userhomepath/showintact1st.tmp"
        jrunnr=`cut -b 4-9 "$userhomepath/showintact1st.tmp" | awk '{printf("%06d",$1)}'` 
        ls -1 $name/Job$jrunnr* 2> "$userhomepath/showintact1st.err" | wc > "$userhomepath/showintact1st.tmp"
        # to test sub paths use only: ls -1 $name/Job$jrunnr* | wc 
        njobfiles=`cat "$userhomepath/showintact1st.tmp" | awk '{printf("%d",$1)}'`
        if [ $njobfiles -gt 0 ]; then
          # - - - - check h1stme only for finished parallel simulations:
          datpattern=`echo $name $jrunnr | awk '{printf("%s/DAT%06d-*",$1,$2)}'`
          # - - - - ignore all file names containing l,a,t:
          ls -1 $datpattern | grep -v l | grep -v a | grep -v t | tail -1 > "$userhomepath/showintact1st.partfile"
          "$userhomepath/showintact1st" < "$userhomepath/showintact1st.partfile"
        fi
      done
      rm -f "$userhomepath/showintact1st.err"
      rm -f "$userhomepath/showintact1st.h1stme"
      rm -f "$userhomepath/showintact1st.partfile"
      # rm -f "userhomepath/showintact1st.lastfile"
      # rm -f "$userhomepath/showintact1st.directs"
      rm -f "$userhomepath/showintact1st.tmp" 
      rm -f "$userhomepath/showintact1st" 
    fi
    rm -f "$userhomepath/showintact1st.prtcls"
    rm -f "$userhomepath/showintact1st.directs"
  fi
fi
# end-of script showintact1st.sh
