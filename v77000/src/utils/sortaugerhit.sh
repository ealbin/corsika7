#!/bin/bash
#
# sortaugerhit.sh:
# ===============
# create list of file names and run 'sortaugerhit' program:
# names of parallel sub paths must always be 'cskiiiiii/';
# ------------------------------------------------------------------------
# gfortran -O0 -fbounds-check sortaugerhit.f -o sortaugerhit
# ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
# ------------------------------------------------------------------------
# usage: ./sortaugerhit.sh
#        ./sortaugerhit.sh 2084
#        ./sortaugerhit.sh 672012
#        ./sortaugerhit.sh DAT672012
#        ./sortaugerhit.sh DAT672012.part
# ------------------------------------------------------------------------
# cd csk002084/
# cp ../sortaugerhit* . 
# ls -1 DAT002084-0????? > sortaugerhit.i002084
# ls -1 DAT002084-* | grep t -v | grep n -v > sortaugerhit.i002084
# ls -1 DAT002084-[0-8]????????-????????? > sortaugerhit.i002084
# ./sortaugerhit < sortaugerhit.i002084 > sortaugerhit.out002084
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
#
# - - parallel path or simply one particle data file:
if [ ! -e "sortaugerhit.f" ] ; then
  echo "          ----- missing source file 'sortaugerhit.f'!"
  exit 1
fi
# - - compile and link fortran code considering the login node:
if [ -e "sortaugerhit.f" ] ; then
  if [ `echo $HOSTNAME | cut -b 1-4` == "iklx" ] ; then
    gfortran -O0 -fbounds-check sortaugerhit.f -o sortaugerhit
  elif [ `echo $HOSTNAME | cut -b 1-4` == "uc1n" ] ; then
    ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
  elif [ `echo $HOSTNAME | cut -b 1-4` == "fh1n" ] ; then
    ifort -C -O0 -check bounds sortaugerhit.f -o sortaugerhit
  fi
fi
# - - test on available tabular file 'sortaugerhit.pmsoutab':
if [ ! -e "sortaugerhit.pmsoutab" ] ; then
  if [ -e "../sortaugerhit.pmsoutab" ] ; then
    sleep 1
    # cp ../sortaugerhit.pmsoutab .
  else
    sleep 1
    # echo " tabular file '../sortaugerhit.pmsoutab' not available. stop."
    # exit 1
  fi
fi
#
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#
# - - test on a given argument applying the script:
jobcode=0
if [ $# -eq 1 ] ; then
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - argument is 'DATiiiiii[.part]' or 'iiiiii' or 'cskiiiiii/':
  scriptarg=`echo $1`
  echo $scriptarg > jobfilearg.tmp
  jbeg=0
  if [ ${#1} -le 6 ] ; then
    # - - - - - - argument should be a run number (up to 6 digits):
    jrunnr=`echo $scriptarg | awk '{printf("%06d",$1)}'`
    if [ -e "DAT$jrunnr" ] ; then
      datfile=`echo $jrunnr | awk '{printf("DAT%06d",$1)}'`
      echo "DAT$jrunnr" > sortaugerhit.i$jrunnr
      jobcode=1
    else
      if [ -e "DAT$jrunnr.part" ] ; then
        datfile=`echo $jrunnr | awk '{printf("DAT%06d.part",$1)}'`
        echo "DAT$jrunnr.part" > sortaugerhit.i$jrunnr 
        jobcode=2
      else
        # - - - - - - test on corsika path 'cskiiiiii/':
        pwd > jobfileerr.tmp  
        datpath=`cat jobfileerr.tmp`
        jcsk=`awk -v a="$datpath" -v b="csk" 'BEGIN{print index(a,b)}'`
        if [ $jcsk -gt 3 ] ; then
          sed -i "/csk/ s/csk/   /" jobfileerr.tmp
          jrunpa=`cat jobfileerr.tmp | awk '{printf("%06d",$2)}'`
          if [ $jrunpa -eq $jrunnr ] ; then 
            jobcode=4
            if [ -e DAT$jrunnr-000001 ] ; then
              ls -1 DAT$jrunnr-0????? > sortaugerhit.i$jrunnr
              ls -l DAT$jrunnr-0????? | awk '{total+=$5}; END { print total > jobfilsize.tmp }'
            else
              ls -1 DAT$jrunnr-[0-8]????????-????????? > sortaugerhit.i$jrunnr
              ls -l DAT$jrunnr-[0-8]????????-????????? | awk '{total+=$5}; END { print total > jobfilsize.tmp }'
            fi
            cat jobfilsize.tmp
            jsize=`cat jobfilsize.tmp | awk '{ printf("%d",$1) }'`
            rm -rf jobfilsize.tmp
            echo "./sortaugerhit < sortaugerhit.i$jrunnr > sortaugerhit.out$jrunnr &"
            echo "          (check elapsed time on created files 'DAT$jrunnr-augerhit*')"
          else # jrunpa != jrunnr
            scriptarg=$jrunpa
            echo "          runnr $jrunnr not available here in 'csk$jrunpa/'"
            jrunnr=$jrunpa
            jobcode=0
          fi
        fi
      fi
    fi
  else # length of argument > 6 characters:
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - argument may be a file name or the current path name:
    datfile=`basename $scriptarg | awk '{printf("%s",$1)}'`
    workfile=`echo $datfile`
    jbeg=`awk -v a="$workfile" -v b="_" 'BEGIN{print index(a,b)}'`
    if [ $jbeg -gt 0 ] ; then
      # set pseudo run number and save name of particle data file:
      jrunnr=987654
      echo $scriptarg > sortaugerhit.i$jrunnr
      # take care on names `prE16M562_41`, other primaries co, gm, he, fe, si.
      jobcode=1
    elif [ `cat "jobfilearg.tmp" | cut -b 1-3` != "csk" ] ; then
      # - - - - - - - no 'cskiiiiii/' path name given:
      echo "$scriptarg" > jobfileerr.tmp
      sed -i "/part/ s/\.part/ part/" jobfileerr.tmp
      sed -i "/DAT/ s/DAT/   /" jobfileerr.tmp
      jrunnr=`cat jobfileerr.tmp | awk '{printf("%06d",$1)}'`
      cat jobfilearg.tmp
      if [ ! -e "$scriptarg" ] ; then
        jrunnr=999999
        scriptarg=0
        echo "          unknown run number or data file not in this path."
        jobcode=0 # 9
      else
        if [ -e "DAT$jrunnr" ] ; then 
          datfile=`echo $jrunnr | awk '{printf("DAT%06d",$1)}'`
          echo "DAT$jrunnr" > sortaugerhit.i$jrunnr
          jobcode=1
        else
          if [ -e "DAT$jrunnr.part" ] ; then
            datfile=`echo $jrunnr | awk '{printf("DAT%06d.part",$1)}'`
            echo "DAT$jrunnr.part" > sortaugerhit.i$jrunnr 
            jobcode=2
          else
            if [ -e "DAT$jrunnr-999999" ] ; then
              datfile=`echo $scriptarg | awk '{printf("%s",$1)}'`
              echo "$datfile" > sortaugerhit.i$jrunnr
              jobcode=3
            else 
              # check conditions??
              jrunnr=999999
              jobcode=0
            fi
          fi
        fi
      fi
    else 
      # - - - - - - - found name 'cskiiiiii/' as argument:
      echo "$scriptarg" > jobfiletxt.tmp
      sed -i "/csk/ s/csk/   /" jobfiletxt.tmp
      jrunna=`cat jobfiletxt.tmp | awk '{printf("%06d",$1)}'`
      pwd > jobfileerr.tmp
      datpath=`cat jobfileerr.tmp`
      jcsk=`awk -v a="$datpath" -v b="csk" 'BEGIN{print index(a,b)}'`
      if [ $jcsk -gt 3 ] ; then
        sed -i "/csk/ s/csk/   /" jobfileerr.tmp
        jrunnr=`cat jobfileerr.tmp | awk '{printf("%06d",$2)}'`
      else
        jrunnr=999999
      fi
      if [ $jrunna -eq $jrunnr ] ; then
        echo "$scriptarg" > sortaugerhit.i$jrunnr
      else
        jrunnr=999999
        scriptarg=0
        echo "          mismatching run number $jrunnr and path 'csk$jrunpa/'"
        jobcode=0 # 9
      fi
    fi
  fi
  # - - - - check run number and jobcode:
  if [ $jrunnr -lt 999999 ] ; then
    if [ $jobcode -gt 0 ] ; then
      if [ $jobcode -eq 1 ] ; then
        if [ $jbeg -eq 0 ] ; then
          ls -l DAT$jrunnr | awk '{ printf("%s\n",$5) }' > jobfilsize.tmp
        else
          ls -l $scriptarg | awk '{ printf("%s\n",$5) }' > jobfilsize.tmp  
        fi
      elif [ $jobcode -eq 2 ] ; then
        ls -l DAT$jrunnr.part | awk '{ printf("%s\n",$5) }' > jobfilsize.tmp
      elif [ $jobcode -eq 3 ] ; then
        ls -l DAT$jrunnr-999999 | awk '{ printf("%s\n",$5) }' > jobfilsize.tmp
      fi
      jsize=`cat jobfilsize.tmp | awk '{ printf("%d",$1) }'`
      rm -rf jobfilsize.tmp
      if [ $jbeg -eq 0 ] ; then
        echo "./sortaugerhit < sortaugerhit.i$jrunnr > sortaugerhit.out$jrunnr &"
      else
        echo "./sortaugerhit < sortaugerhit.i$jrunnr > sortaugerhit.out$scriptarg &"
      fi
      echo "          (check elapsed time on created files 'DAT$jrunnr-augerhit*')"
    fi
  fi
fi
#
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#
# - - no argument means 'cskiiiiii/' path or argument is a number or a file:
if [ $# -eq 0 ] || [ $jobcode -eq 0 ] ; then
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pwd > jobfilearg.tmp
  pwd > jobfileerr.tmp
  datpath=`cat jobfileerr.tmp`
  jcsk=`awk -v a="$datpath" -v b="csk" 'BEGIN{print index(a,b)}'`
  if [ `cat "jobfilearg.tmp" | cut -b 1-3` == "csk" ] ; then
    sed -i "/csk/ s/csk/   /" jobfileerr.tmp
    jrunnr=`cat jobfileerr.tmp | awk '{printf("%06d",$1)}'` 
  else
    if [ $jcsk -gt 3 ] ; then
      sed -i "/csk/ s/csk/   /" jobfileerr.tmp
      jrunnr=`cat jobfileerr.tmp | awk '{printf("%06d",$2)}'`
    else
      jrunnr=999999
      echo "          unknown run number, and being not within 'cskiiiiii/' path."
    fi
  fi
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if [ $jrunnr -lt 999999 ] ; then
    if [ -e DAT$jrunnr-999999 ] ; then
      ls -1 DAT$jrunnr-999999 > sortaugerhit.i$jrunnr
      ls -l DAT$jrunnr-0????? | awk '{total+=$5}; END { print total }' > jobfilsize.tmp
    else
      if [ -e DAT$jrunnr-000001 ] ; then
        ls -1 DAT$jrunnr-0????? > sortaugerhit.i$jrunnr
        ls -l DAT$jrunnr-0????? | awk '{total+=$5}; END { print total }' > jobfilsize.tmp
      else
        ls -1 DAT$jrunnr-[0-8]????????-????????? > sortaugerhit.i$jrunnr
        ls -l DAT$jrunnr-[0-8]????????-????????? | awk '{total+=$5}; END { print total }' > jobfilsize.tmp
      fi
    fi
    jsize=`cat jobfilsize.tmp | awk '{ printf("%d",$1) }'`
    rm -rf jobfilsize.tmp
    echo "./sortaugerhit < sortaugerhit.i$jrunnr > sortaugerhit.out$jrunnr &"
    echo "          (check elapsed time on created files 'DAT$jrunnr-augerhit*')"
  fi
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fi
#
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#
rm -f jobfilearg.tmp
rm -f jobfileerr.tmp
