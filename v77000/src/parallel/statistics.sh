#/bin/bash -l

set -x
let i=0
declare -a R
while read line
do
a=`echo $line | awk '{print $7}'`
b=`echo $line | awk '{print $1}'`
R[$i]=`echo $a`
S[$j]=`echo $b`
#echo ${R[$i]} 
#echo ${S[$j]}
let i=$i+1
let j=$j+1
done < $file

i=0
let k=0
declare -a x;x=()
typeset -i T
M=`head -1 $file | awk '{print $1}'`
N=`tail -1 $file | awk '{print $8}'`

for U in ${S[@]}
do
  for (( T=`echo $M` ; "$T" <= "$N" ; (( T++ )) ))
  do
      if [ "$U" -le "$T" ] && [ "`expr $U + ${R[$i]}`" -gt "$T" ]
     then
       x[$k]=$((x[$k]+1))
     fi
  let k=$k+1
  done
let i=$i+1
k=0 
done

for member in ${x[*]}
do
 echo $member | cat >> stat-$file
done


