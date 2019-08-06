#/bin/bash -l

set -x

let x=1

while read line
do
a=`echo $line | awk '{print $2}'`
a=`echo ${a/.*/}`
b=`date --utc --date "1970-1-1 0:$a" +%s`
c=`echo $line | awk '{print $1}'`
b=`expr $b + $x`
sum=`expr $c + $b`
echo $line $b $sum | cat >> new-job-file
done < job-file

export file=new-job-file
./statistics.sh

