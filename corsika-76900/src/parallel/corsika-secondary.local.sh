
###########  script for secondary local runs  #########

#!/bin/bash -l
#@ environment = $FILE;$NR;$CHECK;$NUM;$EXEC;$JOBINFO;$PARENT;$seed;$id1;$id2,$format;$nofinished;$nostarted;$mpiid;$RNNUM;
#set -x
	t=$(date '+%s');	echo $mpiid $PARENT " 0 " $RNNUM $seed $mpiid $id1 $id2 $t >> ./status_start

	output=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.lst", $1,$2)}'`
	end=`echo DAT$format-$seed- $mpiid | awk '{printf("./%s%09d.end", $1,$2)}'`

	 a=$(/usr/bin/time -f "\t%E Real-time \t%S TotalCPUseconds" 2>&1 ./$EXEC < $FILE > $output);a=$(date '+%s')$a;echo $a $output >> $JOBINFO
         
	endcheck=`grep -c "END OF RUN" $output`
	while [ $endcheck -ne 1 ]
        do
         ./$EXEC < $FILE > $output
        done

	t=$(date '+%s');	echo $mpiid $PARENT " 0 " $RNNUM $seed $mpiid $id1 $id2 $t >> ./status_finish

	
        touch $end
        echo "Secondary Run ID" $mpiid | cat >> $CHECK
	 echo "Secondary Run ID" $mpiid | cat >> $end