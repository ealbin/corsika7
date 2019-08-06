###########  script for primary local run #########

#!/bin/bash -l
#@ environment = $FILE;$EXEC;$JOBINFO;$mpiid;$RNNUM;$CHECK;
   
	output=`echo $RNNUM $mpiid | awk '{printf("./DAT%06d-000000000-%09d.lst", $1,$2)}'`
	end=`echo $RNNUM $mpiid | awk '{printf("./DAT%06d-000000000-%09d.end", $1,$2)}'`

	t=$(date '+%s');	echo $mpiid " 0 1 " $RNNUM " 0 " $mpiid " 0 0 " $t > ./status_start


   a=$(/usr/bin/time -f "\t%E Real-time \t%S TotalCPUseconds" 2>&1 ./$EXEC < $FILE > $output);a=$(date '+%s')$a;echo $a $output > $JOBINFO

	t=$(date '+%s');	echo $mpiid " 0 1 " $RNNUM " 0 " $mpiid " 0 0 " $t > ./status_finish

	endcheck=`grep -c "END OF RUN" $output`
	while [ $endcheck -ne 1 ]
        do
         ./$EXEC < $FILE > $output
        done

	echo "Initial Run" | cat > $end
	echo "Initial Run" | cat > $CHECK
