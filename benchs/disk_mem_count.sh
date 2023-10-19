#!/bin/bash 
EDIR=$( pwd )
start=$( du -ks ${EDIR} 2> /dev/null| cut -f 1 )
#echo "start",$start
max=$start
_user=ppeterlo
max_memory=$( ps hux -U ${_user} | awk -v user=${_user} -v total=$TOTAL '{ sum += $6 } END {print sum;}')
"$@" &
pid=$!

while :
do
       ps p $pid > /dev/null
       status=$?
       sleep 0.5
	if [ $status -eq 1 ]
        then
            break
        fi
        disk=$( du -ks ${EDIR}  2> /dev/null| cut -f 1 )
        memory=$( ps hux -U ${_user} | awk -v user=${_user} -v total=$TOTAL '{ sum += $6 } END {print sum;}')
        if [ $disk -gt $max ] 
        then    
               max=$disk
        fi
        if [ $memory -gt $max_memory ]
        then
               max_memory=$memory
        fi
       
done

disk=$( du -ks ${EDIR}  2> /dev/null| cut -f 1 )
if [ $disk -gt $max ] 
then    
       max=$disk
fi 
du -ks ${EDIR} 2> /dev/null 
# echo "disk",$disk
# echo "start=",$start, "max",$max
echo "************ disk used ( KB )" $(( $max - $start )) "*******************"
echo "************ max memory ( KB )" ${max_memory} "********************"
disk_used=$(( $max - $start ))
disk_GB=$((${disk_used}/1048576))
mem_GB=$((${max_memory}/1048576))
echo "************ disk used ( GB )" ${disk_GB} "*******************"
echo "************ max memory ( GB )" ${mem_GB} "********************"


