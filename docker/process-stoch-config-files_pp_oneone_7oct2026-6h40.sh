#!/bin/bash

# ALLFOLDERS=( /local2/lareri/library/pointToPoint/oneSrcOneSnk/flat/smlInstTigh )

 ALLFOLDERS=( /local2/lareri/library/pointToPoint/oneSrcOneSnk/flat/smlInstLoo )

# ALLFOLDERS=( /local2/lareri/library/pointToPoint/oneSrcOneSnk/flat/larInstTigh )

# ALLFOLDERS=( /local2/lareri/library/pointToPoint/oneSrcOneSnk/flat/larInstLoo )

for dirpath in "${ALLFOLDERS[@]}"
do

	counter=0
	diagfilepath="${dirpath}/diags-gener.txt"
	timingsfilepath="${dirpath}/timings_gener.txt"
	echo "diagfilepath: ${diagfilepath}"
	echo "timingsfilepath: ${timingsfilepath}"
	
	for configfile in "$dirpath"/*.par;
	do
		if [ -s "$configfile" ]; then
		
			echo "configfile: ${configfile}"			
			SECONDS=0
			
			./stoch_gen +v +F ${configfile} &>> ${diagfilepath}
			
			DURATION=${SECONDS}
			
			printf "%s;%d\n" ${configfile} ${DURATION} &>> ${timingsfilepath}
			
			((counter++))
		else			
			echo "ERROR: configuration file ${configfile} could not be opened" >> ${diagfilepath}
			
		fi
	done
	
	echo "counter: ${counter}" >> ${diagfilepath}
	
done