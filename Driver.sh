#!/bin/bash

#Simulations list
Sims=("1-Series1" "2-Series2" "3-Series3" "4-Series4" "5-Series5")
Scripts=("AngleLengthDist.py" )

for sim in "${Sims[@]}"
do
	for script in "${Scripts[@]}"
	do
		cp $script $sim/dump/z-cId
	done
done

for sim in "${Sims[@]}"
do
	cd $sim/dump/z-cId
	for script in "${Scripts[@]}"
	do
		python2 $script & # Remove the "&" sign to turn off the parallel processing
	done
	cd ../../../
done
wait # comment this to turn off the parallel processing


python2 AngleLengthDistAllSims3.py



