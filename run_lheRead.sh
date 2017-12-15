#!/bin/sh

types=("sp0_300" "sp0_400" "sp0_500" "sp0_600" "sp0_700" "sp0_800" "sp0_900" "sp0_1000" "sp0_1100" "sp0_1200" "sp0_1400" "sp1_500" "sp1_600" "sp1_700" "sp1_800" "sp1_900" "sp1_1000" "sp1_1100" "sp1_1200" "sp1_1400" "sp2_500" "sp2_600" "sp2_700" "sp2_800" "sp2_1000" "sp2_1100" "sp2_1200" "sp2_1400") 

for type in ${types[@]}
do
	python lheRead.py $type
	#python lheRead_bkg.py $type
done
