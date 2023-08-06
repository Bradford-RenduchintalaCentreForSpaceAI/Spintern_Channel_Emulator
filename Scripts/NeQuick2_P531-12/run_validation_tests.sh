#!/bin/bash

# =============================================================================
# Bash script to run NeQuick2 P.531 validation input/output data files. 
# The script calls the NeQVal driver for each of the input data files from 
# the reference input folder. Description of validation procedure and file formats 
# is found in TESTING_VALIDATION.txt
# =============================================================================

BASE_DIR="TestingValidation/validation"

# list of files
FILES=$(ls $BASE_DIR/in/*.dat)

echo "Starting validation..."
for FILE in $FILES;
do
	FILENAME=${FILE##*/}
	PREFIX=${FILENAME%%in*} 

	OUTFILENAME=$PREFIX"out_063.dat"
	./NeQVal $FILE 63 > $BASE_DIR"/out/"$OUTFILENAME
	echo $FILENAME " --> " $OUTFILENAME

	OUTFILENAME=$PREFIX"out_128.dat"
	./NeQVal $FILE 128 > $BASE_DIR"/out/"$OUTFILENAME
	echo $FILENAME " --> " $OUTFILENAME

	OUTFILENAME=$PREFIX"out_193.dat"
	./NeQVal $FILE 193 > $BASE_DIR"/out/"$OUTFILENAME
	echo $FILENAME " --> " $OUTFILENAME
done

echo
echo "Results saved in " $BASE_DIR"/out/"

