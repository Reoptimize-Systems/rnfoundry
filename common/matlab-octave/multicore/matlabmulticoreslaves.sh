#!/bin/bash
echo "The multicore directory:"
echo $1 # the multicore directory
echo "The quit date for the slave:"
echo $2 # the quit date for the slave
echo "The start directory for the process"
echo $3 # the start directory for the process
echo "The delay before starting the process"
echo $4 # The delay before starting the process
echo "The slave ID number"
echo $5 # The slave ID number
echo "The nice (process priority) value"
echo $6 # The slave ID number

echo waiting $4 seconds
sleep $4

echo "Running Matlab in batch mode"
# Start the matlab process in batch mode
# and limited to a single thread, we also reduce the process priority a bit 
# using nice to try to make 
nice -n $6 matlab -singleCompThread -nosplash -sd "$3" -batch "disp(pwd);disp('Successfully Completed Startup');mcore.startslave('MulticoreSharedDir', '$1', 'QuitMatlabOnExit', true, 'QuitDatenum', datenum($2), 'SlaveID', $5, 'DebugMode', false);quit;"
