#!/bin/bash
echo "The multicore directory:"
echo $1 # the multicore directory
echo "The quit date for the slave:"
echo $2 # the quit date for the slave
echo "The start directory for the process"
echo $3 # the start directory for the process
echo "The delay before starting the process"
echo $4 # the start directory for the process

echo waiting $4 seconds
sleep $4

echo "Running Maltab without jvm"
# Start the matlab process with no gui, no java (until this is fixed),
# and limited to a single thread
matlab -batch "cd('$3');disp(pwd);ls;startup;disp('Successfully Completed Startup');mcore.startslave('$1', 1, datenum($2));quit;"

