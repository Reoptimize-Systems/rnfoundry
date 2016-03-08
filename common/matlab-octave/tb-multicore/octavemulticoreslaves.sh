#!/bin/bash
# 
# bash script for launching octave slave processes
#
#

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

OCT_COMMAND="fprintf(1, 'System PATH variable is:\\n%s\\n', getenv('PATH')); \
setenv('HOSTNAME','$(hostname)'); \
fprintf(1, 'HOSTNAME environment variable is: %s\\n', getenv('HOSTNAME')); \
cd('$3'); \
disp(pwd); \
ls; \
save_default_options ('-binary'); \
disp('Successfully Completed Startup'); \
startmulticoreslave2('$1', 1, datenum($2)); \
quit;"

# Start the octave process 
octave-cli --eval "${OCT_COMMAND}"
