#!/bin/sh

rm $0

chromium-browser /home/rcrozier/src/rnfoundry-hg/wave/doc/sphinx/_build/html/getting_started.html

echo "

------------------
(program exited with code: $?)" 		


