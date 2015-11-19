#!/bin/bash

rm -rf ./deploy 2>/dev/null

for program in graph_checker redumis ; do 
scons program=$program variant=optimized -j 8 
if [ "$?" -ne "0" ]; then 
        echo "compile error in $program. exiting."
        exit
fi
done

mkdir deploy
cp ./optimized/redumis deploy/
cp ./optimized/graphchecker deploy/

rm -rf ./optimized

