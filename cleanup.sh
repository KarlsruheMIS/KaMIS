#!/bin/bash

# scons program=library variant=optimized -j 4 -c
scons program=redumis variant=optimized -j 4 -c 

# rm extern/KaHIP/libkahip.a
rm -rf deploy
rm -rf optimized
