#!/bin/bash

cd build
make -j8
cd ..
chrispi0 configfiles/chrisPhysics-Pi0.dat
