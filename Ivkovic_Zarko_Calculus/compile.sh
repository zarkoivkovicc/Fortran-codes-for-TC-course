#!/bin/bash
mkdir build
mkdir bin
mkdir bin/Part1
mkdir bin/Part2
cd build
cmake ..
make
cd ..
rm -rf build
