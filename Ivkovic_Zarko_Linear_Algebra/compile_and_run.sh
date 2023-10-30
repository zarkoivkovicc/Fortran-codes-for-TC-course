#!/bin/bash
mkdir build
mkdir bin
mkdir bin/Part1
mkdir bin/Part2
mkdir OUTPUT
mkdir OUTPUT/Part1
mkdir OUTPUT/Part2
cd build
cmake ..
make
cd ..
rm -rf build
# Part1
cd ./bin/Part1/
for i in 1 2 3 4; do
mkdir ../../OUTPUT/Part1/${i}
cp ../../INPUT/Part1/points_${i}.data .
echo $i 3 | ./part1 >> ../../OUTPUT/Part1/r_2
mv *.data ../../OUTPUT/Part1/${i}/.
done
mkdir ../../OUTPUT/Part1/5
cp ../../INPUT/Part1/points_5.data .
echo 5 3 | ./part1 >> ../../OUTPUT/Part1/r_2
mv *.data ../../OUTPUT/Part1/5/.
cp ../../SCRIPTS/part1.py ../../OUTPUT/Part1/.
#Part2
cd ../Part2/
for i in 'butadiene' 'benzene' 'naphthalene' 'pyrene'; do
mkdir ../../OUTPUT/Part2/${i}
cp ../../INPUT/Part2/${i}.xyz .
echo $i 0 | ./part2
mv ${i}.* ../../OUTPUT/Part2/${i}/.
done
mkdir ../../OUTPUT/Part2/cyclopentyl
cp ../../INPUT/Part2/cyclopentyl.xyz .
echo cyclopentyl -1 | ./part2
mv cyclopentyl.* ../../OUTPUT/Part2/cyclopentyl/.
cd ../..
# Plotting Part1
echo "Do you want to make plots for the Part1 ? [Y/n]"
read x
if [ $x = 'Y' ] 
then
	cp SCRIPTS/part1.py OUTPUT/Part1
	cd OUTPUT/Part1
	python part1.py
	cd ../..
fi
echo "Do you want to make plots for the Part2 ? [Y/n]"
read x
if [ $x = 'Y' ] 
then
	cp SCRIPTS/*.m OUTPUT/Part2
	cd OUTPUT/Part2
	for i in 'butadiene' 'benzene' 'naphthalene' 'pyrene' 'cyclopentyl'; do
	mv ${i}_gen_script.m ${i}/.
	mkdir ${i}/plot
	mkdir ${i}/plot/data
	mkdir ${i}/plot/images
	cd ${i}
	octave ${i}_gen_script.m
	gnuplot plot/data/*.gnu
	cd ..
	cp ../../SCRIPTS/part2.py .
	python part2.py
	done
	
fi
