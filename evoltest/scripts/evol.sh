
#1: name of simulation
#2: param file wing growth
#3: param file hinge contraction
#4: wing initial condition (eg, bud2)
#5: wing to compare with (.border)

# bash evol.sh evoltest1 test2steps_1a test2steps_1b bud2 contour1
mkdir $1
cp ../src/VertexSystem.cpp ../src/VertexSystem.h ../src/vertex_2steps.cpp 'param_files/'$2'.vp' 'param_files/'$3'.vp' 'param_files/'$5'.border' $1
cp ./$4/$4.* $1
cd $1

g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c vertex_2steps.cpp
g++ -o vertex vertex_2steps.o VertexSystem.o

python ../scripts/tournament_selection.py -o $1 -m $2'.vp,'$3'.vp' -s $4 -t $5 -r 0.5 -n 10 -k 7 -x 100 -c 30 -w 5000 -z 1000
