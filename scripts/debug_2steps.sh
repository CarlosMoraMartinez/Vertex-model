

#ARGS:
#1: initial conditions name
#2: name of param file (without .vp)
#3: Total number of accepted moves
#4: Write output every N accepted moves

###

mkdir $1
cp src/VertexSystem.cpp src/VertexSystem.h src/vertex_2steps.cpp src/plotOps_loop.py 'param_files/'$2'.vp' 'param_files/'$5'.vp' $1

cd $1
g++ -ggdb -std=c++11 -c VertexSystem.cpp
g++ -ggdb -std=c++11 -c vertex_2steps.cpp
g++ -o vertex vertex_2steps.o VertexSystem.o
gdb --args ./vertex $1 $2 $3 $4 $5 #>$1'.log' 

rm *.o
rm vertex

