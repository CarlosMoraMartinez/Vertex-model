
#ARGS:
#1: initial conditions name
#2: name of param file (without .vp)
#3: Total number of accepted moves
#4: Write output every N accepted moves

###

mkdir $1
cp src/VertexSystemMinimal.cpp src/VertexSystem.h src/vertex_cmm.cpp src/plotOps_loop.py 'param_files/'$2'.vp' $1

cd $1
g++ -ggdb -std=c++11 -c VertexSystemMinimal.cpp
g++ -ggdb -std=c++11 -c vertex_cmm.cpp
g++ -o vertex vertex_cmm.o VertexSystemMinimal.o
gdb --args vertex $1 $2 $3 $4 #>$1'.log' 
