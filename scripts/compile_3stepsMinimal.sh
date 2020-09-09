

#ARGS:
#1: initial conditions name
#2-4: param files (without .vp)

###

mkdir $1
cp src/VertexSystemMinimal.cpp src/VertexSystem.h src/vertex_3steps.cpp src/plotOps_loop.py 'param_files/'$2'.vp' 'param_files/'$3'.vp' 'param_files/'$4'.vp' $1

cd $1
g++ -O2 -std=c++11 -c VertexSystemMinimal.cpp
g++ -O2 -std=c++11 -c vertex_3steps.cpp
g++ -o vertex vertex_3steps.o VertexSystemMinimal.o

./vertex $1 $2 $3 $4 

#python plotOps_loop.py -i $1'_moved_' -s "$((2*$3 / $4 ))" -e -1 -l 1 
python plotOps_loop.py -i $1'_moved_' -s 0 -e 1000 

ffmpeg -framerate 10 -i $1'_moved_'%d.png -codec copy $1'_moved_'.avi

rm *.o
rm vertex

