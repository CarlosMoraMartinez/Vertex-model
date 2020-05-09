

#ARGS:
#1: initial conditions name
#2: name of param file (without .vp)
#3: Total number of accepted moves
#4: Write output every N accepted moves

###

mkdir $1
cp src/VertexSystem.cpp src/VertexSystem.h src/vertex_cmm.cpp src/plotOps_loop.py 'param_files/'$2'.vp' $1

cd $1
g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c vertex_cmm.cpp
g++ -o vertex vertex_cmm.o VertexSystem.o
./vertex $1 $2 $3 $4 #>$1'.log' 
#sleep 20
#python3 plotOps_loop.py -i $1'_moved_' -s 0 -e "$(($3 / $4 + 1))"
python plotOps_loop.py -i $1'_moved_' -s "$(($3 / $4 ))" -e -1 -l 1
ffmpeg -framerate 10 -i $1'_moved_'%d.png -codec copy $1'_moved_'.avi

rm *.o
rm vertex

#mv $1'_moved'* $1

#cd $1
#mkdir t1_transitions
#mv T1_REPORT* t1_transitions
#cd ../


