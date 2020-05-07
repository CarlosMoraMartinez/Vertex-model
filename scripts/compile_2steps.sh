

#ARGS:
#1: initial conditions name
#2: name of param file (without .vp)
#3: Total number of accepted moves
#4: Write output every N accepted moves

###

mkdir $1
cp src/VertexSystem.cpp src/VertexSystem.h src/vertex_2steps.cpp src/plotOps_loop.py 'param_files/'$2'.vp' 'param_files/'$5'.vp' $1

cd $1
g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c vertex_2steps.cpp
g++ -o vertex vertex_2steps.o VertexSystem.o
./vertex $1 $2 $3 $4 $5 #>$1'.log' 
#sleep 20
python plotOps_loop.py -i $1'_moved_' -s "$((2*$3 / $4 ))" -e -1 -l 1 #-v 1 #-s 0 -e "$(($3 / $4 + 1))" 
#newname=$1'_movedspr_'"$(($3 / $4 ))"'_moved_'
#python plotOps_loop.py -i $newname -s "$(($3 / $4 ))" -e -1 -l 1 #-v 1 #-s 0 -e "$(($3 / $4 + 1))" 
ffmpeg -framerate 10 -i $1'_moved_'%d.png -codec copy $1'_moved_'.avi
#ffmpeg -framerate 10 -i $newname%d.png -codec copy $newname.avi


rm *.o
rm vertex

#mv $1'_moved'* $1

#cd $1
#mkdir t1_transitions
#mv T1_REPORT* t1_transitions
#cd ../
