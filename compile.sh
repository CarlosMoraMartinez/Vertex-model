

# mkdir $1
cp VertexSystem.cpp VertexSystem.h vertex_cmm.cpp plotOps_loop.py $1
cd $1
g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c vertex_cmm.cpp
g++ -o vertex vertex_cmm.o VertexSystem.o
./vertex $1 $2 $3 #>$1'.log' 
#sleep 20
python plotOps_loop.py $1'_moved_' 0 "$(($2 / $3))" $4
ffmpeg -framerate 10 -i $1'_moved_'%d.png -codec copy $1'_moved_'.avi

rm *.o
rm vertex

#mv $1'_moved'* $1

#cd $1
mkdir t1_transitions
mv T1_REPORT* t1_transitions
cd ../
