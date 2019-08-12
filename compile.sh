
g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c vertex_cmm.cpp
g++ -o vertex vertex_cmm.o VertexSystem.o
./vertex $1 $2 $3 #>$1'.log' 
#sleep 20
python plotOps_loop.py $1'_moved_' 0 "$(($2 / $3))" $4
rm *.o
rm vertex

mkdir $1
mv $1'_moved'* $1

cd $1
mkdir t1_transitions
mv ../T1_REPORT* t1_transitions
cd ../
