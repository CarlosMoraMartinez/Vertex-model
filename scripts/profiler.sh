



#ARGS:
#1: initial conditions name
#2: name of param file (without .vp)
#3: Total number of accepted moves
#4: Write output every N accepted moves

###




#g++ -std=c++11 -c profile_test.cpp -pg
#g++ -o profiletest profile_test.o -pg
#./profiletest 
#gprof profiletest gmon.out >profile.out
#gprof -l profiletest gmon.out >profile_byline.out

#cat profile.out


mkdir $1
cp src/VertexSystem.cpp src/VertexSystem.h src/vertex_cmm.cpp src/plotOps_loop.py 'param_files/'$2'.vp' $1

cd $1
g++ -std=c++11 -pg -O2 VertexSystem.cpp vertex_cmm.cpp -o vertex
#g++ -std=c++11 -c - pg vertex_cmm.cpp -pg
#g++ -o vertex -pg vertex_cmm.o VertexSystem.o
./vertex $1 $2 $3 $4 #>$1'.log' 
gprof vertex gmon.out >profile.out
gprof -l vertex gmon.out >profile_byline.out

