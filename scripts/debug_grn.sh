

cp src/VertexSystem.cpp src/VertexSystem.h src/basicGRN.cpp src/basicGRN.h src/run_grn.cpp src/GXMatrix.h src/plotOps_loop.py $1
cp param_files/$2.grn $1
cp param_files/$3.vp $1
cd $1


g++ -ggdb -std=c++11 -c VertexSystem.cpp
g++ -ggdb -std=c++11 -c basicGRN.cpp
g++ -ggdb -std=c++11 -c run_grn.cpp
g++ -o vertex run_grn.o basicGRN.o VertexSystem.o
gdb --args vertex $1 $2 $3 #winglike2_s3.25_20x20_n0.4 100000000 1000000



g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c basicGRN.cpp
g++ -std=c++11 -c run_grn.cpp  #test_grn.cpp
#g++ -o matrix test_matrix.o GXMatrix.o
g++ -o grn basicGRN.o run_grn.o VertexSystem.o #test_grn.o

./grn $1 $2 $3

