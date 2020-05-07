g++ -ggdb -std=c++11 -c VertexSystem.cpp
g++ -ggdb -std=c++11 -c vertex_cmm.cpp
g++ -o vertex vertex_cmm.o VertexSystem.o
gdb --args vertex gm0_s3.0_3x3_n0.3 test_params 200000 10000 #winglike2_s3.25_20x20_n0.4 100000000 1000000

