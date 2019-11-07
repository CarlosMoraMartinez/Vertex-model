#g++ -std=c++11 -c GXMatrix.cpp
g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c basicGRN.cpp
g++ -std=c++11 -c test_grn.cpp
#g++ -o matrix test_matrix.o GXMatrix.o
g++ -o testgrn basicGRN.o test_grn.o VertexSystem.o
./testgrn
