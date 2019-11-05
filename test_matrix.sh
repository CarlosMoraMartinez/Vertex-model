#g++ -std=c++11 -c GXMatrix.cpp
g++ -std=c++11 -c test_matrix.cpp
#g++ -o matrix test_matrix.o GXMatrix.o
g++ -o matrix test_matrix.o 
./matrix
