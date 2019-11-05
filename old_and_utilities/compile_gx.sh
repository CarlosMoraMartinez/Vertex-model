

g++ -std=c++11 -c GeneExpression.cpp
g++ -std=c++11 -c test_gx.cpp
g++ -o gx_test test_gx.o GeneExpression.o
./gx_test 
