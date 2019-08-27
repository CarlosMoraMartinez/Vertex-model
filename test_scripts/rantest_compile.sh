g++ -std=c++11 -c randomtest.cpp
g++ -o rantest randomtest.o
./rantest "$1" >"test_random_1.tsv"

python plotRandom.py

rm randomtest.o
rm rantest
rm "test_random_1.tsv"
