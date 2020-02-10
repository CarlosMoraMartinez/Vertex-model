


#g++ -std=c++11 -c profile_test.cpp -pg
#g++ -o profiletest profile_test.o -pg
#./profiletest 
#gprof profiletest gmon.out >profile.out
#gprof -l profiletest gmon.out >profile_byline.out

#cat profile.out


g++ -std=c++11 -c profile_test.cpp
g++ -o profiletest profile_test.o
./profiletest
