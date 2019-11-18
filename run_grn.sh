



# mkdir $1
cp VertexSystem.cpp VertexSystem.h basicGRN.cpp basicGRN.h run_grn.cpp GXMatrix.h plotOps_loop.py $1
cp $2.grn $1
cd $1


g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c basicGRN.cpp
g++ -std=c++11 -c run_grn.cpp  #test_grn.cpp
#g++ -o matrix test_matrix.o GXMatrix.o
g++ -o grn basicGRN.o run_grn.o VertexSystem.o #test_grn.o

./grn $1 $2 $3

python plotOps_loop.py $1'_moved_' $4 $5 $6
ffmpeg -framerate 10 -i $1'_moved_'%d.png -codec copy $1'_moved_'.avi

rm *.o
rm grn

#bash run_grn.sh d1_s3.0_5x6_n0.4 testgrn 0 20 -x


