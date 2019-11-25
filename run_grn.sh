



# mkdir $1
cp VertexSystem.cpp VertexSystem.h basicGRN.cpp basicGRN.h run_grn.cpp GXMatrix.h plotOps_loop.py $1
cp param_files/$2.grn $1
cp param_files/$3.vp $1
cd $1


g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c basicGRN.cpp
g++ -std=c++11 -c run_grn.cpp  #test_grn.cpp
#g++ -o matrix test_matrix.o GXMatrix.o
g++ -o grn basicGRN.o run_grn.o VertexSystem.o #test_grn.o

./grn $1 $2 $3

python plotOps_loop.py -i $1'_moved_' -i2 $2'_moved_'  -g '8'
ffmpeg -framerate 10 -i $1'_moved_'%d.png -codec copy $1'_moved_'.avi
ffmpeg -framerate 10 -i $1'_moved_'%d'g8'.png -codec copy $1'g8'.avi
ffmpeg -i $1'_moved_'.avi -i $1'g8'.avi -filter_complex hstack=inputs=2  stacked.avi

rm *.o
rm grn

#bash run_grn.sh d1_s3.0_5x6_n0.4 testgrn 0 20 -x


