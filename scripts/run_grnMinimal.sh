

#1 wing
#2 .grn file
#3 .vp file


# mkdir $1
cp src/VertexSystemMinimal.cpp src/VertexSystem.h src/basicGRN.cpp src/basicGRN.h src/run_grn.cpp src/GXMatrix.h src/plotOps_loop.py $1
cp param_files/$2.grn $1
cp param_files/$3.vp $1
cd $1


g++ -O2 -std=c++11 -c VertexSystemMinimal.cpp
g++ -O2 -std=c++11 -c basicGRN.cpp
g++ -O2 -std=c++11 -c run_grn.cpp  #test_grn.cpp
#g++ -o matrix test_matrix.o GXMatrix.o
g++ -o grn basicGRN.o run_grn.o VertexSystemMinimal.o #test_grn.o

./grn $1 $2 $3

python plotOps_loop.py -i $1'_moved_' -i2 $2'_moved_'  -g '8'
ffmpeg -framerate 10 -i $1'_moved_'%d.png -codec copy $1'_moved_'.avi
ffmpeg -framerate 10 -i $1'_moved_'%d'g8'.png -codec copy $1'g8'.avi
ffmpeg -i $1'_moved_'.avi -i $1'g8'.avi -filter_complex hstack=inputs=2  stacked.avi

rm *.o
rm grn

#bash run_grn.sh d1_s3.0_5x6_n0.4 testgrn 0 20 -x


