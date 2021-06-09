
#1: ensemble to replot
#2: wing to replot
#3: plot from _moved_$3
#4: plot up to

cd $1
for dir in $(ls | grep $1)
do
 cd $dir/$2 
 cp /home/carmoma/vertex/Vertex-model/src/plotOps_loop.py ./
 ffmpeg -framerate $3 -i $2_moved_%d.png -codec copy $2_moved_.avi
 cd ../../
done

