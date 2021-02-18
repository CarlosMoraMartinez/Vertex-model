
#1: ensemble to replot
#2: wing to replot
#3: plot from _moved_$3
#4: plot up to

cd $1
for dir in $(ls | grep $1)
do
 cd $dir/$2 
 cp /home/carmoma/vertex/Vertex-model/src/plotOps_loop.py ./
 python plotOps_loop.py -i $2'_moved_' -s $3 -e $4 -l 1 -f 1
 rm plotOps_loop.py
 cd ../../
done

