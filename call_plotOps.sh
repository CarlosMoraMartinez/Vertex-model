

fname=$1  #"190726_test"
for i in {0..1000};
do
 python plotOps.py $fname$i
 echo "$i "

done
