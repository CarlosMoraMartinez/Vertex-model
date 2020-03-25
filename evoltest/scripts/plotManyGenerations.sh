

#1: Starting generation
#2: Final generation
#3: First individual
#4: Last individual

for (( generation=$1; generation<=$2; generation++ ))
do
    for (( individual=$3; individual<=$4; individual++ ))
    do
        cat 'plotting: '$generation'-'$individual'_moved_'
        python plotOps_loop.py -i $generation'-'$individual'_moved_' -s 1 -e 1000
    done
done
