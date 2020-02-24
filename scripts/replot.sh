cd $1
wingarray=(bud1) #(wing2E wing2Edumpy) #(wing2C wing2Cdumpy)
params=$(ls | grep $1)

for wing in ${wingarray[@]}
do
    cat $wing
    for p in $params
    do
        cd $p/$wing
        cp ../../../src/plotOps_loop.py ./
        python plotOps_loop.py -i $wing'_moved_' -s $2 -e $3 -l $4
        ffmpeg -framerate 10 -i $wing'_moved_'%d.png -codec copy $wing'_moved_b'.avi
        cd ../../
    done
done 
