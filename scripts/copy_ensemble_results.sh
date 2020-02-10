

cd $1
wingarray=(wing2E wing2Edumpy) #(wing2C wing2Cdumpy)
params=$(ls | grep $1)

for wing in ${wingarray[@]}
do
    cat $wing
    newdir=$wing'_all_final'
    mkdir $newdir
    for p in $params
    do
        cp ./$p/$wing/$wing'_moved_'$2'.png' ./$newdir/$p'_'$wing'.png'
        cp ./$p/$wing/$wing'_moved_.avi' ./$newdir/$p'_'$wing'.avi'

    done
done 
