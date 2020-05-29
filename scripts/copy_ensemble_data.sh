cd $1
wingarray=(wing2E)  #(iso20-3-60dfs) #(wing2E wing2Edumpy) #(wing2C wing2Cdumpy)bud2
params=$(ls | grep $1)

for wing in ${wingarray[@]}
do
    cat $wing
    newdir=$wing'_landmarks'
    mkdir $newdir
    for p in $params
    do
        cp ./$p/$wing/$wing'_moved_'$2'.edges' ./$newdir/$p'_'$wing'.edges'
        cp ./$p/$wing/$wing'_moved_'$2'.celltab' ./$newdir/$p'_'$wing'.celltab'
        cp ./$p/$wing/$wing'_moved_'$2'.points' ./$newdir/$p'_'$wing'.points'


    done
done 
cd $newdir
python3 ../../src/findLandmarks.py $1 $3

