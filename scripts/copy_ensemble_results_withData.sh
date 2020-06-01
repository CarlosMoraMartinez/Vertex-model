

cd $1
wingarray=(strechsmall_3cpv) 
#(wing2E wing2Edumpy) #(wing2C wing2Cdumpy)bud2 budsmall

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
        cp ./$p/$wing/$wing'_moved_'$2'.points' ./$newdir/$p'_'$wing'.points'
        cp ./$p/$wing/$wing'_moved_'$2'.cells' ./$newdir/$p'_'$wing'.cells'
        cp ./$p/$wing/$wing'_moved_'$2'.edges' ./$newdir/$p'_'$wing'.edges'
        cp ./$p/$wing/$wing'_moved_'$2'.celltab' ./$newdir/$p'_'$wing'.celltab'
        cp ./$p/$wing/$wing'_moved_'$2'.spr' ./$newdir/$p'_'$wing'.spr'
    done
done 
