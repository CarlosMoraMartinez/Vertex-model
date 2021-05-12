

cd $1
wingarray=(etournay1_strings15 etournay1_strings10) 
#(wing2E wing2Edumpy) #(wing2C wing2Cdumpy)bud2 budsmall strechsmall_3cpv etournay1_nosprings etournay1_strings3 etournay1_3cpv_3 21grid5
#square1
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
	cp ./$p/$wing/$wing'_moved_'$2'.ptab' ./$newdir/$p'_'$wing'.ptab'
        cp ./$p/$wing/$wing'_moved_'$2'.spr' ./$newdir/$p'_'$wing'.spr'
        cp ./$p/$wing/$wing'_moved_'$2'.sprtab' ./$newdir/$p'_'$wing'.sprtab'
    done
done 
