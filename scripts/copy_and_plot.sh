

cd $1
wingarray=(etournay1_strings10) 
#(wing2E wing2Edumpy) #(wing2C wing2Cdumpy)bud2 budsmall strechsmall_3cpv etournay1_nosprings etournay1_strings3 etournay1_3cpv_3

params=$(ls | grep $1)

for wing in ${wingarray[@]}
do
    cat $wing
    newdir=$wing'_all_final'
    mkdir $newdir
    cp ../ini2plot/* $newdir
    cp ../src/plotOps_loop.py $newdir

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

   cd $newdir
   numsims=$(ls ../ |grep $1 | wc -l)
   echo "Plotting ""$numsims"" conditions"
   echo $params
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -p preferred_area_norm__4 -f 1
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -p area_ratio__4 -f 1
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -p preferred_area__4 -f 1
   python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -p perim_contract__6 -f 1
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -p area__4 -f 1
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -p perimeter__4 -f 1
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -p initial_xpos -f 1
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -c etournayfinal
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -c $wing
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims
   #python plotOps_loop.py -i $1"_*_"$wing -s 0 -e $numsims -f 1

done 
