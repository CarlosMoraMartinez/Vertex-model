

for vp in test_params7 test_params12 test_params13 test_params14 test_params15
do
    mkdir $vp
    cp VertexSystem.cpp VertexSystem.h vertex_cmm.cpp plotOps_loop.py compile.sh $vp
    cd $vp

    for wing in wing1D wingDumpy wingDumpy2
    do
        cp -r ../$wing ./
        mkdir param_files
        cp '../param_files/'$vp'.vp' param_files
        echo $vp ' ' $wing
        bash compile.sh $wing $vp $1 $2 >$vp'-'$wing.out 2>$vp'-'$wing.err &
    done
    cd ../
done
