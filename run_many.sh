

for vp in test_params
do
    for grn in testgrn4 testgrn5 testgrn6 testgrn7 testgrn8
    do
        mkdir $grn'-'$vp
        mkdir $grn'-'$vp/param_files
        cp VertexSystem.cpp VertexSystem.h basicGRN.cpp basicGRN.h run_grn.cpp run_grn.sh GXMatrix.h plotOps_loop.py $grn'-'$vp
        cp param_files/$grn.grn param_files/$vp.vp $grn'-'$vp/param_files
        cd $grn'-'$vp

        for cond in hx5_s2.2_15x20_n0.2 hx6_s2.2_15x20_n0.2 hx7_s2.2_15x20_n0.2 hx8_s2.2_15x20_n0.2
        do
            cp -r ../$cond ./
            bash run_grn.sh $cond $grn $vp >$cond.out&
        done
        cd ../
        #rm $grn'_'$vp/VertexSystem.cpp $grn'_'$vp/VertexSystem.h $grn'_'$vp/basicGRN.cpp $grn'_'$vp/basicGRN.h $grn'_'$vp/run_grn.cpp $grn'_'$vp/GXMatrix.h $grn'_'$vp/plotOps_loop.py
    done
done
