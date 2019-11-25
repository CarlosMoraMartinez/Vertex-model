

for vp in test_params4
do
    for grn in testgrn21 testgrn22 testgrn23 #testgrn4 testgrn5 testgrn6 testgrn7 testgrn8 testgrn9 testgrn10 testgrn11 testgrn12
    do
        mkdir $grn'-'$vp
        mkdir $grn'-'$vp/param_files
        cp VertexSystem.cpp VertexSystem.h basicGRN.cpp basicGRN.h run_grn.cpp run_grn.sh GXMatrix.h plotOps_loop.py $grn'-'$vp
        cp param_files/$grn.grn param_files/$vp.vp $grn'-'$vp/param_files
        cd $grn'-'$vp

        for cond in  hx5_s2.2_15x20_n0.2 hx6_s2.2_15x20_n0.2 hx8_s2.2_15x20_n0.2 hx9_s2.2_15x20_n0.2 #hx13_s2.2_8x20_n0.2 #hx9_s2.2_15x20_n0.2 hx10_s2.2_3x10_n0.2 hx11_s2.2_3x10_n0.2 hx12_s2.2_5x20_n0.2 #hx5_s2.2_15x20_n0.2 hx6_s2.2_15x20_n0.2 hx8_s2.2_15x20_n0.2
        do
            cp -r ../$cond ./
            bash run_grn.sh $cond $grn $vp >$cond.out 2>$cond.err &
        done
        cd ../
        #rm $grn'_'$vp/VertexSystem.cpp $grn'_'$vp/VertexSystem.h $grn'_'$vp/basicGRN.cpp $grn'_'$vp/basicGRN.h $grn'_'$vp/run_grn.cpp $grn'_'$vp/GXMatrix.h $grn'_'$vp/plotOps_loop.py
    done
done
