

for vp in test_params2
do
    for grn in testgrn22 #testgrn4 testgrn5 testgrn6 testgrn7 testgrn8 testgrn9 testgrn10 testgrn11 testgrn12
    do
        mkdir $grn'-'$vp
        mkdir $grn'-'$vp/param_files
        mkdir $grn'-'$vp/src
        cp src/VertexSystem.cpp src/VertexSystem.h src/basicGRN.cpp src/basicGRN.h src/run_grn.cpp scripts/run_grn.sh src/GXMatrix.h src/plotOps_loop.py $grn'-'$vp/src
        cp param_files/$grn.grn param_files/$vp.vp $grn'-'$vp/param_files
        cd $grn'-'$vp

        for cond in hexgrid_s3.0_6x6_n0.3 #hex hx5_s2.2_15x20_n0.2 #hx6_s2.2_15x20_n0.2 hx8_s2.2_15x20_n0.2 hx9_s2.2_15x20_n0.2 #hx13_s2.2_8x20_n0.2 #hx9_s2.2_15x20_n0.2 hx10_s2.2_3x10_n0.2 hx11_s2.2_3x10_n0.2 hx12_s2.2_5x20_n0.2 #hx5_s2.2_15x20_n0.2 hx6_s2.2_15x20_n0.2 hx8_s2.2_15x20_n0.2
        do
            cp -r ../$cond ./
            bash src/run_grn.sh $cond $grn $vp >$cond.out 2>$cond.err &
        done
        cd ../
        #rm $grn'_'$vp/VertexSystem.cpp $grn'_'$vp/VertexSystem.h $grn'_'$vp/basicGRN.cpp $grn'_'$vp/basicGRN.h $grn'_'$vp/run_grn.cpp $grn'_'$vp/GXMatrix.h $grn'_'$vp/plotOps_loop.py
    done
done
