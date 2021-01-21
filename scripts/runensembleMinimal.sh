#$1: Num. accepted movements
#$2: Write every N accepted movements
#$3: Name of ensemble (.vp ensemble file needs to be in ./param_files) 

maxCPU=10
i=0
wingarray=(etournay1_strings8b) #etournay1_strings3
#etournay1_3cpv_3 etournay1_3cpv_4
# etournay1_strings3 etournay1_strings4
#etournay1_vein1d etournay1_vein1a etournay1_nosprings etournay1_nosprings3 etournay1_strings etournay1_nosprings 
#small7-16df s10-3-4 bud5_3y etournay1_unmoveX7
#strechbig1_3cpv (bud0 hex2020_1_s5.0_20x20_n0.4) #(wing2E wing2Edumpy) #(wing1D wingDumpy2) #wing2F wing2Edumpy #bud3b bud2 iso20_3_df gr1

cd param_files
python ../src/vertex_parms_ensemble.py -i $3'.vp' -o $3
cd ../
echo $3' Param Files generated'

mkdir $3
mkdir $3/src
mkdir $3/param_files
cp src/VertexSystemMinimal.cpp src/VertexSystem.h src/vertex_cmm.cpp src/plotOps_loop.py scripts/compileMinimal.sh $3/src

cp ./param_files/$3/* $3/param_files

nsims=$(( ${#wingarray[@]} * ($(ls $3/param_files | wc -l) - 1) ))
echo 'Total number of simulations: '$nsims

cd $3
echo 'cd ../;bash ./scripts/runensemble_nogrn.sh ' $@ > RUN_AGAIN.sh

FILES=./param_files/*.vp

for vp in $FILES
do
    NAME=$(basename $vp .vp)
    mkdir $NAME
    mkdir $NAME/src
    mkdir $NAME/param_files
    cp src/VertexSystemMinimal.cpp src/VertexSystem.h src/vertex_cmm.cpp src/plotOps_loop.py src/compileMinimal.sh $NAME/src
    cd $NAME

    cp '../param_files/'$NAME'.vp' param_files
    echo 'PARAM_FILE: '$NAME

    for wing in "${wingarray[@]}"
    do
        echo 'WING: '$wing
        while [ $i -lt 1 ]
        do
            modelrunning=$(ps -aef | pgrep vertex | wc -l)  #how many CPUs are running?
            plotting=$(ps -aef | pgrep plotOps_loop.py | wc -l)
            cpuRUNNING=$(($modelrunning + $plotting))
            echo '  Counting CPUs...: '$cpuRUNNING
            if [ $cpuRUNNING -lt $maxCPU ]
            then
                echo '      CPUs lower than '$maxCPU
                echo '      sending job: '$NAME ' ' $wing
                cp -r ../../initial_conds/$wing ./
                bash src/compileMinimal.sh $wing $NAME $1 $2 >$NAME'-'$wing.out 2>$NAME'-'$wing.err &
                i=$(($i + 1))
                sleep 3
            else
                sleep 30
            fi
        done
        i=0
    done
    cd ../
done

