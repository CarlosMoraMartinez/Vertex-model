#$1: Num. accepted movements
#$2: Write every N accepted movements
#$3: Name of .json indicating changes (.vp files indicated in it; all need to be in ./param_files) 
###
maxCPU=32
i=0
wingarray=(etournay1_moveX2) 
#small7-16df s10-3-4 bud5_3y
#strechbig1_3cpv (bud0 hex2020_1_s5.0_20x20_n0.4) #(wing2E wing2Edumpy) #(wing1D wingDumpy2) #wing2F wing2Edumpy #bud3b bud2 iso20_3_df gr1

cd param_files
python3 ../src/vertex_parms_ensemble.py -s $3 -o $3 -i $3
cd ../
echo $3' Param Files generated'

mkdir $3
mkdir $3/src
mkdir $3/param_files
cp src/VertexSystemMinimal.cpp src/VertexSystem.h src/vertex_cmm.cpp src/plotOps_loop.py scripts/compileMinimalTemp.sh $3/src

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
    cp src/VertexSystemMinimal.cpp src/VertexSystem.h src/vertex_cmm.cpp src/plotOps_loop.py src/compileMinimalTemp.sh $NAME/src
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
                bash src/compileMinimalTemp.sh $wing $NAME $1 $2 >$NAME'-'$wing.out 2>$NAME'-'$wing.err &
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

