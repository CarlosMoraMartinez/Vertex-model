#$1: Name of file without extension (a .json file that needs to be in ./param_files) 

maxCPU=10
i=0
wingarray=(bud5_3ctfat3) 

cd param_files
python ../src/vertex_parms_ensemble.py -s $1 #must be a .json file
cd ../
echo $1' Param Files generated'

mkdir $1
mkdir $1/src
mkdir $1/param_files
cp src/VertexSystem.cpp src/VertexSystem.h src/vertex_3steps.cpp src/plotOps_loop.py scripts/compile_3steps.sh $1/src

cp ./param_files/$1/* $1/param_files
cp ./param_files/$2'.vp' ./param_files/$3'.vp' ./param_files/$4'.vp' $1/param_files

nsims=$(( ${#wingarray[@]} * ($(ls $/param_files | wc -l) - 1) ))
echo 'Total number of simulations: '$nsims

cd $1
echo 'cd ../;bash ./scripts/runensemble_3steps.sh ' $@ > RUN_AGAIN.sh

FILES1=./param_files/$2*.vp
FILES2=./param_files/$3*.vp
FILES3=./param_files/$4*.vp

for vp in $FILES1
do
    NAME=$(basename $vp .vp)
    mkdir $NAME
    mkdir $NAME/src
    mkdir $NAME/param_files
    cp src/VertexSystem.cpp src/VertexSystem.h src/vertex_3steps.cpp src/plotOps_loop.py src/compile_3steps.sh $NAME/src
    cd $NAME

    cp '../param_files/'$NAME'.vp' '../param_files/'$3'.vp' '../param_files/'$4'.vp' param_files
    echo 'PHASE 1: PARAM_FILE: '$NAME

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
                bash src/compile_3steps.sh $wing $NAME $3 $4 >$NAME'-'$wing.out 2>$NAME'-'$wing.err &
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

for vp in $FILES2
do
    NAME=$(basename $vp .vp)
    mkdir $NAME
    mkdir $NAME/src
    mkdir $NAME/param_files
    cp src/VertexSystem.cpp src/VertexSystem.h src/vertex_3steps.cpp src/plotOps_loop.py src/compile_3steps.sh $NAME/src
    cd $NAME

    cp '../param_files/'$NAME'.vp' '../param_files/'$2'.vp' '../param_files/'$4'.vp' param_files
    echo 'PHASE 2: PARAM_FILE: '$NAME

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
                bash src/compile_3steps.sh $wing $2 $NAME $4 >$NAME'-'$wing.out 2>$NAME'-'$wing.err &
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

for vp in $FILES3
do
    NAME=$(basename $vp .vp)
    mkdir $NAME
    mkdir $NAME/src
    mkdir $NAME/param_files
    cp src/VertexSystem.cpp src/VertexSystem.h src/vertex_3steps.cpp src/plotOps_loop.py src/compile_3steps.sh $NAME/src
    cd $NAME

    cp '../param_files/'$NAME'.vp' '../param_files/'$2'.vp' '../param_files/'$3'.vp' param_files
    echo 'PHASE 3: PARAM_FILE: '$NAME

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
                bash src/compile_3steps.sh $wing $2 $3 $NAME >$NAME'-'$wing.out 2>$NAME'-'$wing.err &
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
