#$1: Num. accepted movements
#$2: Write every N accepted movements
#$3: Name of ensemble (.vp ensemble file needs to be in ./param_files) 

maxCPU=32
i=0
wingarray=(wing2E) #(wing1D wingDumpy2)

cd param_files
python ../vertex_parms_ensemble.py -i $3'.vp' -o $3
cd ../
echo $3' Param Files generated'

mkdir $3
cp VertexSystem.cpp VertexSystem.h vertex_cmm.cpp plotOps_loop.py compile.sh $3
mkdir $3/param_files
cp ./param_files/$3/* $3/param_files

nsims=$((${#wingarray[@]} * $(ls $3/param_files | wc -l)))
echo 'Total number of simulations: '$nsims

cd $3
FILES=./param_files/*

for vp in $FILES
do
    NAME=$(basename $vp .vp)
    mkdir $NAME
    cp VertexSystem.cpp VertexSystem.h vertex_cmm.cpp plotOps_loop.py compile.sh $NAME
    cd $NAME
    mkdir param_files
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
                cp -r ../../$wing ./
                bash compile.sh $wing $NAME $1 $2 >$NAME'-'$wing.out 2>$NAME'-'$wing.err &
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

