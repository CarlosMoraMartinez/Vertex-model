
#1: name of simulation
#2: param file wing growth
#3: param file hinge contraction
#4: wing initial condition (eg, bud2)
#5: wing to compare with (.border)

# Example:
# bash scripts/evol.sh evol1 test2steps_1a test2steps_1b bud3b target1 >evol1.out &
mkdir $1
cp ../src/VertexSystem.cpp ../src/VertexSystem.h ../src/vertex_2steps.cpp 'param_files/'$2'.vp' 'param_files/'$3'.vp' 'param_files/'$5'.border' $1
cp ./$4/$4.* $1
cd $1

g++ -std=c++11 -c VertexSystem.cpp
g++ -std=c++11 -c vertex_2steps.cpp
g++ -o vertex vertex_2steps.o VertexSystem.o

python ../scripts/tournament_selection.py -o $1 -m $2'.vp,'$3'.vp' -s $4 -t $5 -r 0.1 -n 32 -k 10 -e 0.01 -x 10000 -v 2 -c 32 -w 50000000 -z 50000000


#>bash scripts/tournament_selection.py --help
#
#usage: tournament_selection.py [-h] [-m startingparams] [-s startingwing]
#                               [-t targetshapefile] [-r mutation_rate] [-n N]
#                               [-k K] [-x x] [-e e] [-c cores] [-o simname]
#                               [-w w] [-z z]
#
#Process arguments.
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -m startingparams, --StartingParams startingparams
#                        Organism parameters to use to initialize population
#  -s startingwing, --StartingWing startingwing
#                        Initial condition of simulation (same in all
#                        individuals)
#  -t targetshapefile, --TargetShapeFile targetshapefile
#                        Target shape against which error is measured
#  -r mutation_rate, --MutationRate mutation_rate
#                        Mutation rate (recommended: around 0.01)
#  -n N, --PopulationSize N
#                        Population size N (recommended: 100)
#  -k K, --K K           K parameter for tournament selection (recommended: 40
#                        for N = 100)
#  -x x, --MaxGen x      Max number of generations
#  -e e, --Error e       Minimal error threshold
#  -c cores, --Cores cores
#                        Number of cores used by multiprocessing
#  -o simname, --SimName simname
#                        Name given to organisms
#  -w w, --VertexNumSteps w
#                        Max number of generations
#  -z z, --VertexWriteFreq z
#                       Max number of generations
#
