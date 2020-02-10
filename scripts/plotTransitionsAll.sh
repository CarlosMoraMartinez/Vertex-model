
cd $1
wingarray=(wing2E wing2Edumpy) #(wing2C wing2Cdumpy)
params=$(ls | grep $1)

for p in $params
do
    cat $p
    cd $p
    for wing in $(ls | grep .out)
    do
        ##cat $p $wing
        python ../../src/plotRearrangementDistributions.py $wing $2
    done
    cd ../
done 

mkdir transition_plots
mv $1*/*.pdf transition_plots

cd ../
