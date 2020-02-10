
#Useful commands

#Print nth line of all param files in folder

for ff in $(ls *.vp);do sed '6q;d' $ff; done
sed -n '10,33p' < file.txt #prints from line 10 to 33

#
du -sh #used in folder
df -h #all memory used

#Create and run ensemble
cd param_files
python ../vertex_parms_ensemble.py -i ensemble4.vp -o ensemble7
cd ../
bash runensemble_nogrn.sh 500000000 10000000 ensemble9


#wing2C files should be in a folder named wing2C
#ensemble7_2.vp should be in ./param_files
bash compile.sh wing2C ensemble7_2 200000000 10000000
