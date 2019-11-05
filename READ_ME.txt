Carlos Mora, August 2019

Vertex model to simulate Drosophila wing development


- REQUIREMENTS

The program is written in C++ 11 (some language features that are used in the program are not available in older versions)
Compiler used is gcc version 7.2.0 
OS used is Ubuntu 16.04.4 LTS

python scripts for plotting and generating hexagonal grids require python3 and use the following libraries:

numpy
matplotlib
shapely
descartes

Additionally, the compile.sh script uses the ffmpeg utility.



- COMPILE AND RUN:

Use:
bash compile.sh <simulation_name> <max_accepted_moves> <number of accepted moves between outputs> <add vertex numbers to plot (-n) or not (otherwise)>

Example:

bash compile.sh winglike2_s3.25_20x20_n0.4 10000000 100000 -x 


The first three arguments are passed to the vertex model itself. Fourth argument is passed to the plotOps_loop.py script (which comes from Jhon, with slight modifications).

- The first one is the prefix used for input and output files. 
- The second one is the number of accepted moves before simulation stops
- Third argument: every X moves, the program will write the state of the system: .cells, .points and .spr files, necessary to plot, and a .out file with all the information of cells, vertices and edges.


To create a new hexagonal grid, use the make_hexagonal_grid.py script. Usage:

python make_hexagonal_grid.py -o winglike2 -s 3.25 -r 20 -c 20 -n 0.4 -t "10:20;0:3;;" -p ";18,19;18,19;0:20" -l 5 -g 15 -v "5,8, 11, 15"

For explanation of arguments:

python make_hexagonal_grid.py --h

