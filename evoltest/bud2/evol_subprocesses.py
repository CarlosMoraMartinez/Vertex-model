import os
import subprocess

for i in range(10):   
    print("starting", i)
    proc = os.system("./vertex bud2 test2steps_1a 1000000 100000 test2steps_1b")
    print("finished: ", i)


