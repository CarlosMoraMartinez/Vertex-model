import os
import subprocess
import multiprocessing
from string import ascii_lowercase

import numpy as np

PARAM_EXT  ='.vp'
COMMAND_BASE = "./vertex"

class VertexSimulator:
    alphabet = ascii_lowercase
    def __init__(self, basename, initial_cond, fitness_evaluator, ncores, nsteps, write_frequency):
        self.basename = basename
        self.initial_cond = initial_cond
        self.fitness_evaluator_class = fitness_evaluator
        self.ncores = ncores
        self.nsteps = nsteps
        self.write_frequency = write_frequency

    def par_simulate(self, individuals):
        return np.array([self.simulate_single(i) for i in individuals])
    def simulate_single(self, individual):
        paramfiles = self.printIndividualFiles(individual)
        command  = self.getCommandString(individual, paramfiles)
        print("COMMAND: ", command)
        try:
            proc = os.system(command)
            indfile = individual.ind + '_moved_' + str(2*int(self.nsteps/self.write_frequency))
            ind_fitness = self.fitness_evaluator_class.getFitness(indfile)
        except:
            ind_fitness = 0
        return ind_fitness

    def printIndividualFiles(self, individual):
        pl = individual.paramlist
        filenames = []
        for i, pfile in enumerate(pl.keys()):
            letter = pfile.strip(PARAM_EXT)[-1]
            fname = '_'.join([self.basename, individual.ind, letter + PARAM_EXT]) #VertexSimulator.alphabet[i]])
            print('OUT FNAME: ', fname)
            filenames.append(fname)
            outF = open(fname, "w")
            outF.writelines('\n'.join(pl[pfile]))
            outF.close()
        return [f.strip(PARAM_EXT) for f in filenames]
    def getCommandString(self, individual, paramfiles):
        paramfiles.sort() #Since they are hash keys, make sure then are in order
        print(paramfiles)
        return ' '.join([COMMAND_BASE, self.initial_cond, str(paramfiles[0]), str(self.nsteps), str(self.write_frequency), str(paramfiles[1]), str(individual.ind)])



