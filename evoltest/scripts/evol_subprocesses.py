import os
import subprocess
import multiprocessing
from string import ascii_lowercase

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
        proc = os.system(command)
        indfile = individual.ind + '_' + str(int(self.nsteps/self.write_frequency)
        ind_fitness = self.fitness_evaluator_class.getFitness(indfile)
        return ind_fitness

    def printIndividualFiles(self, individual):
        pl = individual.paramlist
        filenames = []
        for i, pfile in enumerate(pl.keys()):
            fname = '_'.join(self.basename, individual.ind, VertexSimulator.alphabet[i]) + PARAM_EXT
            filenames.append(fname)
            outF = open(fname, "w")
            outF.writelines(pl[pfile])
            outF.close()
        return filenames
    def getCommandString(self, individual, paramfiles):
        return ' '.join(COMMAND_BASE, self.initial_cond, paramfiles[0], str(self.nsteps), str(self.write_frequency), paramfiles[1], individual.ind)



