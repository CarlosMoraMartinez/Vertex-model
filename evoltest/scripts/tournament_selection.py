import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

import params_to_vector
from params_to_vector import ParamModelFactory

import evol_getWingBorder
from evol_getWingBorder import WingFitnessEvaluator

import evol_subprocesses
from evol_subprocesses import VertexSimulator





class TournamentPopulation:
    
    """
    This class contains the core of the tournament selection algorithm.
    Holds a population of customized objects (Organisms) and applies a function (simulation) to them
    in order to get their fitness.
    The fittest organisms reproduce sexually
    """
    
    def __init__(self, organism_class, simulator, startingOrg, N, K, mutation_rate, max_generations=100, error=0.01):
        self.organism_class = organism_class
        self.simulator = simulator
        self.N = N #Number of individuals in a population
        self.K = K #Selection is stronger when this parameter is higher
        self.num_genes = len(organism_class.model) #Number of parameters in an individual
        self.mutation_rate = mutation_rate  #You can adjust this dynamically; low mutation rate in the end of the simulation will allow you to fine-tune parameters
        self.max_generations = max_generations
        self.error = error
        self.generation = 0
        self.individuals = self.init_pop(startingOrg)
        self.fitness = np.zeros(self.num_genes)
        self.current = 0
    def __getitem__(self, i):
        return self.individuals[i]
    def __iter__(self):
        self.current = 0
        return self
    def __next__(self):
        if(self.current < self.N):
            ind = self[self.current]
            self.current += 1
            return ind
        else:
            raise StopIteration()
    def init_pop(org):
        chrom = org.chromosome
        orgs = []
        for i in range(self.N):
            new_chromosome = np.array(chrom, copy=True)
            for i in range(self.num_genes):
                if(not self.organism_class.mask[i]):
                    new_chromosome[i] += new_chromosome[i]*np.random.uniform(low=-1*self.mutation_rate, high=self.mutation_rate) 
            orgs.append(self.organism_class(new_chromosome))
    def tournamentSelection(self):        
        """
        Tournament selection algorithm

        RETURNS:
        Mean fitness of the population
        """        
        self.fitness = np.zeros(self.num_genes)
        mean_fitness = [np.mean(self.fitness)]

        for gen in range(self.max_generations):
            self.fitness = self.simulator.par_simulate(self.individuals)  ## This step can be easily parallelized with multiprocessing library
            self.generation += 1 #needs to be before self.getNewIndividual 
            self.individuals = getNewIndividuals() 
            mean_fitness.append(np.mean(self.fitness))
            if(self.generation %1000 == 0):
                print(self.generation, ': ', mean_fitness[-1])
            if(1 - mean_fitness[-1] < self.error):
                break
        return mean_fitness
    def getNewIndividuals(self):
        next_gen=[]
        for i in range(self.N):
            participants = np.random.choice(range(self.N), self.K, replace = False) #Select K random individuals
            parent_index1 = participants[np.argmax(self.fitness[participants])] #The best one will be parent 1
            participants = np.random.choice(range(self.N), self.K, replace = False) #The same for parent 2
            parent_index2 = participants[np.argmax(self.fitness[participants])]
            new_ind = self.combine2individuals(self[parent_index1], self[parent_index2], i)
            next_gen.append(new_ind) #Mix individuals
        return np.array(next_gen)
    def combine2individuals(self, org1, org2, ind=0):
        #prop = [np.mean([org1[gene], org2[gene]])for gene in range(self.num_genes)]
        new_chromosome = [np.random.choice([org1[gene], org2[gene]]) for gene in range(self.num_genes)]
        for i in range(len(new_chromosome)):
            if(not self.organism_class.mask[i]):
                new_chromosome[i] += new_chromosome[i]*np.random.uniform(low=-1*self.mutation_rate, high=self.mutation_rate) 
        return self.organism_class(np.array(new_chromosome), '_'.join(self.generation, ind))




parser = argparse.ArgumentParser(description='Process arguments.')

#PARAMETERS AND INITIAL CONDITIONS (FILE NAMES)
parser.add_argument('-m', '--StartingParams', metavar='startingparams', type=str, default='test2steps_1a.vp,test2steps_1b.vp',
                    help='Organism parameters to use to initialize population')
parser.add_argument('-s', '--StartingWing', metavar='startingwing', type=str, default = 'bud2', 
                    help='Initial condition of simulation (same in all individuals)')
parser.add_argument('-t', '--TargetShapeFile', metavar='targetshapefile', type=str, default = 'contour1', 
                    help='Target shape against which error is measured')

parser.add_argument('-r', '--MutationRate', metavar='mutation_rate', type=float, default = 0.01,
                    help='Mutation rate (recommended: around 0.01)')
parser.add_argument('-n', '--PopulationSize', metavar='N', type=int,  default = 30, 
                    help='Population size N (recommended: 100)')
parser.add_argument('-k', '--K', metavar='K', type=int, default = 7, 
                    help='K parameter for tournament selection (recommended: 40 for N = 100)')
parser.add_argument('-x', '--MaxGen', metavar='x', type=int, default = 100, 
                    help='Max number of generations')
parser.add_argument('-e', '--Error', metavar='e', type=float, default = 100, 
                    help='Minimal error threshold')
parser.add_argument('-c', '--Cores', metavar='cores', type=int, default = 24, 
                    help='Number of cores used by multiprocessing')
#VERTEX PARAMETERS
parser.add_argument('-o', '--SimName', metavar='simname', type=str, default='simul0',
                    help='Name given to organisms')
parser.add_argument('-w', '--VertexNumSteps', metavar='w', type=int, default = 100, 
                    help='Max number of generations')
parser.add_argument('-z', '--VertexWriteFreq', metavar='z', type=int, default = 100, 
                    help='Max number of generations')


#python tournament_selection -o evotest1 -m 'test2steps_1a.vp,test2steps_1b.vp' -s 'bud2' -t 'contour1' -r 0.2 -n 30 -k 7 -x 100 -c 30 -w 50000000 -z 1000000
def main():
    args = parser.parse_args()
    
    simname = args.SimName
    mutation_rate = args.MutationRate
    N = args.PopulationSize
    K = args.K
    max_generations = args.MaxGen
    error = args.Error
    ncores = args.Cores
    vertex_numsteps = args.VertexNumSteps
    vertex_writefreq = args.VertexWriteFreq
   
    starting_wing = args.StartingWing
    initial_params = args.StartingParams.split(',')
    targetShapeFile = args.TargetShapeFile
    
    orgGenerator=ParamModelFactory(initial_params)
    org_class=orgGenerator.getParamClass()
    starting_org = org_class() #By default uses input model 

    fitness_evaluator = WingFitnessEvaluator
    fitness_evaluator.setTargetShape(targetShapeFile) 
    simulator = VertexSimulator(simname, ncores, starting_wing, fitness_evaluator, vertex_numsteps, vertex_writefreq)


    pop = TournamentPopulation(organism_class = org_class, \
                               simulator = simulator, \
                               startingOrg = starting_org, \
                               N = N, \
                               K = K, \
                               mutation_rate = mutation_rate, \
                               max_generations = max_generations, \
                               error = error)

    time_fitness = pop.tournamentSelection()
    print("Ending simulation at generation %d with error %.4f"%(pop.generation, 1 - time_fitness[-1]) )
    print("Best organism was: %d"%(np.argmax(pop.fitness)))

    plt.plot(range(len(time_fitness)), time_fitness)
    plt.xlabel("generation")
    plt.ylabel("mean population error")
    plt.show()

if __name__== "__main__":
  main()
