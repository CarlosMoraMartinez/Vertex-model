
import sys
import os
import argparse

import numpy as np
import pandas as pd

from collections import namedtuple

CELLTYPE_PAR_BEGIN = '>>'
CELLTYPE_PAR_END = '<'
CELL_TYPE_PRODUCT = '%'
PAR_BEGIN = '>'

CELLTYPE_SEP = ':'
DISJOINT_SEP = ','
SEQ_SEP = ';'
RANDOM_VALUE = '?'

OUTPUT_CELLTYPE_BEGIN = '>>'
OUTPUT_CELLTYPE_END = '<'
OUTPUT_PAR_BEGIN = '>'
OUTPUT_EXTENSION = '.vp'

MASK_SYMBOL = '{-}'

CELL_PARAM_TYPE = 1
NUM_PARAM_TYPE = 0

Param = namedtuple("Param", "file name type cell array_index value")

def getNum(string):
    if('.' in string):
        return float(string)
    else: 
        return int(string)

class ParamModelFactory:    
    """
    This class holds a the template for a set of parameters and produces a class translates them from chromosome to simulation parameters
    """
    def __init__(self, fnames):
        self.size = 0
        self.fnames = fnames
        self.params = []
        self.num_params = 0
        self.mask = []
        for i, f in enumerate(fnames):
            self.getFileParams(i, f)
    def getFileParams(self, f_index, fname):      
        rawstrings = self.readFile(fname)
        i = 0
        mascara = False
        while(i < len(rawstrings)):
            if(rawstrings[i].startswith(CELLTYPE_PAR_BEGIN)):
                parname = rawstrings[i].strip(CELLTYPE_PAR_BEGIN)
                i+=1
                cell = 0
                while(not rawstrings[i].startswith(CELLTYPE_PAR_END)):
                    parval = rawstrings[i].split(CELLTYPE_SEP)[1]
                    if(MASK_SYMBOL in parval):
                        mascara = True
                        parval = parval.strip(MASK_SYMBOL)
                    else:
                        mascara = False
                    par = Param(f_index, parname, CELL_PARAM_TYPE, cell, self.num_params, getNum(parval) )
                    i+=1   
                    cell += 1
                    self.num_params+=1
                    self.params.append(par)
                    self.mask.append(mascara)
                i+=1           
            elif(rawstrings[i].startswith(PAR_BEGIN)):
                parval = rawstrings[i + 1]
                if(MASK_SYMBOL in parval):
                    mascara = True
                    parval = parval.strip(MASK_SYMBOL)
                else:
                    mascara = False
                par = Param(f_index, rawstrings[i].strip(PAR_BEGIN), NUM_PARAM_TYPE, -1, self.num_params, getNum(parval))
                self.params.append(par)
                self.mask.append(mascara)
                self.num_params+=1
                i+=2  
            else:
                print("Error: line not recognized")
                i+=1
    def readFile(self, fname):
        f = open(fname, 'r')
        data = [x.strip() for x in f if x[0] != '#' and not x.isspace()]
        f.close()
        return data
    def __str__(self):
        s = '\n'.join([str(i) for i in self.params])
        return s
    def getParamClass(self):
        def __init__(self, values=None, ind=0):
            """
            If no values, uses values from Model. Otherwise, values can be either a list of tuples or a numpy array            
            """

            assert(len(values) == self.__class__.num_params)
            self.ind = ind
            if(values is None):
                self._chromosome = np.array([p.value for p in self.__class__.model])  
            else:
                self._chromosome = values if(type(values) is np.ndarray) else np.array([p.value for p in values])  
        def __str__(self):
            s= self.ind + "\n".join(["*** FILE: " + f +  " ***\n" + '\n'.join(d) for f, d in zip(self.paramlist.keys(), self.paramlist.values())])
            return s
        def __repr__(self):
            self.__str__()
        @property
        def chromosome(self):
            return self._chromosome
        @chromosome.setter
        def chromosome(self, params):
            assert(len(params) == self.__class__.num_params)
            if(type(params[0]) is Param):
                self._chromosome = np.array([p.value for p in params])    
            elif(type(params) is np.ndarray): 
                self._chromosome = params  
            else: 
                raise ValueError("Error initializing param set") 
        @property
        def paramlist(self):
            globalList = dict([(i, []) for i in self.__class__.fnames])
            last_cell = False
            current_file = -1
            for p in self.__class__.model:
                if(p.file != current_file):
                    if(last_cell):
                        plist.append(OUTPUT_CELLTYPE_END)
                        last_cell = False
                    current_file = p.file
                    plist = globalList[ self.__class__.fnames[current_file] ]
                if(p.type == NUM_PARAM_TYPE):
                    if(last_cell):
                        plist.append(OUTPUT_CELLTYPE_END)
                        last_cell = False
                    plist.append(OUTPUT_PAR_BEGIN + p.name)
                    plist.append(str(self._chromosome[p.array_index]))
                else:
                    if(p.cell == 0):
                        last_cell = True
                        plist.append(OUTPUT_CELLTYPE_BEGIN + p.name)
                        plist.append(str(p.cell) + CELLTYPE_SEP + str(self._chromosome[p.array_index]))
                    else:
                        plist.append(str(p.cell) + CELLTYPE_SEP + str(self._chromosome[p.array_index])) 
            if(last_cell):
                plist.append(OUTPUT_CELLTYPE_END)
            return globalList
        newclass = type("Organism", (), dict(model=self.params, mask=self.mask, num_params=self.num_params, fnames=self.fnames, __init__=__init__, __str__=__str__, __repr__=__repr__, chromosome=chromosome, paramlist=paramlist))
        return newclass


def main():
    a=ParamModelFactory(['test2steps_1a.vp', 'test2steps_1b.vp'])
    b=a.getParamClass()
    binst = b(a.params)
    binst2 = b(np.random.rand(b.num_params))
    print("TEST 1:")
    print(binst.chromosome)
    print(binst)
    print(binst2.chromosome)
    print(binst2)

if(__name__ == "__main__"):
    main()



