
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
EQUAL_CELLTYPES_SYMBOL = '{=}'

CELL_PARAM_TYPE = 1
NUM_PARAM_TYPE = 0

Param = namedtuple("Param", "file name type cell array_index value mask all_equal")

def getNum(string):
    if(type(string) is not str):
        return string
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
        self.chrom_size = 0
        self.remove_from_array = 0
        for i, f in enumerate(fnames):
            self.getFileParams(i, f)
    def getFileParams(self, f_index, fname):      
        rawstrings = self.readFile(fname)
        i = 0
        mascara = False  #True if parameter can't evolve
        equal_in_types=False #True if all cell types must keep the same value (implies mask)
        while(i < len(rawstrings)):
            if(rawstrings[i].startswith(CELLTYPE_PAR_BEGIN)):
                parname = rawstrings[i].replace(CELLTYPE_PAR_BEGIN, '')
                i+=1
                cell = 0
                while(not rawstrings[i].startswith(CELLTYPE_PAR_END)):
                    parval = rawstrings[i].split(CELLTYPE_SEP)[1]
                    if(MASK_SYMBOL in parval):
                        mascara = True
                        parval = parval.replace(MASK_SYMBOL, '')
                        self.remove_from_array +=1
                    else:
                        mascara = False
                    if(EQUAL_CELLTYPES_SYMBOL in parval):
                        equal_in_types = True
                        #substract_index = cell
                        if(cell > 0):
                            mascara = True
                            self.remove_from_array +=1
                            parval = self.params[-1].value
                        else:
                            parval = parval.replace(EQUAL_CELLTYPES_SYMBOL, '')
                    else:
                        equal_in_types = False
                        #substract_index = 0
                    array_ind = -1 if mascara and not equal_in_types else self.num_params - self.remove_from_array
                    par = Param(f_index, parname, CELL_PARAM_TYPE, cell, array_ind, getNum(parval), mascara, equal_in_types)#- cell is useful for equal parameters
                    i+=1   
                    cell += 1
                    self.num_params+=1
                    self.params.append(par)
                i+=1           
            elif(rawstrings[i].startswith(PAR_BEGIN)):
                parval = rawstrings[i + 1]
                if(MASK_SYMBOL in parval):
                    mascara = True
                    self.remove_from_array +=1
                    parval = parval.replace(MASK_SYMBOL, '')
                else:
                    mascara = False
                equal_in_types = False
                array_ind =  -1 if mascara else self.num_params - self.remove_from_array
                par = Param(f_index, rawstrings[i].replace(PAR_BEGIN, ''), NUM_PARAM_TYPE, -1, array_ind, getNum(parval), mascara, equal_in_types)
                self.params.append(par)
                self.num_params+=1
                i+=2  
            else:
                print("Error: line not recognized")
                i+=1
        self.chrom_size = self.num_params - self.remove_from_array
    def readFile(self, fname):
        f = open(fname, 'r')
        data = [x.strip() for x in f if x[0] != '#' and not x.isspace()]
        f.close()
        return data
    def __str__(self):
        s = '\n'.join([str(i) for i in self.params])
        return s
    def getParamClass(self):
        def __init__(self, values=None, ind='0'):
            """
            If no values, uses values from Model. Otherwise, values can be either a list of tuples or a numpy array            
            """
            self.ind = ind
            if(values is None):
                self._chromosome = np.array([p.value for p in self.__class__.model if not p.mask])  
            else:
                #assert(len(values) == self.__class__.num_params)
                self._chromosome = values if(type(values) is np.ndarray) else np.array([p.value for p in values if not p.mask])  
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
            #assert(len(params) == self.__class__.num_params)
            if(type(params[0]) is Param):
                self._chromosome = np.array([p.value for p in params if not p.mask])    
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
                if(p.array_index >= 0):
                    pvalue = str(self._chromosome[p.array_index])
                else:
                    pvalue = str(p.value)
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
                    plist.append(pvalue)
                else:
                    if(p.cell == 0):
                        if(last_cell):
                            plist.append(OUTPUT_CELLTYPE_END)
                        last_cell = True
                        plist.append(OUTPUT_CELLTYPE_BEGIN + p.name)
                        plist.append(str(p.cell) + CELLTYPE_SEP + pvalue)
                    else:
                        plist.append(str(p.cell) + CELLTYPE_SEP + pvalue) 
            if(last_cell):
                plist.append(OUTPUT_CELLTYPE_END)
            return globalList
        newclass = type("Organism", (), dict(model=self.params, chrom_size=self.chrom_size, num_params=self.num_params, fnames=self.fnames, __init__=__init__, __str__=__str__, __repr__=__repr__, chromosome=chromosome, paramlist=paramlist))
        return newclass


def main():
    a=ParamModelFactory(['param_files/test2steps_1a.vp', 'param_files/test2steps_1b.vp'])
    for p in a.params:
        print(p)
    print('***')
    b=a.getParamClass()
    binst = b(a.params)
    bc = binst.chromosome
    for p in a.params:
        print(p)
        print("MASK IND" if p.mask and not p.all_equal else bc[p.array_index])
        
    binst2 = b(np.random.rand(b.chrom_size))
    print("TEST 1:")
    print(binst.chromosome)
    print(binst)
    print(binst2.chromosome)
    print(binst2)
    print(binst2.chrom_size, binst2.num_params, binst2.chromosome.shape, binst.chromosome.shape)

if(__name__ == "__main__"):
    main()



