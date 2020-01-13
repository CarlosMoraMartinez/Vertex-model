import sys
import os

import numpy as np


CELLTYPE_PAR_BEGIN = '>>'
CELLTYPE_PAR_END = '<'
CELL_TYPE_PRODUCT = '%'
PAR_BEGIN = '>'

CELLTYPE_SEP = ':'
DISJOINT_SEP = ','
SEQ_SEP = ';'

OUTPUT_CELLTYPE_BEGIN = '>'
OUTPUT_CELLTYPE_END = '<'
OUTPUT_PAR_BEGIN = '>'
OUTPUT_EXTENSION = '.vp'

class paramContainer:

    def __init__(self, fname, outname):
        self.fname = fname
        self.outname = outname
        rawstrings = self.readFile(fname)
        self.params = []
        i = 0
        while(i < len(rawstrings)):
            par = []
            if(rawstrings[i].startswith(CELLTYPE_PAR_BEGIN)):
                par.append(rawstrings[i].strip(CELLTYPE_PAR_BEGIN))
                i+=1
                while(not rawstrings[i].startswith(CELLTYPE_PAR_END)):
                    par.append(rawstrings[i])
                    i+=1   
                parobj = cellTypeParam(par)
                i+=1       
                self.params.append(parobj)      
            elif(rawstrings[i].startswith(PAR_BEGIN)):
                par.append(rawstrings[i])
                par.append(rawstrings[i + 1])
                parobj = param(par)
                i+=2
                self.params.append(parobj)    
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

    def produceOutputs(self):
        from itertools import product
        iters = [iter(p) for p in self.params]
        allnames = [p.name for p in self.params]
        extended = [[val for val in p] for p in self.params]
        combos = [c for c in product(*extended)]
        num = len(combos)
        for i in range(num):
            parset = dict(zip(allnames, combos[i]))
            self.writeFile(parset, self.outname + '_' + str(i))

    def writeFile(self, parset, fname):
        s = '# ' + fname + '\n'
        for k, v in zip(parset.keys(), parset.values()):
            if(type(v) is dict):
                s += OUTPUT_CELLTYPE_BEGIN + k + '\n' 
                for k2, v2 in zip(v.keys(), v.values()):  
                    s += k2 + CELLTYPE_SEP + str(v2) + '\n'            
                s += OUTPUT_CELLTYPE_END + '\n\n' 

            else:
                s += OUTPUT_PAR_BEGIN + k + '\n' + str(v) + '\n'
        if not os.path.exists(self.outname):
            os.makedirs(self.outname)
        outfile = self.outname + '/' + fname
        f = open(outfile + OUTPUT_EXTENSION,"w+")
        f.write(s)
        f.close()


        print(parset)   
        print("\n\n\n")     

class param:
    def __init__(self, parlist):
        self.name = parlist[0].strip(PAR_BEGIN)
        self.atomic, self.value = self.processParam(parlist[1])
        self.iterstatus = 0
        self.finished_at_least_once = False

    def processParam(self, string):
        elements1 = string.split(DISJOINT_SEP)
        if(len(elements1) == 1 and not SEQ_SEP in elements1[0]):
            atomic = True
            value = self.getNum(elements1[0])
        else:
            atomic = False
            value = []
            for i in elements1:
                if(SEQ_SEP in i):
                    i = i.split(SEQ_SEP)
                    for j in np.linspace(self.getNum(i[0]), self.getNum(i[1]), self.getNum(i[2])):
                        if('.' not in i[0]):
                            j = int(j)
                        value.append(j)
                else:
                    value.append(self.getNum(i))
        return(atomic, value)

    def getNum(self, string):
        if('.' in string):
            return float(string)
        else: 
            return int(string)
    def __iter__(self):
        return self
    def __next__(self):
        if(self.finished_at_least_once):
            self.finished_at_least_once = False
            raise StopIteration
        else:
            if(self.atomic):
                item = self.value
                self.finished_at_least_once = True
            else:
                item = self.value[self.iterstatus]
                self.iterstatus += 1
                if(self.iterstatus == len(self.value)): 
                    self.iterstatus = 0
                    self.finished_at_least_once = True
        return item

    def __str__(self):
        s = OUTPUT_PAR_BEGIN + self.name + '\n'
        s += str(self.value)
        return s

               
class cellTypeParam(param):
    def __init__(self, parlist):
        self.name = parlist.pop(0)
        if(self.name.startswith(CELL_TYPE_PRODUCT)):
            self.product = True
            self.len = None
        else:
            self.product = False
            self.len = 1
        self.value = dict()
        self.atomic = dict()
        self.iterstatus = 0
        self.finished_at_least_once = False
        for p in parlist:
            p = p.split(CELLTYPE_SEP)
            at, val = self.processParam(p[1])
            if(at):
                val = [val]
            else: 
                self.len = len(val)
            self.value.setdefault(p[0], val)
            self.atomic.setdefault(p[0], at)

        if(not self.product):
            for k in self.value.keys():
                if(self.atomic[k]):
                    self.value[k] = [self.value[k][0] for i in range(self.len)]
                    self.atomic[k] = False
        self.iterkeys = None
        self.itervals = None
        self.combos = None
    def __iter__(self):
        from itertools import product
        self.iterkeys = [i for i in self.value.keys()]
        self.itervals = [i for i in self.value.values()] 
        if(self.product):
            self.combos = [i for i in product(*self.itervals)]
            self.len = len(self.combos)
        self.iterstatus = 0
        return self
    def __next__(self):
        if(self.finished_at_least_once):
            self.finished_at_least_once = False
            raise StopIteration
        else:
            if(self.product):
                cc = self.combos[self.iterstatus]
            else:
                cc = [v[self.iterstatus] for v in self.itervals]
            self.iterstatus+=1
            if(self.iterstatus == self.len):
                self.iterstatus = 0
                self.finished_at_least_once = True
            return dict(zip(self.iterkeys, cc))

    def __str__(self):
        s = OUTPUT_CELLTYPE_BEGIN + self.name + '\n'
        c = iter(self)
        while(not self.finished_at_least_once):
           dicc = next(c)
           for k, v in zip(dicc.keys(), dicc.values()):
             s += k + CELLTYPE_SEP + str(v) + '\n'
           s += OUTPUT_CELLTYPE_END + '\n'
        self.finished_at_least_once = False    
        return s





class test:
    def __init__(self, a):
        self.value = [i for i in range(a)]
        self.i = 0
        self.iterating = False
    def __iter__(self):
        self.iterating = True
        while (self.iterating):
            yield self.value[self.i]
            self.i = (self.i + 1)%len(self.value)



def main():
    p = paramContainer(sys.argv[1], sys.argv[2])
    p.produceOutputs()


if(__name__ == "__main__"):
    main()
