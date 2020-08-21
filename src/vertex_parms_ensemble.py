import sys
import os
import argparse
import json

import numpy as np
import pandas as pd


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
INPUT_EXTENSION = '.vp'
OUTPUT_EXTENSION = '.vp'

def getNum(string):
    if('.' in string):
        return float(string)
    else: 
        return int(string)

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
        fnames = []
        for i in range(num):
            parset = dict(zip(allnames, combos[i]))
            self.writeFile(parset, self.outname + '_' + str(i))
            fnames.append(self.outname + '_' + str(i))
        self.produceDataFrame(allnames, combos, fnames)

    def produceDataFrame(self, allnames, combos, fnames):
        df = pd.DataFrame(data = dict(zip([s[0:s.find('#')].rstrip() for s in allnames], zip(*combos) )))
        df['name'] = fnames
        df.set_index('name')
        df.to_csv(self.outname + '/' + self.outname + '_allConds.csv', sep = '\t')

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

class paramRandomContainer(paramContainer):
    def __init__(self, fname, outname, ran):
        self.fname = fname
        self.outname = outname
        self.ran = ran
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
                parobj = randomCellPar(par, self.ran)
                i+=1       
                self.params.append(parobj)      
            elif(rawstrings[i].startswith(PAR_BEGIN)):
                par.append(rawstrings[i])
                par.append(rawstrings[i + 1])
                parobj = randomPar(par, self.ran)
                i+=2
                self.params.append(parobj)    
            else:
                print("Error: line not recognized")
                i+=1
    def produceOutputs(self):
        from itertools import product
        allnames = [p.name for p in self.params]
        allparsets = []
        fnames = []
        for i in range(self.ran):
            vals = []
            for p in self.params:
                if(type(p.value) is dict):
                    d = dict()
                    for k, v in zip(p.value.keys(), p.value.values()):
                        if(p.atomic[k]):
                            print(p.name, ", cell ", k, ", atomic")
                            d.setdefault(k, v)
                        else:
                            print(p.name, ", cell ", k, ", not atomic")
                            d.setdefault(k, v[i])
                    vals.append(d)
                else:
                    if(p.atomic):
                        print(p.name, ", atomic")
                        vals.append(p.value)
                    else:
                        print(p.name, ", not atomic")
                        vals.append(p.value[i])
            parset = dict(zip(allnames, vals))
            print("\n\n")
            self.writeFile(parset, self.outname + '_' + str(i))
            allparsets.append(vals)
            fnames.append(self.outname + '_' + str(i))
        self.produceDataFrame(allnames, allparsets, fnames)

class randomPar:
    def __init__(self, parlist, ran):
        self.name = parlist[0].strip(PAR_BEGIN)
        self.ran = ran
        self.atomic, self.value = self.getvalue(parlist[1])
    def getvalue(self, parlist):
        if(RANDOM_VALUE in self.name):
            vv = parlist.split(SEQ_SEP)
            s = getNum(vv[0])
            e = getNum(vv[1])
            if(type(s) is float):
                value = np.random.uniform(s, e, self.ran)
            else:
                value = np.random.random_integers(s, e, self.ran)
            if(self.ran > 1):
                atomic = False
            else:
                atomic = True
        else:
            value = getNum(parlist)
            atomic = True
        return (atomic, value)
    def __str__(self):
        s = OUTPUT_PAR_BEGIN + self.name + '\n'
        s += str(self.value)
        return s

class randomCellPar(randomPar):
    def __init__(self, parlist, ran):
        self.name = parlist.pop(0).strip(PAR_BEGIN)
        self.ran = ran
        self.value = dict()
        self.atomic = dict()
        if(CELL_TYPE_PRODUCT in self.name):
            for p in parlist:
                p = p.split(CELLTYPE_SEP)
                if(SEQ_SEP  in p[1]):
                    at, val = self.getvalue(p[1])
                else:
                    at = True
                    val = getNum(p[1])
                self.atomic.setdefault(p[0], at)
                self.value.setdefault(p[0], val)
        else:
            at0 = None
            val0 = None
            for p in parlist:
                p = p.split(CELLTYPE_SEP)
                if(SEQ_SEP  in p[1] and val0 is None):
                    at0, val0 = self.getvalue(p[1])
                    at, val = (at0, val0)
                elif(SEQ_SEP in p[1]):
                    at, val = (at0, val0)
                else:
                    at = True
                    val = getNum(p[1])
                self.atomic.setdefault(p[0], at)
                self.value.setdefault(p[0], val)
    def __str__(self):
        s = OUTPUT_CELLTYPE_BEGIN + self.name + '\n'
        for k, v in zip(self.value.keys(), self.value.values()):
            s += k + ':' + str(v) + '\n'
        return s

class param:
    def __init__(self, parlist):
        self.name = parlist[0].strip(PAR_BEGIN)
        print(self.name)
        self.atomic, self.value = self.processParam(parlist[1])
        self.iterstatus = 0
        self.finished_at_least_once = False

    def processParam(self, string):
        elements1 = string.split(DISJOINT_SEP)
        if(len(elements1) == 1 and not SEQ_SEP in elements1[0]):
            atomic = True
            value = getNum(elements1[0])
        else:
            atomic = False
            value = []
            for i in elements1:
                if(SEQ_SEP in i):
                    i = i.split(SEQ_SEP)
                    for j in np.linspace(getNum(i[0]), getNum(i[1]), getNum(i[2])):
                        if('.' not in i[0]):
                            j = int(j)
                        value.append(j)
                else:
                    value.append(getNum(i))
        return(atomic, value)
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
        print(self.name)
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


class varyOneByOne:
    """
    Instead of reading an ensemble file with .vp format and performing all-by-all combinations, this reads a .json file which indicates which parameters to vary. 
    In the .json file, in the [param][values] field, two values are indicated: proportion of total variation (a) and proportion of variation between simulations (b).
    Therefore, all parameter files with values X - a*x and X + a*X will be constructed in intervals of b*X.
    Instead of taking the all-by-all approach, only one parameter will be varied each time, and the rest will remain as in the default files (which must be specified 
    in the .json file).
    This class can act on several param files, if specified in the .json file. 
    """
    def __init__(self, iterdef):
        import json
        self.outname = iterdef
        self.tovary = json.load(open(iterdef + '.json', 'r'))
        print(list(self.tovary.keys()))
        self.paramfiles = {i:self.readParamFile(i) for i in list(self.tovary.keys())}
        print("\nParam files read")
    def readParamFile(self, fname):
        return paramContainer(fname + INPUT_EXTENSION, '')
    def produceOutputs(self):
        for part, params in self.tovary.items(): #iter over different files
            for pname, p in params.items(): # Iter over parameters that must be changed in each file
                param_num = [i for i, prm in enumerate(self.paramfiles[part].params) if prm.name[0:prm.name.find('#')].rstrip() == pname]
                if(not param_num):
                    print("ERROR: Parameter not found. In Vary file: %s"%(pname))
                    continue
                param_num  = param_num.pop()
                default_value = self.paramfiles[part].params[param_num].value
                print(pname, default_value)
                newvalues = self.GetParamValues(default_value, p)
                for k, nwv in enumerate(newvalues): #Iter over values of parameter pname
                    self.paramfiles[part].params[param_num].value = nwv
                    fname = '-'.join([part, pname, str(k)])
                    self.writeFile(self.paramfiles[part], fname)
                self.paramfiles[part].params[param_num].value = default_value
    def GetParamValues(self, default, mode):
        if(mode["type"] == "regular"):
            values = np.linspace(default - abs(default)*mode["values"][0], default + abs(default)*mode["values"][0], int(2 / mode["values"][1] + 1))
            if(mode["absol"] == "yes"):
                values = values[values >= 0]
        else:
            this_val = []
            for k,v in default.items():
                v = v[0]
                this_val.append(np.linspace(v - abs(v)*mode["values"][k][0], v + abs(v)*mode["values"][k][0], int(2 / mode["values"][k][1] + 1)))
            if(mode["together"] == "all"):
                values = [{k:l for k, l in zip(default.keys(), vals)} for vals in zip(*this_val)]
            else:   
                values = []       
                for k,vals in zip(default.keys(), this_val):
                    for v in vals:
                        values.append({k2:v2[0] if k2 != k else v for k2, v2 in default.items()})          
        return values      
    def writeFile(self, pcontainer, fname):
        s = '# ' + fname + '\n'       
        for p in pcontainer.params:
            v = p.value
            if(type(v) is dict):
                s += OUTPUT_CELLTYPE_BEGIN + p.name + '\n' 
                for k2, v2 in zip(v.keys(), v.values()): 
                    if(type(v2) is list and v2):
                        v2 = v2[0] 
                    s += k2 + CELLTYPE_SEP + str(v2) + '\n'            
                s += OUTPUT_CELLTYPE_END + '\n\n' 
            else:
                s += OUTPUT_PAR_BEGIN + p.name + '\n' + str(v) + '\n'
        if not os.path.exists(self.outname):
            os.makedirs(self.outname)
        outfile = self.outname + '/' + fname
        f = open(outfile + OUTPUT_EXTENSION,"w+")
        f.write(s)
        f.close()       


parser = argparse.ArgumentParser(description='Generate folder with ensemble of parameter files for VertexSystem. The input file is the standard .vp, except that ranges can be specified with start;end;size, and specific values can be specified with n0,n1,n2, etc. Additionally, for CellType parameters, the symbol % in the name indicates that all combinations of values will be used for this parameter (if ommitted, vectors of values for each cell type should be of equal size or size 1).')
parser.add_argument('-o', '--Outname', metavar='outname', type=str, default = "hexgrid", 
                                        help='Identifier. Used as prefix of all output files. ')
parser.add_argument('-i', '--Inputname', metavar='inputname', type=str, default = "", 
                                        help='Identifier. Used as prefix to read files. ')
parser.add_argument('-r', '--Random', metavar='randompars', type=int, default = 0, 
                                        help='Random set of parameters (1) or specified combinations (0). 0 is default. Only parameters with a ? symbol in their name will be randomly set.')
parser.add_argument('-s', '--Serial', metavar='serialpars', type=str, default = '', 
                                        help='.json file (without format) specifying which param files and parameters vary serially. Overwrites everything else')

def main():
    args = parser.parse_args()
    iname = args.Inputname
    oname = args.Outname
    ran = args.Random
    serial = args.Serial
    if(serial):
        p = varyOneByOne(serial)
    elif(ran == 0):
        p = paramContainer(iname, oname)
    else:
        p = paramRandomContainer(iname, oname, ran)
    p.produceOutputs()


if(__name__ == "__main__"):
    main()
