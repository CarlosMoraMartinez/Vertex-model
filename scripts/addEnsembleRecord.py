import os
import sys
import argparse
import json
from datetime import date

import numpy as np
import pandas as pd


pardir = "./param_files/"
parfile_end = "allConds.csv"
pic_format = '.png'
picsdir_end = "_all_final/"

outdir = "SummaryResults/"
no_plot_cols = ["Unnamed", "name", "Result"]
null_result_value = "-"

def dropIdenticalColumns(data):
    vary = [a for a in data.columns if np.unique(data[a]).shape[0] != 1]
    return data[vary]

def getConditionsThatChange():
    parfile = [pardir + i for i in os.listdir(pardir) if i.endswith(parfile_end)].pop()
    conds=pd.read_csv(parfile, sep="\t")
    c2 = dropIdenticalColumns(conds)
    cdict = {c:list(set(c2[c])) for c in c2.columns if not any([i in c for i in no_plot_cols])}
    return c2, cdict

def getGoodSims(c2):
    if(not "Result" in c2.columns):
        return []    
    return c2.iloc[[i for i in range(c2.shape[0]) if "ok" in c2.loc[i]["Result"] ]]["name"].tolist()

def guessWingName():
    wings = []
    for f in os.listdir():
        if(f.endswith(picsdir_end.replace("/", ""))):
            wings.append( f.replace(picsdir_end.replace("/", ""), ""))
    return wings

def getValuesByResultType(c2, chgconds):
    if(not "Result" in c2.columns):
        print("No Result column", c2.columns)
        print(c2)
        return []    
    res = {}
    vnames = [v for v in c2.columns if not any([i in v for i in no_plot_cols])]
    restypes = set([j for i in c2["Result"] for j in i.split(";")])
    for resval in restypes:
        if(resval == null_result_value):
            continue
        indices = [i for i in range(c2.shape[0]) if resval in c2.iloc[i]["Result"].split(";")]
        c3 = c2.iloc[indices]
        print("\n\n\nResult type: %s"%(resval))
        print(c3)
        resdict = {"total_with_this_result":c3.shape[0]}
        for v in vnames:
            set2 = set(c2[v])
            set3 = set(c3[v])
            if(len(set2) == len(set3)): 
                continue
            present = [i for i in set3]#.tolist() #The .item() is necessary because numpy numbers cannot be converted into json (strings can)
            absent = [i for i in set2 if i not in present]
            count = {str(i):len([0 for j in c3[v] if j == i]) for i in present}
            resdict.setdefault(v.split(" ")[0], {"present_values":present, "absent_values":absent, "count":count})
        try:
            #This is only useful to debug, since some data types are not directly convertible to json
            print("Trying conversion to json...")
            aux = json.dumps(resdict)
            print("ok")
        except:
            print("Failed JSON conversion of result %s"%(resval))
            print("\n****\n", resdict, "\n****\n")
        res.setdefault(resval, resdict) 
    return res   
            

def writeAsJson(res):
    with open(outdir + res['ensemble_name'] + '.json', 'w') as fp:
        json.dump(res, fp, indent=4)



parser = argparse.ArgumentParser(description='Make summary of ensemble. In order to work properly, a tab separated file with the parameters of each simulation must be present in ./ensemble/param_files. The script must be called in the folder in which the ensemble directory (./ensemble) is. Optionally, the conditions file can have an extra column called Result. In this column, standardized descriptions of results can be added. Value ok in the Result column is understood as success. Different values separated by ; can be used to compute summary separated for each of them')
parser.add_argument('-e', '--ensemble', metavar='ensemble', type=str, default = '', 
                    help='Name of ensemble')
parser.add_argument('-p', '--purpose', metavar='purpose', type=str, default = '', 
                    help='Purpose or idea behind this ensemble')
parser.add_argument('-c', '--code_changed', metavar='code_changed', type=str, default = '', 
                    help='Changes in code incorporated to this ensemble')
parser.add_argument('-f', '--from_condition', metavar='from_condition', type=str, default = '', 
                    help='Condition or ensemble from which this ensemble is derived')
parser.add_argument('-r', '--result', metavar='result', type=str, default = '', 
                    help='Some string indicating results')
parser.add_argument('-w', '--wings', metavar='wings', type=str, default = '', 
                    help='Wing names in ensemble separated by ,. If absent, tries to guess (a directory with _all_final ending name must exist)')

def main():
    args = parser.parse_args()
    ensemble = args.ensemble
    try:
        os.chdir(ensemble)
    except:
        print("Ensemble %s not found"%(ensemble))
        return
    wings = args.wings.split(",") if args.wings else guessWingName()
    condtab, chgconds = getConditionsThatChange()
    goodsims = getGoodSims(condtab)
    pbrt = getValuesByResultType(condtab, chgconds)    
    res  = {
        "ensemble_name":ensemble, 
        "size":str(condtab.shape[0]),
        "date":str(date.today()),
        "purpose":args.purpose,
        "code_changed":args.code_changed,
        "from_condition":args.from_condition,
        "wings":wings,
        "result":args.result,
        "conditions_varied":chgconds,
        "good_simulations":goodsims,
        "params_by_result_type":pbrt
    }

    os.chdir("../")
    writeAsJson(res)

if __name__ == '__main__':
    main()

