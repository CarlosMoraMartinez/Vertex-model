import os
import sys
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import imageio


pardir = "./param_files/"
parfile_end = "allConds.csv"
pic_format = '.png'
picsdir_end = "_all_final/"

outfolder_base = "./combined_plots"
no_plot_cols = ["Unnamed", "name", "Result"]


def dropIdenticalColumns(data):
    vary = [a for a in data.columns if np.unique(data[a]).shape[0] != 1]
    return data[vary]

def getImages(picsdir, wing):
    img_files= [picsdir + i for i in os.listdir(picsdir) if i.endswith(pic_format)]
    images = { i.split("/")[1].replace("_" + wing + pic_format, ""):imageio.imread(i) for i in img_files}
    imgsize = [i for i in images.values()][0].shape
    absence = np.zeros(imgsize)
    return (images, absence)

def getParams():
    parfile = [pardir + i for i in os.listdir(pardir) if i.endswith(parfile_end)].pop()
    conds=pd.read_csv(parfile, sep="\t")
    c2 = dropIdenticalColumns(conds)
    return c2

def readData(picsdir, wing):
    files= [picsdir + i.replace(".points", "") for i in os.listdir(picsdir) if i.endswith(".points")]   
    data = {i.split("/")[1].replace("_" + wing, ""):readWing(i) for i in files}
    return (data, None)

def readWing(basename):
    e=pd.read_csv(basename + '.edges', sep="\t") #file with edges, pasted from .out file (need to prepare this first)
    pf=open(basename + '.points', "r")
    p=int(pf.readline()) #discard first line
    p=[[float(j) if '.' in j else int(j) for j in i.split("\t")] for i in pf.read().split("\n") if i != ""]
    p=pd.DataFrame(p,columns=["x", "y", "ind", "movable"])
    return (p, e)

#p=pd.read_csv(basename + '.points', sep="\t") #file with points (like .points but with header (x	y	ind	movable); need to prepare this first)
def plotImage(ax, wing, limits):
    ax.set_xticks([]);
    ax.set_yticks([]);
    ax.set_xticklabels([]);
    ax.set_yticklabels([])
    ax.imshow(wing)

def plotBorders(ax, wing, limits):
    ax.set_xticks([]);
    ax.set_yticks([]);
    ax.set_xticklabels([]);
    ax.set_yticklabels([])
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal")
    if(wing is None):
        ax.plot([limits[0], limits[1]], [limits[2], limits[3]], color="black")
        ax.plot([limits[0], limits[1]], [limits[3], limits[2]], color="black")
        return
    p, e = wing
    colors = ['black', 'green', 'yellow', 'orange', 'blue', 'purple','red', 'peru']
    for i in range(e.shape[0]):
        ind = e["vertices"][i].split(',')
        paux = p[(p['ind'] == int(ind[0])) | (p['ind'] == int(ind[1]))]
        cc = colors[e["type"][i]] if(e["type"][i] < len(colors)) else colors[-1]       
        ax.plot(paux['x'], paux['y'], color = cc)


def plotBordersTension(ax, wing, limits):
    ax.set_xticks([]);
    ax.set_yticks([]);
    ax.set_xticklabels([]);
    ax.set_yticklabels([])
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal")
    if(wing is None):
        ax.plot([limits[0], limits[1]], [limits[2], limits[3]], color="black")
        ax.plot([limits[0], limits[1]], [limits[3], limits[2]], color="black")
        return
    p, e = wing
    mint = np.min(e['tension'][e["type"] != 4])
    maxt = np.max(e['tension'][e["type"] != 4]) #borders always have bigger tension
    e['tension'][e["type"] == 4] = maxt
    normt = (e['tension'] - mint)/(maxt - mint) + 0.1
    normt = normt/np.max(normt)
    colors = ['black', 'green', 'yellow', 'orange', 'blue', 'purple','red', 'peru']
    for i in range(e.shape[0]):
        ind = e["vertices"][i].split(',')
        paux = p[(p['ind'] == int(ind[0])) | (p['ind'] == int(ind[1]))]
        cc = colors[e["type"][i]] if(e["type"][i] < len(colors)) else colors[-1]       
        ax.plot(paux['x'], paux['y'], color = 'black', alpha = normt[i] )


#c2: param table with only variables that change
#c1 and c2: variables that we are plotting
#diff_conds: all conditions to plot this time
#absence: thing to plot instead

def multiplot(c2, diff_conds, ca, cb, current_cond, images, absence, plotFun=plotImage, limits=[-10,80,-10,55], plot=False, outfolder=outfolder_base):
    c_thiscond=c2.iloc[ [np.all([c2[j].iloc[i] == diff_conds[j].iloc[current_cond] for j in diff_conds.columns]) for i in range(c2.shape[0])] ]
    nrows = np.unique(c2[ca]).shape[0]
    ncols = np.unique(c2[cb]).shape[0]
    fig, ax =plt.subplots(nrows=nrows, ncols=ncols)
    fig.text(0.5, 0.04, cb.split(" ")[0], ha='center', va='center')
    fig.text(0.06, 0.5, ca.split(" ")[0], ha='center', va='center', rotation='vertical')
    for xpos, va in enumerate(np.unique(c2[ca])):
        for ypos, vb in enumerate(np.unique(c2[cb])):
            iname = list(c_thiscond["name"].iloc[np.array(c_thiscond[ca] == va) & np.array(c_thiscond[cb] == vb)])[0]
            plotFun(ax[xpos][ypos], images.get(iname, absence), limits)
            if(xpos == nrows - 1): 
                ax[xpos][ypos].set_xlabel(str(vb))
            if(ypos == 0):
                ax[xpos][ypos].set_ylabel(str(va))
            ax[xpos][ypos].set_title(iname)
            #ax[xpos][ypos].set_title( str(va) + ", " + str(vb))
            #ax[xpos][ypos].set_title(ca.split(" ")[0] + "=" + str(va) + ", " + cb.split(" ")[0] + "=" + str(vb))
    if(plot):
        plt.show()
    else:
        out_cond = outfolder + "/" +  ca.split(" ")[0] + ":" + cb.split(" ")[0] + ":cond" + str(current_cond) + '.svg'
        plt.savefig(out_cond, format='svg', dpi=1200)
    plt.close()

def plotAll(c2, images, absence, plotFun=plotImage,  limits=[-10,80,-10,55], outfolder=outfolder_base, plotpars=[]):
    for vnum1 in range(c2.columns.shape[0]):
        ca = c2.columns[vnum1]
        if(any([i in ca for i in no_plot_cols])):
            continue
        if(plotpars and not ca in plotpars):
            continue
        for vnum2 in range(vnum1 + 1, c2.columns.shape[0]):
            cb = c2.columns[vnum2]
            if(any([i in cb for i in no_plot_cols])):
                continue
            itcols = list(map(lambda x: x != ca and x != cb and not any([i in x for i in no_plot_cols]), c2.columns)) #which parameters must be identical
            diff_conds = c2[c2.columns[itcols]].drop_duplicates() #Unique conditions after removing parameters ca and cb (for each row will perform a different plot)
            for current_cond in range(diff_conds.shape[0]):
                multiplot(c2, diff_conds, ca, cb, current_cond,  images, absence, plotFun, limits, False, outfolder)
            diff_conds["multiplot_cond"] = list(range(diff_conds.shape[0]))
            out = outfolder + "/" +  ca.split(" ")[0] + ":" + cb.split(" ")[0] + '.csv'
            diff_conds.to_csv(out, sep="\t")
            print("Printed: " + out)

def guessWingName():
    wings = []
    for f in os.listdir():
        if(f.endswith(picsdir_end.replace("/", ""))):
            wings.append( f.replace(picsdir_end.replace("/", ""), ""))
    return wings

def make_grids_wing(wing, mode, plotTension, limits, plotpars):
    outfolder = '_'.join([outfolder_base, wing])
    try:
        os.mkdir(outfolder)
    except:
        print("folder %s already exists"%(outfolder))
    picsdir = wing + picsdir_end

    print("wings: ", wing)
    print("limits: ", limits)
    print("plot from images: ", int(mode))
    print("plot tension: ", int(plotTension))
    c2 = getParams()
    if(mode):
        images, absence = getImages(picsdir, wing)
        plotFun = plotImage
        plotTension = False
    else:
        images, absence = readData(picsdir, wing)
        if(plotTension):
            plotFun = plotBordersTension
        else:
            plotFun = plotBorders
    plotAll(c2, images, absence, plotFun, limits, outfolder, plotpars) 

parser = argparse.ArgumentParser(description='Plot grid arguments.')
parser.add_argument('-w', '--wingName', metavar='wing', type=str, default = '', 
                    help='Name of wing (for instance wing2E, budsmall, etc). Can provide more than one separated by comma (wing1,wing2). If not provided tries to guess.')
parser.add_argument('-t', '--plotTension', metavar='plot_tension', type=bool, default = False, 
                    help='Plot edge intensity according to tension (overriden by plotFromImages)')
parser.add_argument('-l', '--fixedLimits', metavar='fixed_limits', type=str, default = '-1,100,-1,260', 
                    help='Use these max and min coordinates in all plots. Use this format: x0,xmax,y0,ymax')
parser.add_argument('-i', '--plotFromImages', metavar='plot_imgs', type=bool, default = False, 
                    help='If True, uses .png images to make composites (faster). If False, reads .points and .edges and makes new plots')
parser.add_argument('-p', '--PlotParams', metavar='plot_params', type=str, default = '', 
                    help='If present, only plots combinations with these params. Use exact names separated by ,')
def main():
    args = parser.parse_args()

    mode = args.plotFromImages
    plotTension = args.plotTension
    wings = args.wingName.split(',') if args.wingName != "" else guessWingName()
    limits = [float(i) for i in args.fixedLimits.split(",")]
    plotpars = args.PlotParams
    plotpars = plotpars.split(',') if plotpars else ""

    for wing in wings:
        make_grids_wing(wing, mode, plotTension, limits, plotpars)

if __name__ == '__main__':
    main()


