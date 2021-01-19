import make_hexagonal_grid as hx
import argparse


class argArray:
    """
    Takes arguments identical to make_hexagonal_grid, but everything as STRINGS, each with different param values separated
    by '|'. Generates all combinations of parameters to make a new hexagonal grid.
    """
    iternames = ['Size', 'Rows', 'Cols', 'Noise', 'Pull', 'Strecht', 'Rotate', 'StaticVertices', 'Springs', 'SpringLength', 'Hinge', 'Veins']

    def __init__(self, args):
        self.argumentValues = self.__getArgDictFromInput(args)
        print("ALL: ")
        print(self.argumentValues)
        self.printed = 0
    def getDictFromList(self, args, outname, read):
        argdict = {}
        argdict.setdefault('Size', float(args[0]))
        argdict.setdefault('Rows', int(args[1]))
        argdict.setdefault('Cols', int(args[2]))
        argdict.setdefault('Noise', float(args[3]))
        argdict.setdefault('Pull', float(args[4]))
        argdict.setdefault('Strecht', float(args[5]))
        argdict.setdefault('Rotate', float(args[6]))

        argdict.setdefault('StaticVertices', args[7])
        argdict.setdefault('Springs', args[8])
        argdict.setdefault('SpringLength', float(args[9]))

        argdict.setdefault('Hinge', int(args[10]))
        argdict.setdefault('Veins', args[11])
        argdict.setdefault('Outname', outname)
        argdict.setdefault('Read', read)
        return argdict
    def getArgDictionaries(self):
        from itertools import product
        vallists = [self.argumentValues.get(k) for k in argArray.iternames]
        print("_________\n", [i for i in vallists], "\n_________\n")
        vallists = [self.argumentValues.get(k) for k in argArray.iternames]
        for combo in product(*vallists):
            yield self.getDictFromList(combo, '_'.join([self.argumentValues["Outname"], str(self.printed)]), self.argumentValues["Read"])
            self.printed+=1
    def __getArgDictFromInput(self, args):
        argdict = {}
        argdict.setdefault('Size', args.Size.split("|"))
        argdict.setdefault('Rows', args.Rows.split("|"))
        argdict.setdefault('Cols', args.Cols.split("|"))
        argdict.setdefault('Noise', args.Noise.split("|"))
        argdict.setdefault('Pull', args.Pull.split("|"))
        argdict.setdefault('Strecht', args.Strecht.split("|"))
        argdict.setdefault('Rotate', args.Rotate.split("|"))

        argdict.setdefault('StaticVertices', args.StaticVertices.split("|"))
        argdict.setdefault('Springs', args.Springs.split("|"))
        argdict.setdefault('SpringLength', args.SpringLength.split("|"))

        argdict.setdefault('Hinge', args.Hinge.split("|"))
        argdict.setdefault('Veins', args.Veins.split("|"))
        argdict.setdefault('Outname', args.Outname + '_s'+str(args.Size) + '_' + str(args.Rows) + 'x' + str(args.Cols) + '_n' + str(args.Noise))
        argdict.setdefault('Read', args.Read)
        print('size: ', args.Size, '; num. rows: ', args.Rows, '; num. cols: ', args.Cols, '; noise: ', args.Noise, '; strecht: ', args.Pull, '; output files: ', argdict['Outname'])
        return argdict


parser = argparse.ArgumentParser(description='Hexagonal grid arguments.')
parser.add_argument('-o', '--Outname', metavar='outname', type=str, default = "hexgrid", 
                                        help='Identifier. Used as prefix of all output files. ')
parser.add_argument('-i', '--Read', metavar='Read', type=str, default = "", 
                                        help='Identifier. Read hexagonal grid. ')
parser.add_argument('-s', '--Size', metavar='size', type=str, default = str(hx.size),
                                        help='Separate values with "|". Hexagon size: distance between center and vertices. ')
parser.add_argument('-r', '--Rows', metavar='rows', type=str, default = str(hx.nrow),
                                        help='Number rows. ')
parser.add_argument('-c', '--Cols', metavar='cols', type=str, default = str(hx.ncol),
                                        help='Number of columns.')
parser.add_argument('-n', '--Noise', metavar='noise', type=str, default = str(hx.noise),
                                        help='Maximum proportion of size that each vertex is moved randomly')
parser.add_argument('-t', '--StaticVertices', metavar='static_vertices', type=str, default = ";;;", 
                                        help='Vertices in margin that are static. Sytax: "left;top;bottom;right", where each position can be a series of numbers or ranges (eg. 0:5) separated by commas. Eg: "0,4:6;;;0:10,12,15"')
parser.add_argument('-p', '--Springs', metavar='springs', type=str, default = "", 
                                        help='Vertices in margin that can move but are attached to springs. Sytax: "left;top;bottom;right", where each position can be a series of numbers or ranges (eg. 0:5) separated by commas. Eg: "0,4:6;;;0:10,12,15"')
parser.add_argument('-l', '--SpringLength', metavar='spring_length', type=str, default = str(3.5),
                                        help='Length of springs')
parser.add_argument('-g', '--Hinge', metavar='Hinge_pos', type=str, default = str(-1),
                                        help='Hexagons in columns in range [0:g] will be considered hinge (value 1).')
parser.add_argument('-v', '--Veins', metavar='Veins_pos', type=str, default = '',
                                        help='Hexagons in rows specified (as integers) will be considered veins ( value 2).')
parser.add_argument('-f', '--Pull', metavar='Change relative size of angles', type=str, default = str(0),
                                        help='Value of side angle (default is 60 for a regular hexagon')
parser.add_argument('-k', '--Strecht', metavar='Strecht', type=str, default = str(0),
                                        help='Scale horizontally by a factor of k')
parser.add_argument('-j', '--Rotate', metavar='Rotate', type=str, default = str(0),
                                        help='Rotate the whole grid specified angles')

def main():
        args = parser.parse_args()
        argarray = argArray(args)
        for argdict in argarray.getArgDictionaries():
            print("NUM: ", argarray.printed)
            print(argdict)
            print("\n***\n")
            hex = hx.HexGrid(**argdict)
            hex.writeGrid()
            #hex.plotHex()

if __name__ == '__main__':
        main()
