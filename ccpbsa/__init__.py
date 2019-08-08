from .CCPBSA2 import *
import argparse
import os
import sys


#def gxg(
#    flags,
#    min_mdp,
#    energy_mdp,
#    mdrun_table,
#    pbeparams
#):
#    gxg = GXG(flags, min_mdp, energy_mdp, mdrun_table, pbeparams)
#    gxg()
#    gxg.G.to_csv('gxg.csv')
#    
#    return gxg.G


def main():
    pkgpath = __path__[0]
    cliparser = argparse.ArgumentParser(prog='ccpbsa')
    options = cliparser.add_argument_group("OPTIONS")

    cliparser.add_argument(
        "routine",
        help="Called during installation to setup default parameters.",
        choices={'setup', 'stability'}
    )
    options.add_argument(
        "-w", "--wildtype",
        help=".pdb file of the wildtype protein"
    )
    options.add_argument(
        "-m", "--mutations",
        help="a .txt file with the list of mutations"
    )
    options.add_argument(
        "-f", "--flags",
        help="the flags, which should be passed to CONCOORD and GROMACS. At \
        least the number of structures in CONCOORD and the forcefield and water \
        model should be specified in there",
        default=pkgpath+'/parameters/flags.txt'
    )
    options.add_argument(
        "--chains",
        help="Name of chains in the .pdb file.",
        default='A',
        nargs='+'
    )
    options.add_argument(
        "--stability-parameters",
        help="scaling factors of the ddG stability calculations.",
        default=pkgpath+'/parameters/fit_stability.txt',
    )
    options.add_argument(
        "--affinity-parameters",
        help="scaling factors of the ddG affinity calculations.",
        default=pkgpath+'/parameters/fit_affinity.txt',
    )
    options.add_argument(
        "--energy-mdp",
        default=pkgpath+'/parameters/energy.mdp',
        help=".mdp file for GROMACS energy evaluations. Default: epsilon-r = 1"
    )
    options.add_argument(
        "--gxg-table",
        default=pkgpath+'/parameters/GXG.csv',
        help="GXG table used for stability change calculations"
    )
    options.add_argument(
        "-v",
        help="Print stdout and stderr of the programs",
        action='store_true'
    )

    cliargs = cliparser.parse_args()

    if cliargs.routine == 'setup':
        
        with open(pkgpath+'/parameters/flags.txt', 'a') as paramfile:
            paramfile.write("-tablep="+pkgpath+'/parameters/table4r-6-12.xvg\n')
            paramfile.write("-table="+pkgpath+'/parameters/table4r-6-12.xvg\n')
            paramfile.write("[grompp]\n")
            paramfile.write("-maxwarn=1\n")
            paramfile.write("-f="+pkgpath+'/parameters/min.mdp\n')
            paramfile.write("[gropbe]\n")
            paramfile.write(pkgpath+'/parameters/gropbe.txt\n')
            sys.exit(0)
#    if cliargs.gxg:
#        print("Making a new GXG table.")
#        GXG_ = gxg(
#            cliargs.flags,
#            cliargs.minimization_mdp,
#            cliargs.energy_mdp,
#            cliargs.dielectric_table,
#            cliargs.gropbe
#        )
#        print(GXG_)
#        sys.exit()
    if cliargs.v:
        verbose = 1
    else:
        verbose = 0

    if cliargs.routine == 'stability':
        data = DataGenerator(
            wtpdb = cliargs.wildtype,
            mutlist = cliargs.mutations,
            flags = cliargs.flags,
            calculate = 'stability',
            chains = cliargs.chains,
            spmdp = cliargs.energy_mdp,
            verbosity = verbose
        )
        data.fullrun()
        search = DataCollector(data)
        G = search.search_data()
        print("G values:")
        print(search.G)
        G.to_csv("G.csv")

        with open(cliargs.stability_parameters, 'r') as fit:
            parameters = fit.readlines()
            parameters = [l[:-1] for l in parameters] # Remove newlines
            parameters = [l.split("=") for l in parameters]
            parameters = dict([(l[0], float(l[1])) for l in parameters])

        search.dstability(cliargs.gxg_table)
        print("dG folded values:")
        print(search.dG)
        print("dG unfolded values (GXG):")
        print(search.dG_unfld)
        search.dG.to_csv("dG_fold.csv")
        search.dG_unfld.to_csv("dG_unfold.csv")

        ddG = search.ddstability()
        print("ddG values:")
        print("without fit:")
        print(search.ddG)
        search.ddG.to_csv("ddG.csv")
        ddG_fit = search.fitstability(**parameters)
        print("with fit:")
        print(ddG_fit)
        search.ddG.to_csv("ddG_fit.csv")

    elif cliargs.mode == 'affinity':
        with open(cliargs.affinity_parameters, 'r') as fit:
            parameters = fit.readlines()
            parameters = [l[:-1] for l in parameters] # Remove newlines
            parameters = [l.split("=") for l in parameters]
            parameters = dict([(l[0], float(l[1])) for l in parameters])

        search.daffinity()

        search.ddaffinity(
            parameters['alpha'],
            parameters['beta'],
            parameters['gamma'],
            parameters['c'],
            parameters['pka']
        )

    else:
        print("BAD MODE INPUT")
