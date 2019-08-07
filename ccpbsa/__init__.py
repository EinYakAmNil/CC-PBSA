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
        "--mode",
        help="Calculate either the change of stability or affinity changing \
        effects of mutations.",
        default='stability',
        choices={'stability', 'affinity'}
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
        "--minimization-mdp",
        default=pkgpath+'/parameters/min.mdp',
        help=".mdp file for GROMACS energy minimization. Default: Uses lookup \
        table in \"--dielectric-table\""
    )
    options.add_argument(
        "--energy-mdp",
        default=pkgpath+'/parameters/energy.mdp',
        help=".mdp file for GROMACS energy evaluations. Default: epsilon-r = 1"
    )
    options.add_argument(
        "--dielectric-table",
        default=pkgpath+'/parameters/table4r-6-12.xvg',
        help="Look up table for GROMACS to simulate continuum solvent. \
        Default: epsilon(r) = 4r"
    )
    options.add_argument(
        "--gropbe",
        default=pkgpath+'/parameters/gropbe.txt',
        help="Parameters for gropbe. Defaults are the same as in the paper"
    )
    options.add_argument(
        "--gxg-table",
        default=pkgpath+'/parameters/GXG.csv',
        help="GXG table used for stability change calculations"
    )
    options.add_argument(
        "--gxg",
        help="Creates a new gxg.csv for dG calculation of the unfolded state.",
        action='store_true'
    )
    options.add_argument(
        "-v",
        help="Print stdout and stderr of the programs",
        action='store_true'
    )

    cliargs = cliparser.parse_args()

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

    data = DataGenerator(
        wtpdb = cliargs.wildtype,
        mutlist = cliargs.mutations,
        flags = cliargs.flags,
        calculate = cliargs.mode,
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

#    Very simple parser for reading in alpha, beta...
    if cliargs.mode == 'stability':
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
        search.dG.to_csv.to_csv("dG_fold.csv")
        search.dG_unfld.to_csv("dG_unfold.csv")

        ddG = search.ddstability()
        print("ddG values:")
        print(search.ddG)
        search.ddG.to_csv("ddG.csv")
        ddG_fit = search.fitstability(**parameters)
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
