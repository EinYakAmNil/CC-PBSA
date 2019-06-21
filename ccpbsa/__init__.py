from .CCPBSA import *
import argparse
import os
import sys


def gxg(
    flags,
    min_mdp,
    energy_mdp,
    mdrun_table,
    pbeparams
):
    gxg = GXG(flags)
    gxg()
    gxg.ener_df.to_csv('gxg.csv')
    
    return gxg.ener_df


def main():
    pkgpath = __path__[0]
    cliparser = argparse.ArgumentParser(prog='ccpbsa')

    cliparser.add_argument(
        "-w", "--wildtype",
        help=".pdb file of the wildtype protein"
    )
    cliparser.add_argument(
        "-m", "--mutations",
        help="a .txt file with the list of mutations"
    )
    cliparser.add_argument(
        "-f", "--flags",
        help="the flags, which should be passed to CONCOORD and GROMACS. At \
        least the number of structures in CONCOORD and the forcefield and water \
        model should be specified in there",
        default=pkgpath+'/parameters/flags.txt'
    )
    cliparser.add_argument(
        "--mode",
        help="Calculate either the change of stability or affinity changing \
        effects of mutations.",
        default='stability',
        choices={"stablitiy", "affinity"}
    )
    cliparser.add_argument(
        "--chains",
        help="Name of chains in the .pdb file.",
        default='A',
        nargs='+'
    )
    cliparser.add_argument(
        "--stability-parameters",
        help="scaling factors of the ddG stability calculations.",
        default=pkgpath+'/parameters/fit_stability.txt',
    )
    cliparser.add_argument(
        "--affinity-parameters",
        help="scaling factors of the ddG affinity calculations.",
        default=pkgpath+'/parameters/fit_affinity.txt',
    )
    cliparser.add_argument(
        "--minimization-mdp",
        default=pkgpath+'/parameters/min.mdp',
        help=".mdp file for GROMACS energy minimization. Default: Uses lookup \
        table in \"--dielectric-table\""
    )
    cliparser.add_argument(
        "--energy-mdp",
        default=pkgpath+'/parameters/energy.mdp',
        help=".mdp file for GROMACS energy evaluations. Default: epsilon-r = 1"
    )
    cliparser.add_argument(
        "--dielectric-table",
        default=pkgpath+'/parameters/table4r-6-12.xvg',
        help="Look up table for GROMACS to simulate continuum solvent. \
        Default: epsilon(r) = 4r"
    )
    cliparser.add_argument(
        "--gropbe",
        default=pkgpath+'/parameters/gropbe.txt',
        help="Parameters for gropbe. Defaults are the same as in the paper"
    )
    cliparser.add_argument(
        "--gxg-table",
        default=pkgpath+'/parameters/GXG.csv',
        help="GXG table used for stability change calculations"
    )
    cliparser.add_argument(
        "--gxg",
        help="Creates a new gxg.csv for dG calculation of the unfolded state.",
        action='store_true'
    )

    cliargs = cliparser.parse_args()

#    Very simple parser for reading in alpha, beta...
    if cliargs.mode == 'stability':
        with open(cliargs.stability_parameters, 'r') as fit:
            parameters = fit.readlines()
            parameters = [l[:-1] for l in parameters] # Remove newlines
            parameters = [l.split("=") for l in parameters]
            parameters = dict([(l[0], float(l[1])) for l in parameters])

    if cliargs.mode == 'affinity':
        with open(cliargs.affinity_parameters, 'r') as fit:
            parameters = fit.readlines()
            parameters = [l[:-1] for l in parameters] # Remove newlines
            parameters = [l.split("=") for l in parameters]
            parameters = dict([(l[0], float(l[1])) for l in parameters])

    if cliargs.gxg:
        print("Making a new GXG table.")
        GXG_ = gxg(
            cliargs.flags,
            cliargs.minimization_mdp,
            cliargs.energy_mdp,
            cliargs.dielectric_table,
            cliargs.gropbe
        )
        print(GXG_)
        sys.exit()

    data = DataGenerator(
        wt = cliargs.wildtype,
        mutlist = cliargs.mutations,
        flags = cliargs.flags,
        calculate = cliargs.mode,
        chains = cliargs.chains,
        min_mdp = cliargs.minimization_mdp,
        energy_mdp = cliargs.energy_mdp,
        mdrun_table = cliargs.dielectric_table,
        pbeparams = cliargs.gropbe
    )
    data.fullrun()
    search = DataCollector(data)
    search.search_data()
    print(search.ener_df)

    if cliargs.mode == 'stablility':
        search.ddstability(
            cliargs.gxg_table,
            parameters['alpha'],
            parameters['beta'],
            parameters['gamma'],
            parameters['tau']
        )

    if cliargs.mode == 'stablility':
        search.ddstability(
            parameters['alpha'],
            parameters['beta'],
            parameters['gamma'],
            parameters['c'],
            parameters['pka']
        )
