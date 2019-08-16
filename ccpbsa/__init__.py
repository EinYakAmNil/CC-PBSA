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
    cliparser = argparse.ArgumentParser(prog='ccpbsa', add_help=False)

    cliparser.add_argument(
        "routine",
        help="Called during installation to setup default parameters.",
        choices={'setup', 'stability', 'affinity'}
    )

    cliargs = cliparser.parse_args()

    if cliargs.routine == 'setup':

        options = cliparser.add_argument_group("OPTIONS")

        options.add_argument(
            "-h", "--help",
            help="Creates the default parameter file in the installation \
            directory",
            action='store_true'
        )

        if options.help:
            cliargs.print_help()
            sys.exit(0)
        
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
    if cliargs.routine == 'stability':
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
            "--chains",
            help="Name of chains in the .pdb file.",
            default='A',
            nargs='+'
        )
        options.add_argument(
            "--fit-parameters",
            help="scaling factors of the ddG stability calculations.",
            default=pkgpath+'/parameters/fit_stability.txt',
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

        if cliargs.v:
            verbose = 1
        else:
            verbose = 0

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

        with open(cliargs.fit_parameters, 'r') as fit:
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

    elif cliargs.routine == 'affinity':
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
            "--chains",
            help="Name of chains in the .pdb file for the first protein, \
                second protein group will be determined automatically.",
            default='A',
            nargs='+'
        )
        options.add_argument(
            "--fit-parameters",
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

        if cliargs.v:
            verbose = 1
        else:
            verbose = 0

        with open(cliargs.fit_parameters, 'r') as fit:
            parameters = fit.readlines()
            parameters = [l[:-1] for l in parameters] # Remove newlines
            parameters = [l.split("=") for l in parameters]
            parameters = dict([(l[0], float(l[1])) for l in parameters])

        search.daffinity()
        search.ddaffinity()

        print("G values:")
        print("bound")
        print(search.G_bound)
        search.G_bound.to_csv('G_bound.csv')
        print("unbound")
        print(search.G_grp1)
        search.G_grp1.to_csv('G_grp1.csv')
        print(search.G_grp2)
        search.G_grp2.to_csv('G_grp2.csv')

        print("dG values:")
        print("bound")
        print(search.dG_bound)
        search.dG_bound.to_csv('dG_bound.csv')
        print("unbound")
        print(search.dG_unbound)
        search.dG_bound.to_csv('dG_unbound.csv')

        print("ddG values:")
        print("without fit")
        print(search.ddG)
        search.ddG.to_csv('ddG.csv')

        search.fitaffinity(**parameters)
        print("with fit:")
        print(search.ddG)
        search.ddG.to_csv('ddG_fit.csv')
