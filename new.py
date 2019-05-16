"""Code for CC/PBSA, a fast tool for estimating mutational free energy
differences.
"""
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
import pymol
pymol.finish_launching(['pymol', '-qc'])
cmd = pymol.cmd


def gmx(prog, **kwargs):
    """Uses the subprocess module to run a gmx (GROMACS) program. Returns the
    process object. **kwargs will be passed to subprocess.run
    """
    assert type(prog) == list, "Pass a list of arguments you would use after \
        \"gmx -quiet\"."
    gmx = subprocess.run(['gmx', '-quiet'] + prog, **kwargs)

    return gmx


def makedir(dirname):
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    else:
        print("Directory \"%s\" already exists. Ignoring this function call!" %
            dirname)
        pass


def unpack(lst=[]):
	"""Unpacks a nested list. Does not care about the deepness of nesting since
		it calls itself recursively.
	"""
	unpacked = []
	for i in lst:
		if type(i) is list:
			unpacked.extend(unpack(i))

		else:
			unpacked.append(i)

	return unpacked


class DataGenerator:
    """Main class containing the methods to generate the data for CC/PBSA.
    Creates a directory with the name of the wildtype (wt) protein and
    subdirectories for each structure ensemble (wt and mutant). If the concoord
    method was called, then each subdirectory will have another layer of
    subdirectories for the individual GROMACS runs.The user will have to pass
    the wt protein .pdb file, a list of mutations, a parameter file containing
    the flags for each programm (CONCOORD/GROMACS) and optionally user specific
    .mdp files and tables for GROMACS.
    """
    def __init__(self, wt, mutlist, flags, calculate='stability',
        min_mdp="/Users/linkai/CC_PBSA/min.mdp",
        energy_mdp="/Users/linkai/CC_PBSA/energy.mdp",
        mdrun_table="/Users/linkai/CC_PBSA/table4r-6-12.xvg"):
        """Creates and moves to the main folder upon initialization and copies
        the wt .pdb file into a subdirectory of the working directory.
        self.maindir will be the directory to which each function call returns
        to after completing. Structure (ensembles) will be tracked in
        self.sembles
        """
        self.mode = calculate
        self.wt = wt
        self.flags = self.parse_flags(flags)
        self.mut_df = self.parse_mutlist(mutlist)
        self.e_mdp = energy_mdp
        self.min_mdp = min_mdp
        self.mdrun_table = mdrun_table
        pass


    def parse_flags(self, flags_raw):
        """Parse the list of flags in the .txt flag file for the programs.
        Called upon initialization.
        """
        pass


    def parse_mutlist(self, mutlist_raw):
        """Parse the list of mutations so that .mutate() can perform mutations
        using PyMOL correctly. Called upon initialization. Adds mutants to
        self.sembles.

        The format of the mutation instruction should be (in one letter code
        and without spaces):
            (*Chain_*) *OriginalAA* *ResidueNumber* *NewAA*
        for exapmle:
            - A20G (for a monomer)
            - B_H10I (for a dimer with a chain named \"B\")
        """
        def int_in_str(lst):
            converted = ''.join([item for item in lst if item.isdigit()])
            return int(converted)


#        One letter AA code to three letter code as presented on PyMOL Wiki.
        aa1 = list("ACDEFGHIKLMNPQRSTVWY")
        aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU \
            MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
        aa123 = dict(zip(aa1,aa3))

        mut_lst = list(open(mutlist_raw, 'r'))
        mut_lst = [mut.split('\n')[0] for mut in mut_lst]

#        Returned object is a pandas dataframe.
        parsed_mut = pd.DataFrame(
            columns=["Chain", "AA", "Residue", "Mutation"], index=mut_lst)

        parsed_mut['Chain'] = [i[0] if '_' in i else '' for i in mut_lst.copy()]
        parsed_mut['AA'] = [i[0] if '_' not in i else i[i.index('_')+1] \
            for i in mut_lst.copy()]
        parsed_mut['Residue'] = [int_in_str(i) for i in mut_lst]
        parsed_mut['Mutation'] = [aa123[i[-1]] for i in mut_lst]

        return parsed_mut


    def mutate(self):
        """Uses PyMOL and the list of mutations to mutate the wildtype protein.
        Does not accuratly mutate if input structure or mutation instructions
        are flawed. WARNING: No specific message is given if that happens. Best
        to check if the residues in the .pdb file are correctly numbered.
        """
        pass


    def log(self, fname, proc_obj):
        """Used to log some of GROMACS output
        """
        log_file = open(fname, 'a')
        log_file.write(proc_obj.stdout.decode('utf-8'))
        log_file.close()

    
    def concoord(self):
        """Performs the CONCOORD procedure to generate protein structure
        ensembles. Takes additional flags from "flag_parse" as input (pass it
        as "flag_parse_output['CC']"). Make sure that \"CONCOORDRC.bash\" is
        sourced.
        """
        pass


    def minimize(self):
        """Minimzes all structures in self.structures
        """
        pass


    def electrostatics(self):
        pass


    def lj(self):
        """calculates the single point Lennard Jones Energy (1-4 and
        shortrange) of all self.structures
        """
        pass


    def area(self):
        """Calculate the solvent accessible surface area and saves it to
        area.xvg. If the mode is set to affinity, only the wt protein structure
        ensemble will be used and the values for the interaction surface will
        be written into the .xvg file
        """
        pass


    def schlitter(self):
        """Calculates an upper limit of the entropy according to Schlitter's
        formula. Used in .fullrun() if the mode is stability
        """
        pass


    def fullrun(self):
        """Performs a full run based on the list of mutations, number of
        CONCOORD structures and other parameters. Methods used will depend on
        the mode chosen (stability/affinity).
        """
        pass


if __name__ == '__main__':
    x = DataGenerator("1pga.pdb", "mut.txt", "param.txt")
    print(x.mut_df)
