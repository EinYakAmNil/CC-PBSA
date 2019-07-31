#import glob
import os
import shutil
#import subprocess
import pandas as pd
import numpy as np
import pymol
pymol.finish_launching(['pymol', '-qc'])
cmd = pymol.cmd

def int_in_str(*strings):
    """Takes strings as input and tries to find integers in them. Used to find
    the residue number in the mutation file.
    Returns a list of the integers found.
    """
    converted = [''.join([item for item in str_ if item.isdigit()]) \
        for str_ in strings]
    return converted
#    return list(map(int, converted))


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


def parse_flags(file_):
    """Parse the file containing the flags for CONCOORD and GROMACS. Returns a
    list of dictionaries with each program as key for the flags.
    """
    raw = open(file_, 'r').readlines()
    
#    Strip lines of comments.
    ignc = [i[:i.find(";")] for i in raw]
    ignc = [i.strip() for i in ignc if len(i) > 0]
    ignc = ["-"+i if i[0] != "[" else i for i in ignc ]

#    Find where programs are defined.
    progs = [i for i in ignc if i[0] == "["]
    idx = [ignc.index(i) for i in progs]

#    Pack program flags together with their program...
    parsed = [ignc[i:idx[idx.index(i)+1]] for i in idx[:-1]]
    parsed.append(ignc[idx[-1]:])

#    ...and make them into a dictionary. For easy reference later on.
    parsed = [{i[0][1:-1]: unpack([j.split("=")  \
        for j in i[1:]])} for i in parsed]

    return dict(i for dict_ in parsed for i in dict_.items())


def parse_mutations(file_):
    """Returns a 3 dimensional pandas DataFrame of Mutations. Containing:
        - chain
        - amino acid (AA)
        - residue(number)
        - new amino acid (Mutation)
    Each new protein stores the information in a new row.
    Multiple mutations extend into the third dimension.
    """
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU \
            MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    aa123 = dict(zip(aa1,aa3))
    raw = open(file_, 'r').readlines()

#    Remove whitespaces. Like BLM tries to do.
    raw = [i.replace(" ", "") for i in raw]
    raw = [i[:-1].split(",") for i in raw]
    maxmuts = max([len(i) for i in raw])
    muts = [tuple(i + [""] * (maxmuts - len(i))) for i in raw]

#    Index definition for the DataFrame.
    n_muts = list(np.arange(len(raw)))
    headers = ["Chain", "AA", "Residue", "Mutation"]

    mut_df = pd.DataFrame(
        columns=headers,
        index=pd.MultiIndex.from_tuples(muts)
    )

#    Fill the Chain column.
    for i, j in enumerate(raw):
        ext = maxmuts - len(j)

        chains = [k.split("_")[0] if len(k.split("_")) == 2 else "" for k in j]
        chains.extend([""] * ext)
        mut_df["Chain"][i] = chains

        aa = [k.split("_")[0][0] if len(k.split("_")) == 1 \
                else k.split("_")[1][0] for k in j]
        aa.extend([""] * ext)
        mut_df["AA"][:][i] = aa

        res = int_in_str(*j)
        res.extend([""] * ext)
        mut_df["Residue"][i] = res

        muts = [aa123[k[-1]] for k in j]
        muts.extend([""] * ext)
        mut_df["Mutation"][i] = muts

    return mut_df


class DataGenerator:
    """Main class to make CC/PBSA work. Creates a directory with the name of
    the wildtype (wt) protein and subdirectories for each structure ensemble.
    """
    def __init__(
        self,
        wtpdb, # wild type .pdb file
        mutlist, # list of mutations
        flags, # file specifying the flags for CONCOORD and GROMACS
        calculate, # Output either stability or affinity
        chains, # specify chains for affinity calculation
        min_mdp,
        mdrun_table,
        pbeparams
    ):
        self.wtpdb = os.getcwd() + "/" + wtpdb
        self.wt = self.wtpdb[:self.wtpdb.find(".pdb")]
        self.flags = parse_flags(flags)
        self.mut_df = parse_mutations(mutlist)
        self.calculate = calculate
        self.chains = chains
        self.min_mdp = min_mdp
        self.mdrun_table = mdrun_table
        self.pbeparams = pbeparams
        self.dirs = [self.wt]
        self.dirs.extend([", ".join([j for j in i if len(j) > 0]) \
            for i in self.mut_df.index])

#        os.mkdir(self.wt)
#        os.chdir(self.wt)
#        os.mkdir(self.wt)
#        self.maindir = os.getcwd()
#        shutil.copy("../" + self.wtpdb, self.wt)


    def mutate(self):
        """Create directories for each mutation and save the .pdb file in
        there. Requires the mut_df attribute.
        """
        for i in range(len(self.mut_df.index)):
            cmd.load(self.wtpdb)
            cmd.wizard('mutagenesis')

            for j in range(len(self.mut_df["Residue"][i])):
                cmd.get_wizard().do_select('///%s/%s' % (
                    self.mut_df["Chain"][i][j],
                    self.mut_df["Residue"][i][j],
                ))
                cmd.get_wizard().set_mode(self.mut_df["Mutation"][i][j])
                cmd.get_wizard().apply()

            os.mkdir(self.dirs[i+1])
            cmd.save(self.dirs[i+1] + "/" + self.dirs[i+1] + ".pdb")
            cmd.reinitialize()


    def concoord(self):
        for fpf in [1, 2, 3]:
            dist_input = [
                'dist',
                '-p', '%s' % fpf+'.pdb',
                '-op', '%s_dist.pdb' % fpf,
                '-og', '%s_dist.gro' % fpf,
                '-od', '%s_dist.dat' % fpf,
            ]

            disco_input = [
                'disco',
                '-d', '%s_dist.dat' % fpf,
                '-p', '%s_dist.pdb' % fpf,
                '-op', '',
                '-or', '%s_disco.rms' % fpf,
                '-of', '%s_disco_Bfac.pdb' % fpf
            ]

            dist_input.extend(self.flags["dist"])
            disco_input.extend(self.flags["disco"])


x = DataGenerator(
    "1ayi.pdb",
    "muts",
    "flags.txt",
    0,
    0,
    0,
    0,
    0
)
print(x.mut_df)
x.mutate()
