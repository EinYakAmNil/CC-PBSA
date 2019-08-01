#import glob
import os
import shutil
import subprocess
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


def flatten(lst=[]):
    """Unpacks a nested list. Does not care about the deepness of nesting since
    it calls itself recursively.
    """
    unpacked = []
    for i in list(lst):
        if type(i) is not str and type(i) is not int:
            unpacked.extend(flatten(list(i)))

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
    ignc = [i if i[0] != "[" else i for i in ignc ]

#    Find where programs are defined.
    progs = [i for i in ignc if i[0] == "["]
    idx = [ignc.index(i) for i in progs]

#    Pack program flags together with their program...
    parsed = [ignc[i:idx[idx.index(i)+1]] for i in idx[:-1]]
    parsed.append(ignc[idx[-1]:])

#    ...and make them into a dictionary. For easy reference later on.
    parsed = [{i[0][1:-1]: flatten([j.split("=")  \
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


def concoord(verbosity=1, *pdb, **flags):
    """Takes arbitrarily many .pdb files as input to generate structure
    ensembles using CONCOORD. Additional flags can be passed via a dictionary
    with the keys corresponding to the programs, i.e. \"dist\" and \"disco\".
    Verbosity can be controlled via the keyword \"verbosity\" which takes
    either 0 (do not show any messages) or 1 (show all messages).
    Returns the two process objects of dist and disco.
    """
    if verbosity == 0:
        pipe = subprocess.PIPE

    elif verbosity == 1:
        pipe = None

    else:
        raise ValueError("Choose verbosity as an integer between 0 and 1.")

    for p in pdb:
        os.chdir("/".join([i for i in p.split("/")[:-1]]+["."]))
        
        dist_input = [
            'dist',
            '-p', p,
        ]
        if 'dist' in flags.keys():
            dist_input.extend(flags['dist'])

        disco_input = [
            'disco'
        ]
        if 'disco' in flags.keys():
            disco_input.extend(flags['disco'])
        
        if 'input' in flags.keys():
            input_ = flags['input']

        else:
            input_ = None

        dist = subprocess.run(
            dist_input,
            input=input_,
            stdout=pipe,
            stderr=pipe
        )
        disco =subprocess.run(disco_input, stdout=pipe, stderr=pipe)

        return dist, disco


def gmx(prog, **kwargs):
    """Run a GROMACS program with its flags by passing them in a list object.
    kwargs are passed to the subprocess.run method.
    """
    gmx = subprocess.run(['gmx', '-quiet'] + prog, **kwargs)

    return gmx


def minimize(
    pdb2gmx,
    editconf,
    grompp,
    mdrun,
    verbosity=1
):
    """Run a energy minimization starting from a .pdb file.
    """
    if verbosity == 0:
        pipe = {
            "stdout": subprocess.PIPE,
            "stderr": subprocess.PIPE
        }

    elif verbosity == 1:
        pipe = {
            "stdout": None,
            "stderr": None
        }

    else:
        raise ValueError("Choose verbosity as an integer between 0 and 1.")

    gmx(['pdb2gmx'] + pdb2gmx, **pipe)
    gmx(['editconf'] + editconf, **pipe)
    gmx(['grompp'] + grompp, **pipe)
    gmx(['mdrun'] + mdrun, **pipe)
 

def gro2pdb(gro, pdb, verbosity=1):
    if verbosity == 0:
        pipe = {
            "stdout": subprocess.PIPE,
            "stderr": subprocess.PIPE
        }

    elif verbosity == 1:
        pipe = {
            "stdout": None,
            "stderr": None
        }

    else:
        raise ValueError("Choose verbosity as an integer between 0 and 1.")

    trjconv = gmx(
        ['trjconv', '-f', gro, '-s', gro, '-o', pdb],
        input=b'0',
        **pipe
    )
    
    return trjconv


def filecheck(*strs):
    """Checks whether or not a string is a valid path to a file. Relative or
    absolute path does not matter. Yields the absolute path if the file was
    found
    """
    start = os.getcwd()

    for s in strs:
        splt = s.split("/")
        f = splt[-1]
        p = "/".join(splt[:-1]+["."])

        try:
            os.chdir(p)

            if f in os.listdir():
                yield s, os.getcwd() + "/" + f

        except FileNotFoundError:
            pass

        os.chdir(start)


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
    ):
        self.verb = 1
        self.wtpdb = os.getcwd() + "/" + wtpdb
        self.wt = wtpdb[:wtpdb.find(".pdb")]
        self.flags = parse_flags(flags)
        self.flags.setdefault("disco", []).extend(["-op", ""])
        self.mut_df = parse_mutations(mutlist)
        self.calculate = calculate
        self.chains = chains
        self.wds = [self.wt]
        self.wds.extend([", ".join([j for j in i if len(j) > 0]) \
            for i in self.mut_df.index])

#        Create the directory for data generation.
        try:
            os.mkdir(self.wt)
            os.chdir(self.wt)

        except FileExistsError:
            newdir = input("Directory already exists. Enter a new name:\n")
            os.mkdir(newdir)
            os.chdir(newdir)

        self.maindir = os.getcwd()
        os.mkdir(self.wt)
        shutil.copy(self.wtpdb, self.wt)

#        Copy all parameter files into the main directory. And update new file
#        locations.
        os.chdir('..')
        for k, v in self.flags.items():
            files = list(i for i in filecheck(*v))

            for ori, abs_ in files:
                shutil.copy(abs_, self.maindir)
                v[v.index(ori)] = abs_
                self.flags[k] = v

        os.chdir(self.maindir)

#        Create mutated structures during initialization of object. This is
#        done the be consistent with the working directories (wds) attribute
#        self.do_mutate()


    def __len__(self):
        """Returns the number of structures that should be generated by
        CONCOORD.
        """
        try:
            n_idx = self.flags['disco'].index('-n')
            return int(self.flags['disco'][n_idx+1])

        except ValueError:
            return 300


    def do_mutate(self):
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

            os.mkdir(self.wds[i+1])
            cmd.save(self.wds[i+1] + "/" + self.wds[i+1] + ".pdb")
            cmd.reinitialize()


    def do_concoord(self):
        """Goes into the directories listed in self.wds and generates
        structure ensmebles using CONCOORD. Each generated structure has its on
        directory one level down the directoy tree. Updates self.wds with the
        new structures
        """
        for d in self.wds:
            os.chdir(d)
            pdb = d.split("/")[-1] + ".pdb"
            concoord(self.verb, pdb, **self.flags, input=b'1\n1')

            for i in range(1, len(self)+1):
                os.mkdir(str(i))
                shutil.move(str(i) + '.pdb', str(i))

            os.chdir(self.maindir)

        self.wds = [d+'/'+str(i) for d in self.wds \
                for i in range(1, len(self)+1)]


    def chaindx(self, tpr='topol.tpr', **kwargs):
        """Create an .ndx file using an .tpr file to select chains for affinity
        calculation.
        """
        chainselec = [b'chain %b\n' % bytes(c, 'utf-8') for c in self.chains]
        chainselec = b''.join(chainselec) + b'q\n'
        ndx = gmx(
            ['make_ndx', '-f', tpr],
            input=chainselec,
            **kwargs
        )
        return ndx


    def do_minimization(self):
        """Do energy minimization on the structures in self.wds. If the
        affinity is to be calculated, then an index file for the chains will be
        made and the chains specified in self.chains are minimized.
        """
        for d in self.wds:
            os.chdir(d)
            pdb = d.split("/")[-1] + ".pdb"
            minimize(
                self.flags['pdb2gmx'] + ['-f', pdb],
                self.flags['editconf'],
                self.flags['grompp'],
                self.flags['mdrun']
            )

            if self.calculate == 'affinity':
                self.chaindx()
                gmx

            os.chdir(self.maindir)


    def update_structs(self):
        """Use gro2pdb to convert (minimized) confout.gro files into .pdb files
        with the same name as the directory.
        """
        for d in self.wds:
            os.chdir(d)
            gro2pdb(
                'confout.gro',
                d.split("/")[-1] + ".pdb",
                verbosity=self.verb
            )

            os.chdir(self.maindir)


    def electrostatics(self, tpr="sp"):
        """Calculate the Coulomb and Solvation Energies of structures of the
        specified .tpr file prefix in the directories. By default the file
        should be named sp (Single Point).
        Parameters for this are stored in a separate file for gropbe. Reference
        it through the main parameter file for CC/PBSA.
        """
        if self.verb == 1:
            pipe = {
                "stdout": subprocess.PIPE,
                "stderr": None,
            }

        else:
            pipe = {
                "stdout": subprocess.PIPE,
                "stderr": subprocess.PIPE,
            }
        chainselec = ",".join(str(i) for i in range(len(flatten(self.chains))))

        for d in self.wds:
            os.chdir(d)
            gropbe = subprocess.run(
                ["gropbe", self.flags["gropbe"][0]],
                input=bytes(chainselec),
                **pipe
            )

            if self.calculate == 'affinity':
                selidx = flatten(self.chains)
                selidx.sort()
                sel1 = ",".join([str(selidx.index(i)) for i in x.chains[0]])
                sel2 = ",".join([str(selidx.index(i)) for i in x.chains[1]])


    def fullrun(self):
        """Use all of the default behaviour to generate the data.
        """
        if self.calculate == 'stability':
            pass

        if self.calulate == 'affinity':
#            Temporarily disable minimization of chains.
            self.calculate = 'naffinity'
            self.do_minimization()
            self.do_concoord()
#            And turn it back on.
            self.calculate = 'affinity'
            self.do_minimization()
            self.electrostatics()
                


x = DataGenerator(
    "/home/linkai/test/1bxi.pdb",
    "mutations_1bxi.txt",
    "flags.txt",
    "affinity",
    [("A"), ("B")],
)
x.do_mutate()
x.do_minimization()
x.update_structs()
x.do_concoord()
x.do_minimization()
