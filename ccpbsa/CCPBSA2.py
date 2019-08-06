#import glob
import os
import shutil
import subprocess
import pandas as pd
import numpy as np
import pymol
pymol.finish_launching(['pymol', '-qc'])
cmd = pymol.cmd

def log(fname, proc_obj):
    """Write stdout and stderr of a process object to a file.
    """
    log_file = open(fname, 'a')

    if hasattr(proc_obj.stdout, 'decode'):
        log_file.write(proc_obj.stdout.decode('utf-8'))

    if hasattr(proc_obj.stderr, 'decode'):
        log_file.write(proc_obj.stderr.decode('utf-8'))

    log_file.close()


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


def concoord(pipe, *pdb, **flags):
    """Takes arbitrarily many .pdb files as input to generate structure
    ensembles using CONCOORD. Additional flags can be passed via a dictionary
    with the keys corresponding to the programs, i.e. \"dist\" and \"disco\".
    Verbosity can be controlled via the keyword \"verbosity\" which takes
    either 0 (do not show any messages) or 1 (show all messages).
    Returns the two process objects of dist and disco.
    """
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
    pipe
):
    """Run a energy minimization starting from a .pdb file.
    """
    gmx(['pdb2gmx'] + pdb2gmx, stdout=pipe, stderr=pipe)
    gmx(['editconf'] + editconf, stdout=pipe, stderr=pipe)
    gmx(['grompp'] + grompp, stdout=pipe, stderr=pipe)
    gmx(['mdrun'] + mdrun, stdout=pipe, stderr=pipe)
 

def gro2pdb(gro, tpr, pdb, pipe):
    trjconv = gmx(
        ['trjconv', '-f', gro, '-s', tpr, '-o', pdb],
        input=b'0',
        stdout=pipe,
        stderr=pipe
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
        spmdp,
        verbosity=0
    ):
        if verbosity == 0:
            self.pipe = subprocess.PIPE

        elif verbosity == 1:
            self.pipe = None

        else:
            raise ValueError

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

        shutil.copy(spmdp, self.maindir)
        self.spmdp = self.maindir + "/" + spmdp.split("/")[-1]
        print(self.spmdp)

        os.chdir(self.maindir)

#        Create mutated structures during initialization of object. This is
#        done the be consistent with the working directories (wds) attribute
        self.do_mutate()


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
            concoord(self.pipe, pdb, **self.flags, input=b'2\n1')

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
        chainselec = [b'chain %b\n' % bytes(c, 'utf-8') \
                for c in flatten(self.chains)]
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
                self.flags['mdrun'],
                pipe=self.pipe
            )

            if self.calculate == 'affinity':
                self.chaindx()

                for cn, c in enumerate(flatten(self.chains)):
                    gmx(
                        [
                            'editconf', '-f', 'conf.gro',
                            '-n', 'index.ndx',
                            '-o', 'chain%s.pdb' % c
                        ],
                        stdout=self.pipe,
                        stderr=self.pipe,
                        input=bytes(str(10+cn), 'utf-8')
                    )
                    minimize(
                        self.flags['pdb2gmx'] + \
                            [
                                '-f', 'chain%s.pdb' % c,
                                '-o', 'chain%s.gro' % c,
                                '-p', 'chain%s.top' % c,
                                '-i', 'chain%s.itp' % c
                            ],
                        self.flags['editconf'] + \
                                [
                                    '-f', 'chain%s.gro' % c,
                                    '-o', 'chain%s.gro' % c
                                ],
                        self.flags['grompp'] + \
                            [
                                '-c', 'chain%s.gro' % c,
                                '-p', 'chain%s.top' % c,
                                '-o', 'chain%s.tpr' % c
                            ],
                        self.flags['mdrun'] + \
                            [
                                '-deffnm', 'chain%s' % c,
                                '-c', 'chain%sout.gro' % c
                            ],
                        pipe=self.pipe
                    )

            os.chdir(self.maindir)


    def update_structs(self):
        """Use gro2pdb to convert (minimized) confout.gro files into .pdb files
        with the same name as the directory.
        """
        for d in self.wds:
            os.chdir(d)
            gro2pdb(
                'confout.gro',
                'topol.tpr',
                d.split("/")[-1] + ".pdb",
                pipe=self.pipe
            )

            os.chdir(self.maindir)


    def single_point(self):
        """Creates a single point .tpr file.
        """
        for d in self.wds:
            os.chdir(d)
            gmx([
                'grompp', '-f', self.spmdp,
                '-c', 'confout.gro',
                '-o', 'sp.tpr',
                '-maxwarn', '1'
            ])

            if self.calculate == 'affinity':

                for c in flatten(self.chains):
                    gmx([
                        'grompp', '-f', self.spmdp,
                        '-c', 'chain%sout.gro' % c,
                        '-o', 'chain%ssp.tpr' % c,
                        '-maxwarn', '1'
                    ])

            os.chdir(self.maindir)


    def electrostatics(self, tpr="sp"):
        """Calculate the Coulomb and Solvation Energies of structures of the
        specified .tpr file prefix in the directories. By default the file
        should be named sp (Single Point).
        Parameters for this are stored in a separate file for gropbe. Reference
        it through the main parameter file for CC/PBSA.
        """
        chainselec = ",".join(str(i) for i in range(len(flatten(self.chains))))

        for d in self.wds:
            os.chdir(d)
            shutil.copy(self.flags['gropbe'][0], 'gropbe.prm')

            with open('gropbe.prm', 'a') as params:
                params.write("in(tpr,sp.tpr)")

            gropbe = subprocess.run(
                ["gropbe", 'gropbe.prm'],
                input=bytes(chainselec, 'utf-8'),
                stdout=subprocess.PIPE,
                stderr=self.pipe
            )
            log("solvation.log", gropbe)

            if self.calculate == 'affinity':
                selidx = flatten(self.chains)
                selidx.sort()
                sel1 = ",".join([str(selidx.index(i)) for i in \
                    range(len(self.chains[0]))])
                sel2 = ",".join([str(selidx.index(i)) for i in \
                    range(len(self.chains[1]))])

                with open('gropbe.prm', 'a') as params:
                    params.append("in(tpr,sp.tpr)")

                gropbe = subprocess.run(
                    ["gropbe", self.flags['gropbe'][0]],
                    input=bytes(chainselec),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                log("solvation.log", gropbe)

                with open('gropbe.prm', 'a') as params:
                    params.append("in(tpr,sp.tpr)")

                gropbe = subprocess.run(
                    ["gropbe", self.flags['gropbe'][0]],
                    input=bytes(chainselec),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                log("solvation.log", gropbe)

            os.chdir(self.maindir)


    def lj(self):
        """Calculate the Lennard-Jones Energy based on a sp.tpr
        """
        for d in self.wds:
            os.chdir(d)

            gmx([
                'mdrun', '-s', 'sp.tpr',
                '-rerun', 'confout.gro',
                '-deffnm', 'sp'
            ])
            lj = gmx(
                ['energy', '-f', 'sp.edr', '-sum', 'yes'],
                input=b'5 7',
                stdout=subprocess.PIPE,
                stderr=None
            )
            log("lj.log", lj)

            os.chdir(self.maindir)


    def area(self):
        """Calculate the solvent accessible surface area and saves it to
        area.xvg. If the mode is set to affinity, only the wt protein structure
        ensemble will be used and the values for the interaction surface will
        be written into the .xvg file
        """
        sasa = ['sasa', '-s', 'confout.gro']

        for d in self.wds:
            os.chdir(d)

            if self.mode == 'affinity':
                gmx(sasa + ['-n', 'index.ndx', '-output', '10', '11'],
                    input=b'0')

            else:
                gmx(sasa, input=b'0')

            os.chdir(self.maindir)


    def schlitter(self):
        """Calculates an upper limit of the entropy according to Schlitter's
        formula. Used in .fullrun() if the mode is stability
        """
        trjcat = ['trjcat', '-cat', 'yes', '-f']
        covar = ['covar', '-f', 'trajout.xtc', '-nofit', '-nopbc', '-s']
        anaeig = ['anaeig', '-v', 'eigenvec.trr', '-entropy']
        ensembles = [i.split('/')[-2] for i in self.wds]

        for en in ensembles:
            os.chdir(en)
            trrs = [self.maindir+'/'+en+'/%d/traj.trr' % (n+1)
                for n in range(len(self))]

            gmx(trjcat+trrs)
            gmx(covar+[self.maindir+'/'+en+"/1/confout.gro"], input=b'0')
            entropy = gmx(anaeig, stdout=subprocess.PIPE)
            log('entropy.log', entropy)

            os.chdir(self.maindir)


    def fullrun(self):
        """Use all of the default behaviour to generate the data.
        """
        if self.calculate == 'stability':
            self.do_minimization()
            self.update_structs()
            self.do_concoord()
            self.do_minimization()
            self.single_point()
            self.electrostatics()

        if self.calculate == 'affinity':
#            Temporarily disable minimization of chains.
            self.calculate = 'naffinity'
            self.do_minimization()
#            And turn it back on.
            self.calculate = 'affinity'
            self.update_structs()
            self.do_concoord()
            self.do_minimization()
            self.single_point()
            self.electrostatics()
                

x = DataGenerator(
    "1ayi.pdb",
    "mutations_1ayi.txt",
    "flags.txt",
    "stability",
    [("A",)],
    "/home/linkai/CC-PBSA/ccpbsa/parameters/energy.mdp",
    verbosity=1
)
x.fullrun()
