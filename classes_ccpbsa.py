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


def makedir(dirname):
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    else:
        print("Directory \"%s\" already exists. Ignoring this function call!" %
            dirname)
        pass


def unpack(lst):
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


class Process_Logger:
    """Very simple logger, which only writes stdout of given subprocess to a
    file. Inherited by all classes which log something.
    """
    def __init__(self, logname):
        """self. logname is the name of the file the stdout will be printed to.
        self.process is the object of the process to be logged.
        """
        self.logname = logname
        refresh = open(self.logname, 'w')
        refresh.close()

    def log(self, proc_obj):
        log = open(self.logname, 'a')
        log.write(proc_obj.stdout.decode('utf-8'))
        log.close()
    

class Mutater:
    """Takes a .pdb file and mutates one residue using PyMOL. Can take a .txt
    file with a list of mutations to generate many versions.
    """
    def __init__(self, pdb_file, txt):
        assert type(res) is list, "res, mut and chain must be of type list."
        assert type(mut) is list, "res, mut and chain must be of type list."
        assert type(chain) is list, "res, mut and chain must be of type list."
        self.wt = pdb_file
        self.mut_lst = None
        self.chain = None
        self.aa = None
        self.res = None
        self.mut = None
        self.parse_mutations(txt)


    def mutate(self):
        """Uses PyMOL and the list of mutations to mutate the wildtype protein.
        """
        cmd.wizard('mutagenesis')

        for m in range(len(self)):

            cmd.get_wizard().do_select('///%s/%s' %
                (self.chain[m], self.res[m]))
            cmd.get_wizard().set_mode(self.mut[m])
            cmd.get_wizard().apply()
            cmd.save("%s_%s.pdb" % (self.wt[:-4], self.mut_lst[m]))
            cmd.reinitialize()
        

    def parse_mutations(self, txt):
        """Modifies the attributes of this object to make them passable to
        \"mutate()\".
        """
        def int_in_str(lst):
            converted = ''.join([item for item in lst if item.isdigit()])
            return int(converted)


#        One letter AA code to three letter code as presented on PyMOL Wiki.
        aa1 = list("ACDEFGHIKLMNPQRSTVWY")
        aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU \
            MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
        aa123 = dict(zip(aa1,aa3))
        
        self.mut_lst = list(open(mutations_txt, 'r'))
        self.mut_lst = [mut.split('\n')[0] for mut in self.mut_lst]
        self.chain = [i[0] if '_' in i else '' for i in self.mut_lst.copy()]
        self.aa = [i[0] if '_' not in i else i[i.index('_')+1] \
            for i in self.mut_lst.copy()]
        self.res = [int_in_str(i) for i in self.mut_lst]
        self.mut = [i[-1] for i in self.mut_lst]


    def __len__(self):
        """Returns the number of mutations stored in \"mut_lst\"
        """
        return len(self.mut_lst)

        
class FileOrganizer:
    """This class can split up GROMACS trajectory files into their seperate
    directories using \"gmx trjconv\". Also contains methods to search for
    files either exactly or partly matching the entered string.
    Many methods are relative to the attribute \"cwd\".
    """
    def __init__(self, wd=os.getcwd()):
        self.cwd = wd
        self.files_split = 0
        
        try:
            os.chdir(self.cwd)
        
        except FileNotFoundError as err:
            raise type(err)("Could not change to directory \"%s\"" % self.cwd)


    def find(self, file_, match=False):
        """This method will try to find a file matching either part or exactly
        the input relative to the current working directory.
        It searches from the current directory up to all parent and
        child directories (no sibling trees). If match is set to True only the
        first result matching the filename exactly will be returned.
        Search Order:
        1. First searches from the current directory up to its parents.
        2. Then goes down to all the child directories.
        """
        floc = [self.cwd+"/"+f for f in os.listdir() if f == file_]

        if len(floc) > 0 and match:
            return floc[0]

        parents = [i for i in self.cwd.split('/')][1:]
        parents[0] = '/' + parents[0]

        for i in range(1, len(parents)):
            parents[i] = parents[i-1] + '/' + parents[i]

#        Invert list to start search from current directory.
        parents = parents[::-1]

        for i in parents:

            if file_ in os.listdir(i):

                if match:
                    floc = i + '/' + file_
                    return floc

                else:
                    floc.append(i + '/' + file_)

        floc.extend([f for f in self.list_children() if file_ in f])

        if len(floc) == 0:
            raise FileNotFoundError('Could not find file %s' % file_)

#        Special case for when user wants only the string to one file matching
#        the input exactly.
        if match:

            try:
#                Any matching path will be stripped down to its filename again
                files = [''.join(i.split('/')[-1]) for i in floc.copy()]
#                And the first element matching the input filename is returned
                return floc[files.index(file_)]

#            Make the error message more understandable, since it is not
#            actually a IndexError.
            except IndexError:
                raise FileNotFoundError('Could not find file %s' % file_)

        else:
#            Remove the list of any duplicates and returns it
            return list(dict.fromkeys(floc))



    def list_children(self):
        """List all files with the full path of a tree.
        """
        files = [os.getcwd()+'/'+i for i in os.listdir() if os.path.isfile(i)]
        dirs = [i for i in os.listdir() if os.path.isdir(i)]

        for dir_ in dirs:
            os.chdir(dir_)
            files.extend(self.list_children())
            os.chdir('..')

        return files
        

    def parse_proc_flags(self, flags_raw):
        """Takes a file with flags for \"gmx\" and \"CONCOORD\" programs as
        input. Returns a dictionary in the following format:

        parsed_flags = {
            'CC': {
                'DIST FLAGS': [list of flags],
                'DISCO FLAGS': [list of flags]
            },
            'gmx': {
                'PDB2GMX FLAGS': [list of flags],
                'EDITCONF FLAGS': [list of flags],
                'GROMPP FLAGS': [list of flags],
                'MDRUN FLAGS': [list of flags],
            }
        }
        
        The list can then be extended directly to other input lists for
        subprocess.run().
        """
        parsed_flags = {
            'CC': {
                'DIST FLAGS': [],
                'DISCO FLAGS': []
            },
            'gmx': {
                'PDB2GMX FLAGS': [],
                'EDITCONF FLAGS': [],
                'GROMPP FLAGS': [],
                'MDRUN FLAGS': [],
            }
        }

#        Search for this file just in case it is not in the current directory of
#        the class.
        flags_raw = self.find(flags_raw, match=True)
        flag_file = list(open(flags_raw, 'r'))
        content = [line[:line.index('\n')] for line in flag_file]
        uncommented = []

        for line in content:
            if ';' in line:
                line = line[:line.index(';')]
    
            if len(line) > 0:
                uncommented.append(' '.join(line.split()))
    
    #    Indexes the file for where the flags are defined
        idx = [uncommented.index(j) for i in parsed_flags.values() for j in i]
        idx.append(len(uncommented))
    
        i = 0
        for keys, vals in parsed_flags.items():
    
            for prog in vals:
                parsed_flags[keys][prog] = unpack([i.split('=') for i in \
                    uncommented[idx[i]+1:idx[i+1]]])
                i += 1
    
        return parsed_flags


class CC(Process_Logger):
    """Interacts with the programs of CONOCOORD (\"from CONstraints to
    COORDinates\"). Also parses the \"flags.txt\" file for gmx programs.
    """
    def __init__(self, pdb_files, flags):
        """The attribute flag_file is the file containing the flags for all
        following programs. n_structures will be 0 until the \"concoord()\"
        method has been called.
        """
        self.flags = flags
        self.n_structures = 0
        self.pdb_files = unpack([pdb_files])
        self.logname = "CONCOORD.log"
        self.start = os.getcwd()

    def concoord(self):
        """Performs the CONCOORD procedure to generate protein structure
        ensembles. Takes additional flags from "flag_parse" as input (pass it as
        "flag_parse_output['CC']"). Make sure that \"CONCOORDRC.bash\" is
        sourced.
        """
        for f in self.pdb_files:
            fname = f[:-4]
            path = ('/'.join(fname.split('/')[:-1]))
            os.chdir(path)
    
            dist_input = [
                'dist',
                '-p', '%s' % f,
                '-op', '%s_dist.pdb' % fname,
                '-og', '%s_dist.gro' % fname,
                '-od', '%s_dist.dat' % fname,
            ]
    
            disco_input = [
                'disco',
                '-d', '%s_dist.dat' % fname,
                '-p', '%s_dist.pdb' % fname,
                '-op', '',
                '-or', '%s_disco.rms' % fname,
                '-of', '%s_disco_Bfac.pdb' % fname
            ]
            
            dist_input.extend(self.flags['DIST FLAGS'])
            disco_input.extend(self.flags['DISCO FLAGS'])
            self.n_structures = int(disco_input[disco_input.index('-n')+1])
    
            dist_proc = subprocess.run(
                dist_input,
                input=b'1\n1',
                stdout=subprocess.PIPE
            )
            self.log_CC(dist_proc)
            disco_proc = subprocess.run(
                disco_input,
                stdout=subprocess.PIPE
            )
            self.log_CC(disco_proc)
    
            for struct in range(self.n_structures):
                nr = str(struct+1)
                os.mkdir(path+'/'+nr)
                shutil.move(nr+'.pdb', nr)
            
            print("Generated %d structures for %s" % (self.n_structures, f))

        os.chdir(self.start)


    def __len__(self):
        """Returns the number of input structures used to initialize this
        object.
        """
        return len(self.pdb_files)


    def log_CC(self, proc_object):
        """Logs stout of all CONCOORD subprocesses into \"CONCOORD.log\".
        """
        log = open(self.logname, 'a')
        log.write(proc_object.stdout.decode('utf-8'))
        log.close()


class Energy_Minimizer(Process_Logger):
    """This class contains the functions and attributes to perform energy
    minimizations on structure (ensembles).
    """
    def __init__(
        self, pdb_files, mdp_file,
        gmx_flags, logname="gmx_minimization.log"
    ):
        """Make a class using a single structure or a structure ensemble saved
        as single structures as a list. Also the .mdp file for energy
        minimization and flags for the gmx programs which don't involve naming
        in-/output need to be specified.
        """
        super().__init__(logname)
        self.pdb_files = unpack([pdb_files])
        self.fn = [name[:-4] for name in self.pdb_files]
        self.mdp_file = str(mdp_file)
        self.gmx_flags = gmx_flags
        self.gmx_procedure = {
            'pdb2gmx_input': [
                'gmx', '-quiet', 'pdb2gmx',
                '-f', 'FILENAME.pdb',
                '-o', 'FILENAME.gro' ,
                '-p', 'FILENAME.top' 
        ],
            'editconf_input': [
                'gmx', '-quiet', 'editconf',
                '-f', 'FILENAME.gro',
                '-o', 'FILENAME_box.pdb' 
        ],
            'grompp_input': [
                'gmx', '-quiet', 'grompp',
                '-f', mdp_file,
                '-c', 'FILENAME_box.pdb' ,
                '-p', 'FILENAME.top' ,
                '-o', 'FILENAME.tpr' 
        ],
            'mdrun_input': [
                'gmx', '-quiet', 'mdrun',
                '-s', 'FILENAME.tpr' ,
                '-deffnm', 'FILENAME_min' ,
                '-nt', '1',
                '-table', 'table4r-6-12.xvg',
                '-tablep', 'table4r-6-12.xvg'
        ]}


    def minimize(self):
        """Method that does the actual energy minization. Calls the \"log()\"
        method to record the GROMACS terminal output. Shows a progress bar based
        on the number of structures passed to this method.
        """
        for f in range(len(self)):

            x = 0
            for fname in self.gmx_procedure.copy().values():
                proc_input = [name.replace('FILENAME', self.fn[f]) for name in fname]
                proc_input.extend(list(self.gmx_flags.values())[x])
                proc = subprocess.run(proc_input, stdout=subprocess.PIPE)
                self.log(proc)

                if proc.returncode > 0:
                    raise Exception("Something went wrong!")

                x += 1


    def __len__(self):
        """Returns the number of structures in self.min_gro.
        """
        return len(self.pdb_files)


class LJ_Energy(Energy_Minimizer):
    """This class contains the output from GROMACS when running an energy
    minimzation and extracting the energy from a single point MD-run.
    """
    def __init__(self, pdb_files, mdp_file, gmx_flags, logname="gmx_energy.log"):
        self.logname = logname
        self.lj = 0
        super().__init__(pdb_files, mdp_file, gmx_flags, logname=self.logname)

        self.gmx_procedure = {
            'pdb2gmx_input': [
                'gmx', '-quiet', 'pdb2gmx',
                '-f', 'FILENAME.gro',
                '-o', 'FILENAME_energy.gro',
                '-p', 'FILENAME_energy.top' 
        ],
            'editconf_input': [
                'gmx', '-quiet', 'editconf',
                '-f', 'FILENAME_energy.gro',
                '-o', 'FILENAME_energy_box.pdb' 
        ],
            'grompp_input': [
                'gmx', '-quiet', 'grompp',
                '-f', mdp_file,
                '-c', 'FILENAME_energy_box.pdb',
                '-p', 'FILENAME_energy.top',
                '-o', 'FILENAME_energy.tpr' 
        ],
            'mdrun_input': [
                'gmx', '-quiet', 'mdrun',
                '-s', 'FILENAME_energy.tpr',
                '-deffnm', 'FILENAME_energy',
                '-nt', '1'
        ]}
        self.energy_input = [
            'gmx', '-quiet', 'energy',
            '-f', 'FILENAME_energy.edr',
            '-o', 'FILENAME_energy.xvg'
        ]

    def calc_lj(self):
        """Calculate the Lennard-Jones energy (1-4 and SR) using the arguments to
        initialize this structure (ensemble). Uses the energy minimization
        method of Energy_Minimizer before extracting the energies.
        """
        self.minimize()

        for f in range(len(self)):
            proc = subprocess.run(
                [i.replace('FILENAME', self.fn[f]) \
                for i in self.energy_input.copy()],
                input=b'5 7', stdout=subprocess.PIPE
            )
            self.log(proc)

            if proc.returncode > 0:
                raise Exception("Something went wrong!")
            
            energy = proc.stdout.decode('utf-8').split()

            shortrange = energy[energy.index('(SR)')+1]
            onefour = energy[energy.index('LJ-14')+1]
            self.lj += float(onefour) + float(shortrange)

        return self.lj/len(self)


class Schlitter_Entropy(Process_Logger):
    """Calulate the Entropy of a .xtc file according to Schlitter's formula.
    Has the method \"concat\" to concatenate .trr files regardless of time into
    a trajectory file.
    """
    def __init__(self, logname="schlitter.log"):
        """Choose name of the logfile here.
        """
        super().__init__(logname)

    def calc(self, xtc_file, topol_file):
        """Calculate the entropy according to Schlitter's formula
        """
        eigenvec = xtc_[:-4] + '.trr'

        entropy = {
            'covar_input': [
                'gmx', '-quiet', 'trjcat',
                '-f', xtc,
                '-pbc', 'no',
                '-fit', 'no',
                '-v', eigenvec
            ],
            'anaeig_input': [
                'gmx', '-quiet', 'trjcat',
                '-v', eigenvev,
                '-s', topol_file,
                '-entropy', 'no'
            ]
        }
        proc = subprocess.run(entropy['covar_input'], input=b'0')
        self.log(proc) 
        proc = subprocess.run(entropy['anaeig_input'], stdout=subprocess.PIPE)
        self.log(proc)
        schlitter = proc.stdout.decode('utf-8').split('\n')[0]
        schlitter = schlitter[schlitter.index('is ')+len('is ') \
            :schlitter.index(' J')]
        
        return float(schlitter)
        

    def concat(self, trr_files):
        """Use \"gmx trjcat\" to concatenate .trr files of MD-runs.
        """
        trjcat_input = [
            'gmx', '-quiet', 'trjcat',
            '-cat', '-f'
        ]
        trjcat_input.extend(trr_files)
        proc = subprocess.run(trjcat_input)
        self.log(proc)
        
        if proc.returncode > 0:
            raise Exception("Could not concatenate files. Check your files")


class Protein:
    """Contains all the attributes and methods to create the CC/PBSA workflow.
    Requires the programs \"CONCOORD\" and \"GROMACS 2019.2\". Python 3.7
    module dependencies are numpy, pandas and python3-pymol.
    """
    def __init__(self, pdb_file, mut_list, flags,
        min_mdp="/home/linkai/CC_PBSA/min.mdp",
        energy_mdp="/home/linkai/CC_PBSA/energy.mdp", calculate='stability',
        mdrun_table="/home/linkai/CC_PBSA/table4r-6-12.xvg"):
        """Mandatory inputs:
        - The protein in a .pdb file to be used.
        - List of single amino acid (aa) mutations to intoduce to the protein
        - A .txt file for mandatory and optional flags depending of the
          program.
        Optional inputs:
        - Parameter files to be used for structure generation, energy
          minimization, energy calculation, user-defined tables for \"gmx
          mdrun\", surface area calculation and the entropy according to
          Schlitter's formula.
        - Mode of calculation: \"stability\" or \"affinity\".
        Description of the attributes:
        - mut_df is a pandas dataframe
        """
        self.pdb_file = pdb_file
        self.structures = [pdb_file]
        self.name = pdb_file[:-4]
        self.mut_df = self.parse_mutations(mut_list)
        self.parsed_flags = self.parse_flags(flags)
        self.n_structs = int(self.parsed_flags['CC']['DISCO FLAGS'] \
                [self.parsed_flags['CC']['DISCO FLAGS'].index('-n')+1])
        self.e_mdp = energy_mdp
        self.min_mdp = min_mdp
        self.mdrun_table = mdrun_table
#        self.energy_df
#        self.cpt
        
        makedir(self.name)
        os.chdir(self.name)
        makedir(self.name)
        shutil.copy("../" + self.pdb_file, self.name)
        self.maindir = os.getcwd()
        

    def __repr__(self):
        """Returns the filename of the .pbd file. This string will also be used
        for all protein mutations as prefix.
        """
        return self.name


    def parse_flags(self, flags_raw):
        """This function will be called upon object creation. Takes a file
        with flags for \"gmx\" and \"CONCOORD\" programs as input. Returns a
        dictionary in the following format:

        parsed_flags = {
            'CC': {
                'DIST FLAGS': [list of flags],
                'DISCO FLAGS': [list of flags]
            },
            'gmx': {
                'PDB2GMX FLAGS': [list of flags],
                'EDITCONF FLAGS': [list of flags],
                'GROMPP FLAGS': [list of flags],
                'MDRUN FLAGS': [list of flags],
            }
        }
        
        The list can then be extended directly to other input lists for
        subprocess.run().
        """
        parsed_flags = {
            'CC': {
                'DIST FLAGS': [],
                'DISCO FLAGS': []
            },
            'gmx': {
                'PDB2GMX FLAGS': [],
                'EDITCONF FLAGS': [],
                'GROMPP FLAGS': [],
                'MDRUN FLAGS': [],
            }
        }

#        Search for this file just in case it is not in the current directory of
#        the class.
        flag_file = list(open(flags_raw, 'r'))
        content = [line[:line.index('\n')] for line in flag_file]
        uncommented = []

        for line in content:
            if ';' in line:
                line = line[:line.index(';')]
    
            if len(line) > 0:
                uncommented.append(' '.join(line.split()))
    
    #    Indexes the file for where the flags are defined
        idx = [uncommented.index(j) for i in parsed_flags.values() for j in i]
        idx.append(len(uncommented))
    
        i = 0
        for keys, vals in parsed_flags.items():
    
            for prog in vals:
                parsed_flags[keys][prog] = unpack([i.split('=') for i in \
                    uncommented[idx[i]+1:idx[i+1]]])
                i += 1
    
        return parsed_flags


    def parse_mutations(self, txt):
        """Modifies the attributes of this object to make them passable to
        \"mutate()\".
        """
        def int_in_str(lst):
            converted = ''.join([item for item in lst if item.isdigit()])
            return int(converted)


#        One letter AA code to three letter code as presented on PyMOL Wiki.
        aa1 = list("ACDEFGHIKLMNPQRSTVWY")
        aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU \
            MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
        aa123 = dict(zip(aa1,aa3))

        mut_lst = list(open(txt, 'r'))
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
        """
        self.structures.extend(["{0}_{1}.pdb".format(
            self, i) for i in self.mut_df.axes[0]])

        for m in range(len(self.mut_df)):
            cmd.load(self.name+'/'+self.pdb_file)
            cmd.wizard('mutagenesis')
            name = "%s_%s" % (self, self.mut_df.axes[0][m])
            cmd.get_wizard().do_select('///%s/%s' %
                (self.mut_df["Chain"][m], str(self.mut_df["Residue"][m])))
            cmd.get_wizard().set_mode(self.mut_df["Mutation"][m])
            cmd.get_wizard().apply()
            cmd.save(name + ".pdb")
            cmd.reinitialize()
            makedir(name)
            shutil.move(name + ".pdb", name)
        
        
    def log(self, fname, proc_obj):
        """Used to log some of GROMACS output
        """
        log_file = open(fname, 'a')
        log_file.write(proc_obj.stdout.decode('utf-8'))
        log_file.close()


    def concoord(self):
        """Performs the CONCOORD procedure to generate protein structure
        ensembles. Takes additional flags from "flag_parse" as input (pass it as
        "flag_parse_output['CC']"). Make sure that \"CONCOORDRC.bash\" is
        sourced.
        """
        for s in self.structures:
            os.chdir(s[:-4])
            dist_input = [
                'dist',
                '-p', '%s' % s,
                '-op', '%s_dist.pdb' % s[:-4],
                '-og', '%s_dist.gro' % s[:-4],
                '-od', '%s_dist.dat' % s[:-4],
            ]
    
            disco_input = [
                'disco',
                '-d', '%s_dist.dat' % s[:-4],
                '-p', '%s_dist.pdb' % s[:-4],
                '-op', '',
                '-or', '%s_disco.rms' % s[:-4],
                '-of', '%s_disco_Bfac.pdb' % s[:-4]
            ]
            dist_input.extend(self.parsed_flags['CC']['DIST FLAGS'])
            disco_input.extend(self.parsed_flags['CC']['DISCO FLAGS'])
            subprocess.run(dist_input, input=b'1\n1')
            subprocess.run(disco_input)
            
            for n in range(1, self.n_structs+1):
                nr = str(n)
                makedir(nr)
                shutil.move(nr+'.pdb', nr)

            os.chdir('..')

        self.structures = [s[:-4]+'/'+str(nr)+'/'+str(nr)+'.pdb' \
            for s in self.structures for nr in range(1, self.n_structs+1)]

    
    def minimize(self, energy=False):
        """Standard GROMACS routine to minimize a molecule to avoid bad
        starting structures. If the argument \"energy\" is set to True, then
        \"gmx grompp" will use by default a .mdp file for single point energy
        evaluation without user supplied tables.
        This method will by default be run on all files of the structures generated by
        CONCOORD. The order in which to run those structures or adding/removing
        structures can be achieved by modifying the attribute \"structures\".
        """
        if energy:
            mdp = self.e_mdp

        else:
            mdp = self.min_mdp

        for s in self.structures:
            
            if "/" in s:
                os.chdir('/'.join(s.split('/')[:-1]))
                struct = s.split('/')[-1][:-4]

            else:
                struct = s[:-4]

            gmx_procedure = {
                'pdb2gmx_input': [
                    'gmx', '-quiet', 'pdb2gmx',
                    '-f', '%s.pdb' % struct,
                    '-o', '%s.gro' % struct,
                    '-p', '%s.top' % struct
            ],
                'editconf_input': [
                    'gmx', '-quiet', 'editconf',
                    '-f', '%s.gro' % struct,
                    '-o', '%s_box.pdb' % struct
            ],
                'grompp_input': [
                    'gmx', '-quiet', 'grompp',
                    '-f', mdp,
                    '-c', '%s_box.pdb' % struct,
                    '-p', '%s.top' % struct,
                    '-o', '%s.tpr' % struct
            ],
                'mdrun_input': [
                    'gmx', '-quiet', 'mdrun',
                    '-s', '%s.tpr' % struct,
                    '-deffnm', '%s_min' % struct,
                    '-nt', '1',
                    '-table', self.mdrun_table,
                    '-tablep', self.mdrun_table
            ]}

            proc_count = 0
            proc_ext = list(self.parsed_flags['gmx'].values())
            
            for proc in gmx_procedure.values():
                proc_input = proc + proc_ext[proc_count]
                subprocess.run(proc_input)
                proc_count += 1

            os.chdir(self.maindir)

if __name__ == '__main__':
    x = Protein("1pga.pdb", "mut.txt", "param2.txt")
    x.mutate()
    x.concoord()
    x.minimize()
