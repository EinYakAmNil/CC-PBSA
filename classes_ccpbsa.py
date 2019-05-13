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
        self.energy_df = pd.DataFrame(0,
            columns=['SOLV', 'COUL', 'LJ', 'SAS', 'S'],
            index=[self.name])
        self.energy_df.index.name = 'MUTATION'
#        self.cpt
        
        makedir(self.name)
        os.chdir(self.name)
        makedir(self.name)
        shutil.copy("../" + self.pdb_file, self.name)
        self.maindir = os.getcwd()
        self.subdir = []
        self.trr = []
        

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
        The format of the mutation instruction should be (in one letter code
        and without spaces):
            (*Chain_*) *OriginalAA* *ResidueNumber* *NewAA*
        for exapmle:
            - A20G (for a monomer)
            - B_H10I (for a dimer with a chain named \"B\")
        """
        self.structures.extend([i+".pdb" for i in self.mut_df.axes[0]])
        new_idx = list(self.energy_df.index)
        new_idx.extend([m for m in self.mut_df.axes[0]])
        self.energy_df = self.energy_df.reindex(new_idx, fill_value=0)

        for m in range(len(self.mut_df)):
            cmd.load(self.name+'/'+self.pdb_file)
            cmd.wizard('mutagenesis')
            mut_key = self.mut_df.axes[0][m]
#            name = "%s_%s" % (self, mut_key)
            cmd.get_wizard().do_select('///%s/%s' %
                (self.mut_df["Chain"][m], str(self.mut_df["Residue"][m])))
            cmd.get_wizard().set_mode(self.mut_df["Mutation"][m])
            cmd.get_wizard().apply()
            cmd.save(mut_key + ".pdb")
            cmd.reinitialize()
            makedir(mut_key)
            shutil.move(mut_key + ".pdb", mut_key)
        
        
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
                self.subdir.append(os.getcwd() + '/' + nr)

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
        Updates the \"structures\" attribute to the minimized structure files.
        """
        if energy:
            mdp = self.e_mdp
            deffnm = "energy"

        else:
            mdp = self.min_mdp
            deffnm = "min"

        for s in self.structures:
            
            if "/" in s:
                os.chdir('/'.join(s.split('/')[:-1]))
                struct = s.split('/')[-1][:-4]

            else:
                struct = s[:-4]

            gmx_procedure = {
                'pdb2gmx_input': [
                    'gmx', '-quiet', 'pdb2gmx',
                    '-f', struct+s[-4:],
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
                    '-deffnm', '%s_%s' % (struct, deffnm),
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

        if energy:
            self.structures = [s[:-4]+"_energy.gro" for s in self.structures]

        else:
            self.structures = [s[:-4]+"_min.gro" for s in self.structures]

        self.trr = [d+'/'+trr for d in self.subdir for trr in os.listdir(d) \
                if '.trr' == trr[-4:]]


    def lj(self):
        """Calculate the Lennard-Jones Potential
        """
        self.minimize(energy=True)

        def tryfloat(anything):
            try:
                return float(anything)
            except ValueError:
                return 0
        
        for s in self.structures:

            if "/" in s:
                splt = s.split('/')
                df_idx = splt[-3]
                os.chdir('/'.join(splt[:-1]))
                struct = splt[-1][:-4]

            else:
                struct = s[:-4]
                df_idx = s[:-4]

            gmx_energy = ['gmx', '-quiet', 'energy', '-f', '%s.edr' % struct]
            lj_proc = subprocess.run(gmx_energy, input=b'5\n7',
                stdout=subprocess.PIPE)
            self.log("lj.txt", lj_proc)

#            The last 2 lines contain the energies. Find the values, convert
#            them to float and add them to the total amount.
            lje = lj_proc.stdout.decode('utf-8').split('\n')[-2:]
            lje = [j for i in lje for j in i.split('  ')]
            lje = np.array(list(map(tryfloat, lje))).sum()
            self.energy_df['LJ'][df_idx] += lje

            os.chdir(self.maindir)

#        Calculate the mean.
        self.energy_df['LJ'] /= self.n_structs
        self.energy_df.to_csv("Values.csv")


    def electrostatics(self):
        """Calculate the solvation free energy and Coulomb energy using a
        programm from Franziska Kranz
        """
        pass
        self.energy_df['SOLV'] /= self.n_structs
        self.energy_df['COUL'] /= self.n_structs
        self.energy_df.to_csv("Values.csv")


    def sasa(self):
        """Calculate the solvent accessible surface area using \"gmx sasa\"
        """
        sasa_input = ['gmx', '-quiet', 'sasa',
            '-surface', 'system',
            '-pbc', 'no',
            '-f', '-s']

        for s in self.structures:

            if "/" in s:
                splt = s.split('/')
                df_idx = splt[-3]
                os.chdir('/'.join(splt[:-1]))
                struct = splt[-1]

            else:
                struct = s
                df_idx = s[:-4]

            s_input = sasa_input[:-1] + [struct] + [sasa_input[-1]] + [struct]
            subprocess.run(s_input)
            area = list(open("area.xvg"))[-1]
            area = area.split()[-1]
            self.energy_df['SAS'][df_idx] += float(area)

            os.chdir(self.maindir)

        self.energy_df['SAS'] /= self.n_structs
        self.energy_df.to_csv("Values.csv")


    def schlitter(self):
        """Calculate a upper limit of the entropy according to Schlitter's
        formula. Uses the .trr files it can find in the respective directory in
        the attribute \"subdir\".
        """
        trjcat = ['gmx', '-quiet', 'trjcat', '-cat', 'yes', '-f']
        trrs = ['{0}/{0}_min.trr'.format(i) for i in range(1, self.n_structs+1)]
        trjcat += trrs
        covar = ['gmx', '-quiet', 'covar', '-f', 'trajout.xtc',
            '-s', '1/1.gro', '-fit',
            'no', '-pbc', 'no']
        anaeig = ['gmx', '-quiet', 'anaeig', '-v', 'eigenvec.trr', '-entropy']
        for d in self.energy_df.axes[0]:
            os.chdir(d)
            subprocess.run(trjcat)
            subprocess.run(covar, input=b'0')
            entropy = subprocess.run(anaeig, stdout=subprocess.PIPE)
            entropy = entropy.stdout.decode('utf-8')
            valstart = entropy.index('is ')+3
            valend = entropy.index(' J/mol')
            entropy = entropy[valstart:valend]
            self.energy_df['S'][d] += float(entropy)
            os.chdir('..')

        self.energy_df.to_csv("Values.csv")


    def fullrun(self):
        """Calculate the total mutational free energy according to the CC/PBSA
        method. This function will suffice fully, if no other behaviour is
        required
        """
        self.mutate()
        self.concoord()
        self.minimize()
        self.schlitter()
        self.sasa()
        self.electrostatics()
        self.lj()
        print(self.energy_df)


class GXG(Protein):
    """Subclassed from Protein. Does not need specification of input
    structures. Creates four .csv look-up-tables for the Thermodynamic Cycle.
    It is recommended to use the same parameters for this calculation as for
    the Protein.
    """
    def __init__(self, flags,
        min_mdp="/home/linkai/CC_PBSA/min.mdp",
        energy_mdp="/home/linkai/CC_PBSA/energy.mdp",
        mdrun_table="/home/linkai/CC_PBSA/table4r-6-12.xvg"):
        """Only the parameters need to be specified. Structures are the 20 GXG
        amino acids.
        """
        self.parsed_flags = self.parse_flags(flags)
        self.n_structs = int(self.parsed_flags['CC']['DISCO FLAGS'] \
            [self.parsed_flags['CC']['DISCO FLAGS'].index('-n')+1])
        self.e_mdp = energy_mdp
        self.min_mdp = min_mdp
        self.mdrun_table = mdrun_table
        self.structures = []
        self.subdir = []
        
        makedir('GXG')
        os.chdir('GXG')
        self.maindir = os.getcwd()


    def makepeps(self):
        """Create the tripeptides using PyMOL. Each peptide is saved in a
        separate directory similar to mutations in Protein.
        """
        aa1 = list("ACDEFGHIKLMNPQRSTVWY")

        for x in aa1:
            gxg = "G%sG" % x
            cmd.fab(gxg, gxg)
            cmd.save(gxg+'.pdb')
            cmd.reinitialize()
            makedir(gxg)
            shutil.move(gxg+'.pdb', gxg)
            self.structures.append(gxg+'.pdb')


if __name__ == '__main__':
    x = Protein("1pga.pdb", "mut.txt", "param.txt")
    print(x.__dict__)
#    x.fullrun()
