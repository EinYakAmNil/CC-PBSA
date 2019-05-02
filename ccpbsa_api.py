"""
This module contains the various functions used to create the workflow of
CC/PBSA. A method to compute mutational free energies.
Requirements (version numbers can probably be newer):
    - PyMOL 2.3
    - Python 3.7:
        - numpy
        - pandas
        - pymol (Easiest if you change your standard interpreter to PyMOL's)
    - CONCOORD 2.1.2
    - GROMACS 2019.1
"""
import numpy as np
import os
import shutil
import subprocess
import sys
import pandas as pd


def calc_entropy(directory, tpr_file=None):
    """
    Search all subfolders for minimized .gro structures (marked by \"_min\")
    and concatenate them into a trajectory file using \"gmx trjcat\".
    The output can then be analyzed by \"gmx covar\" to get the input for
    \"gmx anaeig\" which can calculate the entropy according to Schlitter's
    formula.
    \"gmx covar\" requires a Structure+mass file as input.
    If none is specified for this function, the first one found will be taken.
    """
    all_dirs = os.walk(directory)
    all_dirs = [i[0] for i in all_dirs] 
    start = os.getcwd()

    for dr in all_dirs:

        if tpr_file is None:
            os.chdir(dr)

            for f in os.listdir():

                if '.tpr' in f:
                    tpr_file = os.getcwd() + '/' + f
                    break

            os.chdir(start)

    files = []

    for dr in all_dirs:
        os.chdir(dr)

        for f in os.listdir():
            
            if '_min.trr' in f:
                check = subprocess.run([
                    'gmx', '-quiet', 'check',
                    '-f', os.getcwd() + '/' + f
                ], stderr=subprocess.PIPE).stderr.decode('utf-8').split()

                if not int(check[check.index('Step')+1]) == 0:
                    files.append(os.getcwd() + '/' + f)

        os.chdir(start)

    trjcat_input = [
        'gmx', '-quiet', 'trjcat', '-cat',
        '-f'
    ]
    trjcat_out = ['-o', '%s.xtc' % directory]
    trjcat_input.extend(files)
    trjcat_input.extend(trjcat_out)

    covar_input = [
        'gmx', '-quiet', 'covar',
        '-f', '%s.xtc' % directory,
        '-s', tpr_file,
        '-v', '%s_eigenvec.trr' % directory,
        '-fit', 'no'
    ]

    anaeig_input = [
        'gmx', '-quiet', 'anaeig',
        '-v', '%s_eigenvec.trr' % directory,
        '-s', tpr_file,
        '-entropy'
    ]

    os.chdir(directory)
    subprocess.run(trjcat_input, stdout=subprocess.PIPE)
    subprocess.run(covar_input, input=b'0', stdout=subprocess.PIPE)
    entropy = subprocess.run(anaeig_input, stdout=subprocess.PIPE)
    os.chdir('..')
    entropy = entropy.stdout.decode('utf-8')
    valstart = entropy.index('is')+2
    valend = entropy.index('J/mol')-1
    entropy = entropy[valstart:valend]

    return float(entropy)/1000 # J --> kJ


def extract_energies(min_gro, emdp_file, gmx_flags):
    """
    Performs a zero step energy minimization on the minimized .gro files at
    epsilon-r = 1 and infinite Cut-off. Then extract the Coulomb and Lennard
    Jones using "gmx energy".
    Returns two values:
        1. Coulomb Potential
        2. Lennard-Jones Potential
    """
    gro = [name for name in min_gro.split('/')][-1][:-4]

    energies = gmx_flags['ENERGIES'][0]

    gmxflags = gmx_flags.copy()
    gmxflags.pop('ENERGIES')

    pdb2gmx_input = [
        'gmx', '-quiet', 'pdb2gmx',
        '-f', min_gro,
        '-o', '%s_energy.gro' % gro,
        '-p', '%s_energy.top' %gro,
    ]

    editconf_input = [
        'gmx', '-quiet', 'editconf',
        '-f', min_gro,
        '-o', '%s_energy_box.pdb' % gro
    ]

    grompp_input = [
        'gmx', '-quiet', 'grompp',
        '-f', emdp_file,
        '-c', '%s_energy_box.pdb' % gro,
        '-p', '%s_energy.top' % gro,
        '-o', '%s_energy.tpr' % gro
    ]

    mdrun_input = [
        'gmx', '-quiet', 'mdrun',
        '-s', '%s_energy.tpr' % gro,
        '-deffnm', '%s_energy' % gro,
        '-nt', '1'
    ]

    energy_input = [
        'gmx', '-quiet', 'energy',
        '-f', '%s_energy.edr' % gro,
        '-o', '%s_energy.xvg' % gro
    ]

    gmx_procedure = [
        pdb2gmx_input,
        editconf_input,
        grompp_input,
        mdrun_input,
    ]

    i = 0
    for flags in gmxflags.values():
        gmx_procedure[i].extend(flags)
        subprocess.run(gmx_procedure[i])
        i += 1

#    if os.stat('%s_energy.edr' % gro).st_size == 0:
#        raise Exception("Could not calculate energy.")

    energy = subprocess.run(energy_input, input=energies, encoding='utf-8',
            stdout=subprocess.PIPE)

    energy = energy.stdout.split('\n')[-5:-1]
    energy = [i.split('  ') for i in energy]

    for i in range(len(energy)):
        energy[i] = [j for j in energy[i] if len(j) > 0]

    energy = [float(i[1]) for i in energy]
    lj = energy[0] + energy[2]
    coulomb = energy[1] + energy[3]

    return coulomb, lj


def sasa_mean(directory):
    """
    Calculates the mean solvent-accessible surface area of all minimized \".trr\"
    and \".tpr\" files (marked by \"_min_energy\") using \"gmx sasa\".
    The calculated area in \"*input directory*_area.xvg\" is the return value of
    this function.
    """
    all_dirs = os.walk(directory)
    all_dirs = [i[0] for i in all_dirs] 
    start = os.getcwd()
    tpr_files = []
    trr_files = []

    for dr in all_dirs:
        os.chdir(dr)

        for f in os.listdir():

            if '_min_energy.tpr' in f:
                tpr_files.append(os.getcwd() + '/' + f)

            if '_min.trr' in f:
                check = subprocess.run([
                    'gmx', '-quiet', 'check',
                    '-f', os.getcwd() + '/' + f
                ], stderr=subprocess.PIPE).stderr.decode('utf-8').split()

                if not int(check[check.index('Step')+1]) == 0:
                    trr_files.append(os.getcwd() + '/' + f)

        os.chdir(start)

    assert len(tpr_files) == len(trr_files), \
        "Something might have went wrong during energy minimization. The \
        amount of .tpr and .trr files are not equal."

    sasa_input = [
        'gmx', '-quiet', 'sasa',
        '-surface', '0',
        '-pbc', 'no'
    ]

    os.chdir(directory)
    for f in range(len(tpr_files)):
        sasa_input.extend([
            '-f', trr_files[f],
            '-s', tpr_files[f],
            '-o', '%d_area.xvg' % f
        ])
        subprocess.run(
            sasa_input,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        sasa_input = sasa_input[:-6]

    area_files = [f for f in os.listdir() if '_area.xvg' in f]
    sasas = []

    for f in area_files:
        time_area = list(open(f, 'r'))[-1]
        area = time_area.split()[-1]
        sasas.append(float(area))

    mean_sasa = np.array(sasas)
    os.chdir('..')

    return np.mean(mean_sasa)


def concoord(pdb_file, flags):
    """
    Performs the CONCOORD procedure to generate protein structure ensembles.
    Takes additional flags from "flag_parse" as input (pass it as
    "flag_parse_output['CC']").
    Returns the name of the structure ensemble.
    """
    pdb = pdb_file[:-4]

    dist_input = [
        'dist',
        '-p', '%s' % pdb_file,
        '-op', 'dist_%s' % pdb_file,
        '-og', 'dist_%s.gro' % pdb,
        '-od', 'dist_%s.dat' % pdb,
    ]

    disco_input = [
        'disco',
        '-d', 'dist_%s.dat' % pdb,
        '-p', 'dist_%s' % pdb_file,
#        '-op', '',
        '-on', 'disco_%s' % pdb,
        '-or', 'disco_%s.rms' % pdb,
        '-of', 'disco_Bfac_%s' %pdb_file
    ]
    
    dist_input.extend(flags['DIST FLAGS'])
    disco_input.extend(flags['DISCO FLAGS'])
    n = int(disco_input[disco_input.index('-n')+1])

    subprocess.run(dist_input, input=b'1\n1')
    subprocess.run(disco_input)
    
#    for f in range(1, n+1):
#        os.mkdir(str(f))
#        shutil.move("%s.pdb" % f, str(f))
    return 'disco_%s' % pdb_file


def manage_CC_files(pdb_ensemble):
    """
    Splits the protein structure ensemble of the CONCOORD program into
    individual enumerated .pdb structure files.
    If the structure ensembles are not split into seperate files, pdb2gmx will
    only convert the first structure in the file and stop.
    Returns the number of structures it counted.
    """
#    Use "gmx trjconv" to split the structure ensemble.
    trjconv_input = [
        'gmx', '-quiet', 'trjconv',
        '-f', pdb_ensemble,
        '-s', pdb_ensemble,
        '-o', '.pdb',
        '-sep'
    ]

    subprocess.run(trjconv_input, input=b'0')

#    Count the number of split structures and puts them into their directory.
    n_structures = os.listdir()
    n_structures = [pdb for pdb in n_structures if pdb[-4:] == '.pdb']
    n_structures = [int(num[:-4]) for num in n_structures if num[:-4].isdigit()]
    n_structures = max(n_structures)+1

    for new_dir in range(n_structures):
        os.mkdir(str(new_dir))
        shutil.move('%d.pdb' % new_dir, str(new_dir))

    return n_structures

def flag_parse(filename):
    """
    Takes the .txt file in which the flags of CC/PBSA are saved in as input.
    NOT viable (already defined) flags for "emin" are:
    - CONCOORD:
        - dist: -p, -op, -og, -od
        - disco: -d, -p, -on, -or, -of
    - Gromacs: -quiet
        - pdb2gmx: -f, -p, -o
        - editconf: -f, -o
        - grompp: -f, -c, -p, -o
        - mdrun: -s, -deffnm, -table, -tablep, -nt
        - energy: -f, -o

    Contains a subfunction "file_parse" which turns the lines of the flag file into a list of strings.
    Newlines (\"\\n\") and comments (";") are removed.
    """

    def file_parse(f_raw):
        parsed = []

        for line in f_raw:
            line = line.split('\n')
            line = [j for j in line if j is not '']
            line = ''.join(line)

            if ';' in line:
                line = line[:line.index(';')]
            
            if len(line) > 0:
                parsed.append(' '.join(line.split()))

        return parsed


    file_raw = list(open(filename))
    file_parsed = file_parse(file_raw)
    flags = {
        'CC': {
            'DIST FLAGS': [],
            'DISCO FLAGS': []
        },

        'gmx': {
            'PDB2GMX FLAGS': [],
            'EDITCONF FLAGS': [],
            'GROMPP FLAGS': [],
            'MDRUN FLAGS': [],
            'ENERGIES': []
        }
    }

    indices = [file_parsed.index(k) for val in flags.values() for k in val.keys()]
    indices.append(len(file_parsed))
    idx = 0

    for key, val in flags.items():
        
        for k in val.keys():
            flags[key][k] = [
                j for i in \
                file_parsed[indices[idx]+1:indices[idx+1]] \
                for j in i.split('=')
            ]
            idx += 1

    return flags


def interpret_flags(sys_argvs):
    """
    Interprets the list of sys.argv as the flag input required for CC/PBSA:
    -p: .pdb file of the protein
    -m: .txt file for the list of mutations
    -d: .mdp file used by GROMACS.
    -f: .txt file defining additional flags for CONCOORD and GROMACS.
    Returns a dictioary with the file name (value) to the corresponding flag
    (key).
    """
    if '-h' in sys_argvs:
        help(sys.modules[__name__])

    flags = {
        '-p': '.pdb',
        '-m': '.txt',
        '-d': '.mdp',
        '-f': '.txt'
    }

    for key, val in flags.items():
        file_input = sys_argvs[sys_argvs.index(key)+1]
        
        if key in sys_argvs and val == file_input[-4:]:
            flags[key] = file_input

        else:
            raise RuntimeError("%s requires a %s file" % (key, val))

    return flags


def file_find(filename):
    """
    This function will try to find a file matching the input.
    It searches from the current directory up to all parent directories.
    Closer results are prioritized.
    """
    floc = False
    cwd = os.getcwd()
    parents = [i for i in cwd.split('/')][1:]
    parents[0] = '/' + parents[0]

    for i in range(1, len(parents)):
        parents[i] = parents[i-1] + '/' + parents[i]

#    Invert list to start search from current directory.
    parents = parents[::-1]

    for i in parents:

        if filename in os.listdir(i):
            floc = i + '/' + filename
            return floc

    if floc is False:
        raise FileNotFoundError('Could not find file %s' % filename)


def parse_input(mutations_txt):
    """
    Takes a file with mutation descriptions in the form of (without the spaces):
    *(Chain_)* *Original Amino Acid* *Residue Number* *New Amino Acid*
    Each mutation must be one line.
    Returns a dictionary with the mutation description as the key and the
    information contained in the file split up in that order.
    """
#    Convert one letter amino acid code to three letter code
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    aa123 = dict(zip(aa1,aa3))
    
    mut_lst = open(mutations_txt, 'r')
#    Remove newlines from mutation file.
    mut_lst = [str(i[:-1]) for i in mut_lst]
    mut_loc = {}

    for i in mut_lst:
#        Detect if a chain is specified
        chain = None
        if '_' in i:
            chain = i[0]
            aa = i[2]

        else:
            aa = i[0]

#        Detect residue number
        n = [int(s) for s in i if s.isdigit()] 
        resi = 0
        for j in range(len(n)):
            resi += n[j] * 10 ** (len(n)-j-1)

        mut = aa123[i[-1]]

        mut_loc[i] = [chain, aa, resi, mut]

#    Output is a dictionary with the original line in mutations_txt as key
    return mut_loc


def make_gxg(aa):
    """
    Use PyMOL to generate GXG tripeptides for further energy evaluation of the
    unfolded state of proteins in the thermodynamic cycle.
    """
    import pymol
    pymol.finish_launching(['pymol', '-qc'])
    cmd = pymol.cmd

    cmd.fab("G%sG" % aa, "G%sG" % aa)
    cmd.remove('hydrogens')
    cmd.save("G%sG.pdb" % aa)
    cmd.reinitialize()


#NOTE:
#Residue and chain selection must be correct since PyMOL does not
#raise Errors if selection is not valid.
#Otherwise there might be no mutation in the output files.
def mutate(pdb_file, residue, mutation, save_as, chain=None):
    """
    Takes a .pdb file and mutates it using PyMOL.
    Additional arguments are the residue number, the mutation,
    and how the mutated protein should be named.
    Sometimes the chain must be specified to mutate correctly.
    """
#    quietly import PyMOL
    import pymol
    pymol.finish_launching(['pymol', '-qc'])
    cmd = pymol.cmd

    prot_name = [i for i in pdb_file.split('/')][-1][:-4]
    save_as = prot_name + '_' + save_as

    pymol.cmd.load(pdb_file)
    cmd.wizard("mutagenesis")

    try:
#        Mutation procedure using PyMOL's Wizard module
        if chain:
            cmd.get_wizard().do_select("///%s/%s" % (chain, residue))
    
        else:
            cmd.get_wizard().do_select("////%s" % (residue))
        
        cmd.get_wizard().set_mode(mutation)
        cmd.get_wizard().apply()
        print("mutated %s successfully" % save_as)
#        Save Mutation as "(original filename)_(change as written in mutation file).pdb"
        cmd.save("%s.pdb" % save_as)
    
    except pymol.CmdException:
        print("\nfailed to mutate %s" % save_as)

    print('\n')
    cmd.reinitialize()


def emin(pdb_file, mdp_file, gmx_flags, table=None):
    """
    Minimizes the energy of a structure according to the gromacs standard procedure
    Input arguments are:
    -the .pdb file to be minimized
    -the .mdp file used for grompp
    -additional parameters like forcefield or box type.
     Parameters affecting naming input and output are already taken care of.
    -User defined potential can also be supplied if specified in the .mdp file.
     Then the location of the .xvg file which mdrun should use must also be given.
    """
    assert type(gmx_flags) is dict, \
        "Pass the flags for gmx programs as dictionary."
    pdb_name = [i for i in pdb_file.split('/')][-1][:-4]

#    Function will raise a KeyError since second argument of pop is not defined.
#    This is a sort of control to check whether or not the input is correct.
    gmxflags = gmx_flags.copy()
    gmxflags.pop('ENERGIES')

#    Linkage of in-/output names to ensure flow of the gmx procedure.
    pdb2gmx_input = [
        'gmx', '-quiet', 'pdb2gmx',
        '-f', pdb_file,
        '-o', '%s.gro' % pdb_name,
        '-p', '%s_topol.top' % pdb_name
    ]

    editconf_input = [
        'gmx', '-quiet', 'editconf',
        '-f', '%s.gro' % pdb_name,
        '-o', '%s_box.pdb' % pdb_name
    ]

    grompp_input = [
        'gmx', '-quiet', 'grompp',
        '-f', mdp_file,
        '-c', '%s_box.pdb' % pdb_name,
        '-o', '%s_min.tpr' % pdb_name,
        '-p', '%s_topol.top' % pdb_name
    ]

    if table is not None:
        mdrun_input = [
            'gmx', '-quiet', 'mdrun',
            '-table', table,
            '-tablep', table,
            '-deffnm', '%s_min' % pdb_name,
            '-s', '%s_min.tpr' % pdb_name,
            '-nt', '1'
        ]

    else:
        mdrun_input = [
            'gmx', '-quiet', 'mdrun',
            '-deffnm', '%s_min' % pdb_name,
            '-s', '%s_min.tpr' % pdb_name,
            '-nt', '1'
        ]
 
#    Still need to find a good way to add all the addtitional flags to the
#    procedure.
#    Right now it the function is unflexible to changes in the output of 'flag'.
    gmx_procedure = [
        pdb2gmx_input,
        editconf_input,
        grompp_input,
        mdrun_input
    ]

    i = 0
    for flags in gmxflags.values():
        gmx_procedure[i].extend(flags)
        x = subprocess.run(gmx_procedure[i])
        i += 1

#    Finally we can pass the parameters to their respective function.
    return os.stat('%s_min.edr' % pdb_name).st_size

#If this file is called as a script instead of import.
if __name__ == '__main__':
    run_files = interpret_flags(sys.argv)
    pdb_file = file_find(run_files['-p'])
    prot_name = pdb_file[:-4]
    cc_pbsa = prot_name
    mut_file = file_find(run_files['-m'])
    mdp = file_find(run_files['-d'])
    energy_mdp = file_find('energy.mdp')
    flags = file_find(run_files['-f'])
    tab = file_find('table4r-6-12.xvg')
    cc_flags = flags['CC']
    gmx_flags = flags['gmx']

    if os.path.isdir(prot_name):
        save_in_dir_ans = input("\"%s\" directory already exists.\n\
            Should all data be stored in there? (y/n)")

        while save_in_dir_ans != 'y' or save_in_dir_ans != 'n':

            if save_in_dir_ans == 'y':
                break

            elif save_in_dir_ans == 'n':
                mkdir_ans = input("Enter name of the to create directory:\n")
                os.mkdir(mkdir_ans)
                cc_pbsa = mkdir_ans

            else:
                save_in_dir_ans = input("\"%s\" directory already exists.\n\
                    Should all data be stored in there? (y/n)")

    else:
        os.mkdir(prot_name)
        print("Made directory %s" % prot_name)

    os.chdir(cc_pbsa)

#    Make directory for unmodified protein.
    os.mkdir(prot_name)
    os.chdir(prot_name)
    cc_name = concoord(pdb_file, cc_flags)
    n_struct = manage_CC_files(cc_name)
    all_data = pd.DataFrame(columns=list(range(n_struct)))

    G_wt_coulomb_mean = 0
    G_wt_lj_mean = 0

    for n in range(n_struct):
        os.chdir(str(n))
        emin('%d.pdb' % n, mdp, gmx_flags, table=tab)
        G_wt_coulomb, G_wt_lj = \
            extract_energies('%d_min.gro' % n, energy_mdp, gmx_flags)
        G_wt_coulomb_mean += G_wt_coulomb
        G_wt_lj_mean += G_wt_lj
        os.chdir('..')

    G_wt_coulomb_mean /= n_struct
    G_wt_lj_mean /= n_struct
    os.chdir('..')
    G_wt_entropy = calc_entropy(prot_name)
    G_wt_sasa = sasa_mean(prot_name)
    
    mutations = parse_input(mut_file)
    n_muts = len(mutations.keys())
    G_mut_coulomb = dict(zip(mutations.keys(), np.zeros(n_muts)))
    G_mut_lj = dict(zip(mutations.keys(), np.zeros(n_muts)))
    G_mut_S = dict(zip(mutations.keys(), np.zeros(n_muts)))
    G_mut_sasa = dict(zip(mutations.keys(), np.zeros(n_muts)))

    for mut_name in mutations:
        os.mkdir(mut_name)
        os.chdir(mut_name)
        
        chain = mutations[mut_name][0]
        resi = mutations[mut_name][2]
        mut = mutations[mut_name][3]

        mutate(pdb_file, resi, mut, mut_name, chain)
        cc_name = concoord('%s_%s.pdb' % (prot_name, mut_name), cc_flags)
        n_struct = manage_CC_files(cc_name)

        for n in range(n_struct):
            os.chdir(str(n))
            emin('%d.pdb' % n, mdp, gmx_flags, tab)
            coulomb, lj = \
                extract_energies('%d_min.gro' % n, energy_mdp, gmx_flags)
            G_mut_coulomb[mut_name] += coulomb
            G_mut_lj[mut_name] += lj
            os.chdir('..')

        G_mut_coulomb[mut_name] /= n
        G_mut_lj[mut_name] /= n
        os.chdir('..')
        G_mut_S[mut_name] = calc_entropy(mut_name)
        G_mut_sasa[mut_name] = sasa_mean(mut_name)

    ddG_stability = pd.DataFrame(
        columns=['CALC', 'SOLV', 'COUL', 'LJ', 'SAS', '-TS'],
        index=mutations
    )

    alpha = 0.224
    beta = 0.217
    gamma = 0.0166*4.184*100 # kcal mol^-1 Å^-2 --> kJ mol^-1 nm^-2
    tau = 0.0287
    T = 298.15 # K
    gxg_coulomb = pd.read_csv('GXG_Coulomb.csv', index_col=0)
    gxg_LJ = pd.read_csv('GXG_LJ.csv', index_col=0)
    gxg_sasa = pd.read_csv('GXG_Area.csv', index_col=0)
    gxg_S = pd.read_csv('GXG_S.csv', index_col=0)

    for mut_name in mutations:
        aa1 = mutations[mut_name][1]
        aa2 = mut_name[-3]
        ddG_fold_stability[mut_name]['COUL'] = \
            alpha * (G_wt_coulomb_mean - G_mut_coulomb[mut_name] - \
            gxg_coulomb[aa2][aa1])
        ddG_fold_stability[mut_name]['LJ'] = \
            beta * (G_wt_lj_mean - G_mut_lj[mut_name] - \
            gxg_lj[aa2][aa1])
        ddG_fold_stability[mut_name]['SAS'] = \
            gamma * (G_wt_sasa_mean - G_mut_sasa[mut_name] - \
            gxg_sasa[aa2][aa1])
        ddG_fold_stability[mut_name]['-TS'] = \
            tau * T * (G_wt_S_mean - G_mut_S[mut_name] - \
            gxg_S[aa2][aa1])

    ddG_fold_stability.to_csv('%s.csv' % prot_name)
