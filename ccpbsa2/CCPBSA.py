import os
import pandas as pd
import shutil
import subprocess as sub
import sys
import pymol
pymol.finish_launching(['pymol', '-Qc'])
cmd = pymol.cmd

def flatten_list(lst=[]):
    """Reduce arbitrarily deeply nested lists to a unnested list and returns
    it.
    """
    unpacked = []

    for i in list(lst):

        if type(i) is list:
            unpacked.extend(flatten_list(list(i)))

        else:
            unpacked.append(i)

    return unpacked
    

def parse_mutations(*raws):
    """Parse strings of mutations.
    """
    MutFrame = pd.DataFrame(
        [],
        index=raws,
        columns=['chain', 'aa old', 'res number', 'aa new']
    )

    for i in raws:
        j = i.replace(" ", "")
        j = j.split(",")
        MutFrame.loc[i, 'chain'] = [k[0] if '_' in k else '' for k in j]
        MutFrame.loc[i, 'aa old'] = [k[2] if '_' in k else k[0] for k in j]
        MutFrame.loc[i, 'res number'] = \
            [int(k[3:-1]) if '_' in k else int(k[1:-1]) for k in j]
        MutFrame.loc[i, 'aa new'] = [k[-1] for k in j]

    return MutFrame


def mutate(pdb, mut_db):
    """Mutates a .pdb file according to the mutation(s) provided. The mutation
    database should be a pandas.DataFrame object. The mutated protein will be
    saved to a separate file named after the mutation performed.
    """
#    From the PyMOL Wiki (https://pymolwiki.org/index.php/Aa_codes)
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    aa123 = dict(zip(aa1,aa3))

    for i in mut_db.index:
        cmd.load(pdb)
        cmd.wizard('mutagenesis')

        for j in range(len(mut_db.loc[i, 'aa new'])):
            cmd.get_wizard().do_select('///{}/{}'.format(
                mut_db.loc[i, 'chain'][j], mut_db.loc[i, 'res number'][j]
            ))

            if mut_db.loc[i, 'aa new'][j] in aa1:
                cmd.get_wizard().set_mode(aa123[mut_db.loc[i, 'aa new'][j]])

            elif mut_db.loc[i, 'aa new'][j] in aa1:
                cmd.get_wizard().set_mode(mut_db.loc[i, 'aa new'][j])

            else:
                raise SyntaxError("Mutation notation can't be parsed")

            cmd.get_wizard().apply()
        
        cmd.save(i + '.pdb')
        cmd.reinitialize()


def gmx(prog, **subkwargs):
    """Tries to execute the GROMACS program specified in the first argument.
    Keyword arguments are passed to the subprocess.run function.
    """
    while True:

        try:
           gmx_sub = sub.run(['gmx', '-quiet'] + prog, **subkwargs)

        except PermissionError:
            sleep(120)
            continue

        break

    if gmx_sub.returncode == 0:
        return gmx_sub

    else:
        print(gmx_sub.stderr)
        raise Exception("Something in GROMACS went wrong!")


def parse_parameters(paramsfile):
    """Additional flags and interactive input for each subprogram called by
    CC/PBSA can be specified in here. Returns a dictionary wtih the proper
    flags parsed and a dictionary for interactive inputs.
    """
    with open(paramsfile, 'r') as raw:
#        Remove newlines and empty lines
        params = [i[:-1] for i in raw.readlines()]
        params = [i for i in params if len(i) > 0]

#    Determine programs
    prog_idx = [
        i for i in range(len(params)) if "[" in params[i] and "]" in params[i]
    ] + [len(params)]
    params = [params[prog_idx[j]:prog_idx[j+1]] for j in range(len(prog_idx)-1)]

    inputs = {
        i[0][1:-1]: "\n".join([
            j.replace("<<<", "").strip() for j in i if "<<<" in j
        ]) for i in params
    }

    flags = {
        i[0][1:-1]: flatten_list([
            j.strip().split() for j in i \
            if "<<<" not in j and "[" not in j
        ]) for i in params
    }

    return inputs, flags


def concoord(struct_file, flags, inputs, **subkwargs):
    """Create structure ensembles using the CONCOORD algorithm. Accepts .pdb
    and .gro files.
    """
    if struct_file[-3:] == "pdb":
        dist = sub.run(
            ['dist', '-p', struct_file] + flags['dist'],
            input=inputs['dist'],
            **subkwargs,
        )

    elif struct_file[-3:] == "gro":
        dist = sub.run(
            ['dist', '-g', struct_file] + flags['dist'],
            input=inputs['dist'],
            **subkwargs,
        )

    else:
        raise FileNotFoundError("dist.exe won't be able to determine file type.")

    disco = sub.run(
        ['disco', "-op", ""] + flags['disco'],
        **subkwargs
    )

    return dist, disco


def energy_minimization(struct_file, flags, inputs, **subkwargs):
    """Performs energy minimization on a .pdb or .gro file with the given
    flags. Makes Use of the gmx function and returns all the used subprocesses.
    """
    if struct_file[-3:] == 'pdb':
        pdb2gmx = gmx(
            ['pdb2gmx', '-f', struct_file] + flags['pdb2gmx'],
            input=inputs['pdb2gmx'],
            **subkwargs
        )
        editconf = gmx(['editconf'] + flags['editconf'], **subkwargs)

    else:
        pdb2gmx = None
        editconf = None

    grompp = gmx(['grompp'] + flags['grompp'], **subkwargs)
    mdrun = gmx(['mdrun'] + flags['mdrun'], **subkwargs)

    return pdb2gmx, editconf, grompp, mdrun


def log(process, logname):
    """Write stdout and stderr of a subprocess object to a file.
    """
    with open(logname, 'a') as logfile:
        
        if hasattr(process.stdout, 'decode'):
            logfile.write(process.stdout.decode(sys.stdin.encoding))
            logfile.write(process.stderr.decode(sys.stdin.encoding))

        else:
            logfile.write(process.stdout)
            logfile.write(process.stderr)


def parse_xvg(xvg):
    """Takes an energy.xvg file generated by \"gmx energy\" and reads out the
    latest energy values. Returns a dictionary containing the assigned values.
    If the \"-sum\" flag used to generate the file, then all energies will be
    combined into one key.
    """
    with open(xvg, 'r') as raw:
        lines = raw.readlines()
        energies = [i.split('\"')[1] for i in lines if "@ s" in i]

        if "Sum" in energies:
            energies.pop(energies.index("Sum"))
            energies = [" + ".join(energies)]

        values = [float(i) for i in lines[-1].split()[1:]]

    return dict(zip(energies, values))


def gropbe(params, nchains, logname="gropbe.log", **subkwargs):
    """Uses GroPBE to calculate the Coulomb and Solvation energy of a .tpr
    file. The number of chains in the .tpr file needs to be specified. Logs the
    output to \"gropbe.log\" or any other name is not specifically turned off
    (set as None). Returns Coulomb and Solvation energy as floats.
    """
    chains = ",".join(str(i) for i in range(nchains))
    pbe = sub.run(
        ['gropbe', params],
        input=chains,
        **subkwargs
    )

    if logname is not None:
        log(pbe, logname)

    coul = float(next(
        i.split()[-2] for i in pbe.stdout.split("\n") if "Coulomb" in i
    ))
    solv = float(next(
        i.split()[-2] for i in pbe.stdout.split("\n") if "Solvation" in i
    ))

    return coul, solv


def schlitter_entropy(
    structs, flags, inputs, logname="entropy.log", **subkwargs
):
    """Uses gmx trjcat, covar and anaeig to return the entropy according to
    Schlitter's and the quasiharmonic formula.
    """
    trjcat = gmx(['trjcat', '-f'] + structs + flags['trjcat'], **subkwargs)
    covar = gmx(['covar'] + flags['covar'], input=inputs['covar'], **subkwargs)
    anaeig = gmx(['anaeig'] + flags['anaeig'], **subkwargs)

    if logname is not None:
        log(trjcat, logname)
        log(covar, logname)
        log(anaeig, logname)
    
    schlitter = float(next(
        i.split()[-3] for i in anaeig.stdout.split('\n') if 'Schlitter' in i
    ))
    quasi = float(next(
        i.split()[-3] for i in anaeig.stdout.split('\n') if 'Quasiharmonic' in i
    ))

    return schlitter, quasi


if __name__ == '__main__':
    traj = 'traj.trr'
    subkwargs = {
        'encoding': sys.stdin.encoding,
        'capture_output': True
    }
    inputs, flags = parse_parameters("params")
    
    with open("muts", 'r') as raw:
        mut_df = parse_mutations(*[i[:-1] for i in raw.readlines() if len(i) > 1])
        mutate("1ayi.pdb", mut_df)

    if "test" not in os.listdir():
        os.mkdir("test")
        shutil.copy("1ayi.pdb", "test")
        os.chdir("test")
        energy_minimization("1ayi.pdb", flags, inputs, **subkwargs)
        concoord("confout.gro", flags, inputs, **subkwargs)

        for i in range(3):
            pdb = "{}.pdb".format(i+1)
            d = str(i+1)
            os.mkdir(d)
            shutil.move(pdb, d)
            os.chdir(d)
            energy_minimization(pdb, flags, inputs, **subkwargs)
            os.chdir('..')

        structs = [i+"/"+ traj for i, j, k in os.walk('.') if traj in k][1:]
        schlitter, quasi = schlitter_entropy(
            structs, flags, inputs, **subkwargs
        )
        print(schlitter, quasi)

#        structs = [i+"/"+ traj for i, j, k in os.walk('.') if traj in k][1:]

#        energykwargs = {**subkwargs, **{'input': "5\n7"}}
#        gmx(["energy"], **energykwargs)
#        gmx(["sasa", "-f", "confout.gro"], **subkwargs)

    else:
        os.chdir("test")

    
