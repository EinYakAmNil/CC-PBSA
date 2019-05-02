import ccpbsa_api as api
import os
import numpy as np
import sys
import pandas as pd # pandas indexes DataFrames differently from numpy!
import shutil

mdp = api.file_find('min.mdp')
energy_mdp = api.file_find('energy.mdp')
flags = api.file_find('param2.txt')
flags = api.flag_parse(flags)
tab = api.file_find('table4r-6-12.xvg')
mutations = api.file_find('mutations_1pga.txt')
mutations = api.parse_input(mutations)

gmx_flags = flags['gmx']

os.chdir('1pga')
omitted = open("omitted.log", 'w')

ddG_fold_stability = pd.DataFrame(
    columns=['CALC', 'SOLV', 'COUL', 'LJ', 'SAS', '-TS'],
    index=mutations
)

val_index = ['1pga']
val_index.extend(mutations)
values = pd.DataFrame(columns=['COUL', 'LJ', 'SAS', '-TS'],
        index=val_index)


alpha = 0.224
beta = 0.217
gamma = 0.0166*4.184*100 # 0.0166 kcal mol^-1 Å^-2 --> kJ mol^-1 nm^-2
tau = 0.0287
T = 298.15 # K

os.chdir('1pga')
n_struct = api.manage_CC_files('disco_1pga.pdb')
n_used = n_struct

G_wt_coulomb_mean = 0
G_wt_lj_mean = 0

for n in range(n_struct):
    os.chdir(str(n))
    min_filesize = api.emin('%d.pdb' % n, mdp, gmx_flags, table=tab)

    if min_filesize > 0:
        G_wt_coulomb, G_wt_lj = \
            api.extract_energies('%d_min.gro' % n, energy_mdp, gmx_flags)
        G_wt_coulomb_mean += G_wt_coulomb
        G_wt_lj_mean += G_wt_lj

    else:
        omitted.write("{0}, Nr. {1}".format('1pga', n))
        n_used -= 1

    os.chdir('..')

os.chdir('..')

G_wt_coulomb_mean /= n_used
G_wt_lj_mean /= n_used
G_wt_S_mean = api.calc_entropy('1pga')
G_wt_sasa_mean = api.sasa_mean('1pga')

values['COUL']['1pga'] = G_wt_coulomb_mean
values['LJ']['1pga'] = G_wt_lj_mean
values['SAS']['1pga'] = G_wt_sasa_mean
values['-TS']['1pga'] = G_wt_S_mean

n_muts = len(mutations.keys())
G_mut_coulomb = dict(zip(mutations.keys(), np.zeros(n_muts)))
G_mut_lj = dict(zip(mutations.keys(), np.zeros(n_muts)))
G_mut_S = dict(zip(mutations.keys(), np.zeros(n_muts)))
G_mut_sasa = dict(zip(mutations.keys(), np.zeros(n_muts)))

gxg_coulomb = pd.read_csv('../GXG/GXG_Coulomb.csv', index_col=0)
gxg_lj = pd.read_csv('../GXG/GXG_LJ.csv', index_col=0)
gxg_sasa = pd.read_csv('../GXG/GXG_Area.csv', index_col=0)
gxg_S = pd.read_csv('../GXG/GXG_S.csv', index_col=0)

for mut_name in mutations:
    f_name = '1pga_' + mut_name
    os.chdir(f_name)

    cc_name = 'disco_%s.pdb' % f_name
    n_struct = api.manage_CC_files(cc_name)
    n_used = n_struct

    for n in range(n_struct):
        os.chdir(str(n))
        min_filesize = api.emin('%d.pdb' % n, mdp, gmx_flags, tab)
        
        if min_filesize > 0:
            coulomb, lj = \
                api.extract_energies('%d_min.gro' % n, energy_mdp, gmx_flags)

        else:
            omitted.write("{0}, Nr. {1}".format(f_name, n))
            n_used -= 1

        G_mut_coulomb[mut_name] += coulomb
        G_mut_lj[mut_name] += lj
        os.chdir('..')

    G_mut_coulomb[mut_name] /= n_used
    G_mut_lj[mut_name] /= n_used
    os.chdir('..')
    G_mut_S[mut_name] = api.calc_entropy(f_name)
    G_mut_sasa[mut_name] = api.sasa_mean(f_name)

    values['COUL'][mut_name] = G_mut_coulomb[mut_name]
    values['LJ'][mut_name] = G_mut_lj[mut_name]
    values['SAS'][mut_name] = G_mut_sasa[mut_name]
    values['-TS'][mut_name] = G_mut_S[mut_name]
    values.to_csv('vals.csv')

    aa1 = mutations[mut_name][1]
    aa2 = mut_name[-1]
    print(G_wt_coulomb_mean, G_mut_coulomb[mut_name], gxg_coulomb[aa1][aa2]) 
    print(G_wt_lj_mean, G_mut_lj[mut_name], gxg_lj[aa1][aa2]) 
    print(G_wt_sasa_mean, G_mut_sasa[mut_name], gxg_sasa[aa1][aa2]) 
    print(G_wt_S_mean, G_mut_S[mut_name], gxg_S[aa1][aa2]) 
    ddG_fold_stability['COUL'][mut_name] = \
        alpha * (G_mut_coulomb[mut_name] - G_wt_coulomb_mean - \
        gxg_coulomb[aa2][aa1])
    ddG_fold_stability['LJ'][mut_name] = \
        beta * (G_mut_lj[mut_name] - G_wt_lj_mean - \
        gxg_lj[aa2][aa1])
    ddG_fold_stability['SAS'][mut_name] = \
        gamma * (G_mut_sasa[mut_name] - G_wt_sasa_mean - \
        gxg_sasa[aa2][aa1])
    ddG_fold_stability['-TS'][mut_name] = \
        tau * T * (G_mut_S[mut_name] - G_wt_S_mean - \
        gxg_S[aa2][aa1])

    ddG_fold_stability.to_csv('testing.csv')
    print(ddG_fold_stability)


omitted.close()
os.chdir('..')
