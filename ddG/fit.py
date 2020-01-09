from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

help(leastsq)

def fit_func(
    X,
    alpha,
    beta,
    gamma,
    tau
):
    a, b, c, d, e, f = X
    return alpha*(a+b) + beta*(c+d) + gamma*e + tau*f


p0 = np.array([0.224, 0.217, 6.9388, 0.0287])
prots = ['1ayi', '1pga', '1hz6', '2lzm', '3chy', '1ypc', '1stn']
calc_dfs = {}
exp_dfs = {}

for i in prots:
    calc_dfs[i] = pd.read_csv('%s.csv' % i, index_col=0)
    calc_dfs[i] = calc_dfs[i].iloc[:, 1:]
    calc_dfs[i].sort_index(inplace=True)
    exp_dfs[i] = pd.read_csv('%s-compare.csv' % i, index_col=0)
    exp_dfs[i] = exp_dfs[i].loc[calc_dfs[i].index]
    exp_dfs[i] = np.array(exp_dfs[i]['EXP'])
    exp_dfs[i] *= 4.18

exp = tuple(np.array(i) for i in exp_dfs.values())
exp = np.concatenate(exp)
calc = tuple(np.array(i) for i in calc_dfs.values())
calc = np.concatenate(calc)
calc = (calc[:, 0], calc[:, 1], calc[:, 2], calc[:, 3], calc[:, 4], calc[:, 5])
#coeffs, _ = leastsq(
#    fit_func,
#    calc,
#    exp,
#    *p0,
##    bounds=np.array((0, 10))
#)
#print(coeffs)
#x = fit_func(calc, *coeffs)
#print(np.corrcoef(exp, x))
#x = fit_func(calc, *p0)
#print(np.corrcoef(exp, x))
