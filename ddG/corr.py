import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csvs = ['1ayi', '1ypc', '3chy', '1pga', '2lzm', '1hz6', '1stn']

for c in csvs:
    x = pd.read_csv(c + '-compare.csv', index_col=0)
    x *= 4.18
    y = pd.read_csv(c + '.csv', index_col=0)
    y['LJ'] = y['LJ (SR)'] + y['LJ (SR)']
    x.sort_index(inplace=True)
    y.sort_index(inplace=True)

    for col in x.columns[1:]:
        corrcoef = np.corrcoef(x[col], y[col])[0][1]
        line = np.arange(
            min(min(x[col]), min(y[col])),
            max(max(x[col]), max(y[col]))
        )
        plt.plot(x[col], y[col], '.')
        plt.plot(line, line)
        plt.text(
            line[0],
            line[-1],
            "Correlation coefficient of %s %s: %f" % (c, col, float(corrcoef))
        )
        plt.ylabel("$\Delta \Delta G_{Linkai}\\,[\\frac {kJ} {mol}]$")
        plt.xlabel("$\Delta \Delta G_{Paper}\\,[\\frac {kJ} {mol}]$")
        plt.show()

#y = pd.read_csv('1ayi' + '.csv', index_col=0)
#print(max(y['LJ (1-4)']))
