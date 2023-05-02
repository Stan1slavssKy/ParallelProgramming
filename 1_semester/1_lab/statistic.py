import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from pyparsing import *
from statistics import mean

vals = 0
T_vals = np.empty(8)

def parse_file(vals, T_vals):
    X_param = 0
    T_param = 0
    with open('./build/statistic.txt') as f:
        counter = 0
        number = 0
        for line in f:
            line = line.split()
            
            if (line[0] == 'n:'):
                if (counter != 0):
                    T_vals[number - 1] = vals / counter
                    counter = 0
                    vals = 0
                number = int(line[-1])
            elif (line[0] == 'exec_time:'):
                vals = vals + float(line[-1])
                counter = counter + 1
            elif (line[0] == 'X:'):
                X_param = int(line[-1])
            elif (line[0] == 'T:'):
                T_param = int(line[-1])

        if (counter != 0):
            T_vals[number - 1] = vals / counter
    return X_param, T_param

X_param, T_param = parse_file(vals, T_vals)

p = np.arange(1, len(T_vals) + 1)
T0 = T_vals[0]
S = T0 / T_vals
E = S / p

plt.figure(figsize=[18, 10])
plt.title('Зависимость ускорения от числа процессов, X = ' + str(X_param) + ', T = ' + str(T_param))
plt.plot(p, S)
plt.xlabel("p")
plt.ylabel("S")
plt.grid()
plt.savefig('/home/stanislav/ParallelProgramming/1_lab/images/statistic_x_' + str(X_param) + '_t_' + str(T_param) + '.png')
plt.show()
