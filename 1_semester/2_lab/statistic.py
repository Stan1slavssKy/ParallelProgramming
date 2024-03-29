import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from pyparsing import *
from statistics import mean

vals = 0
T_vals = np.empty(16)

def parse_file(vals, T_vals):
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
            elif (line[-1] == 'us'):
                vals = vals + float(line[0])
                counter = counter + 1

        if (counter != 0):
            T_vals[number - 1] = vals / counter
    return

parse_file(vals, T_vals)

p = np.arange(1, len(T_vals) + 1)
T0 = T_vals[0]
S = T0 / T_vals
E = S / p

plt.figure(figsize=[18, 10])
plt.title('Зависимость ускорения от числа процессов для функции sin(1 / (x + 20)) на промежутке [-19.99, 10.0] шаг 1e-11')
plt.plot(p, S)
plt.xlabel("p")
plt.ylabel("S")
plt.grid()
plt.savefig('/home/stanislav/ParallelProgramming/1_semester/2_lab/images/statistic_2.png')
plt.show()
