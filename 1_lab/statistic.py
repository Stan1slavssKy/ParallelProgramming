import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from pyparsing import *
from statistics import mean
vals = 0
T_vals = np.empty(8)

def parse_file(vals, T_vals):
    with open('statistic.txt') as f:
        counter = 0
        number = 0
        for line in f:
            line = line.split()
            
            if (line[0] == 'n:'):
                if (counter != 0):
                    T_vals[number - 1] = vals / counter
                    print(number - 1, vals / counter)
                    counter = 0
                    vals = 0
                number = int(line[-1])
                print("number = ", number)
            elif (line[0] == 'exec_time:'):
                vals = vals + float(line[-1])
                print(vals)
                counter = counter + 1

        if (counter != 0):
            T_vals[number - 1] = vals / counter
            print(number - 1, vals / counter)
    return

parse_file(vals, T_vals)
print(T_vals)

p = np.arange(1, len(T_vals) + 1)
T0 = T_vals[0]
S = T0 / T_vals
E = S / p

plt.figure(figsize=[18, 10])
plt.title("Зависимость ускорения от числа процессов")
plt.plot(p, S)
plt.xlabel("p")
plt.ylabel("S")
plt.grid()
plt.savefig("/home/stanislav/ParallelProgramming/1_lab/images/statistic.png")
plt.show()