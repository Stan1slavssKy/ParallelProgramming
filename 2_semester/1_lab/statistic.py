import argparse
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from pyparsing import *
from statistics import mean

def draw_graph(thread_time_arr, plot_name):
    p = np.arange(1, len(thread_time_arr) + 1)
    T0 = thread_time_arr[0]
    S = T0 / thread_time_arr
    E = S / p

    plt.figure(figsize=[18, 10])
    plt.title('Зависимость ускорения от числа исполнителей для ' + plot_name)
    plt.plot(p, S)
    plt.xlabel("p")
    plt.ylabel("S")
    plt.grid()
    plt.savefig('/home/stanislav/ParallelProgramming/2_semester/1_lab/pictures/statistic_' + plot_name + '.png')
    plt.show()

def parse_file(thread_number, filename):
    vals = 0
    T_vals = np.empty(int(thread_number))

    with open(f'./build/{filename}') as f:
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

        if (counter != 0):
            T_vals[number - 1] = vals / counter
    return T_vals

def parse_args():
    parser = argparse.ArgumentParser(prog='statistic.py')

    parser.add_argument('thread_number', help='max number of threads')
    parser.add_argument('statistic_filename',  help='name of the statistic file')
    parser.add_argument('picture_name',  help='name of the statistic picture')

    return parser.parse_args()

def main():
    args = parse_args()
    T_vals = parse_file(args.thread_number, args.statistic_filename)
    draw_graph(T_vals, args.picture_name)

if __name__ == '__main__':
    main()
