import csv
import numpy as np

def get_time_data(folder):
    file = f"{folder}/time_data.csv"

    energy = []
    time = []
    bins = []

    with open(file, newline="") as fo:
        reader = csv.reader(fo, delimiter=",")
        for row in reader:
            if len(energy) == 0:
                for i in range(1, len(row)):
                    energy.append(int(row[i]))
            else:
                time.append(float(row[0]))
                bins.append([])
                for i in range(1, len(row)):
                    bins[-1].append(float(row[i]))
    
    energy = np.asarray(energy)
    time = np.asarray(time)
    bins = np.asarray(bins)

    return time, energy, bins


def get_boltzmann(folder):
    file = f"{folder}/boltzmann.csv"

    energy = []
    bins = []

    with open(file, newline="") as fo:
        reader = csv.reader(fo, delimiter=",")
        next(reader)
        for row in reader:
            energy.append(int(row[0]))
            bins.append(float(row[1]))
    
    energy = np.asarray(energy)
    bins = np.asarray(bins)

    return energy, bins

def get_initial_condition(folder):
    file = f"{folder}/initial_condition.csv"

    energy = []
    bins = []

    with open(file, newline="") as fo:
        reader = csv.reader(fo, delimiter=",")
        next(reader)
        for row in reader:
            energy.append(int(row[0]))
            bins.append(float(row[1]))
    
    energy = np.asarray(energy)
    bins = np.asarray(bins)

    return energy, bins

