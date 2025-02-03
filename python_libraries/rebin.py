import argparse as ap
import numpy as np
from scipy.interpolate import CubicSpline
import os
import csv


def main():
    parser = ap.ArgumentParser()
    parser.add_argument("-min", type=int, help="Energy minimum in cm^-1", required=True)
    parser.add_argument("-max", type=int, help="Energy maximum in cm^-1", required=True)
    parser.add_argument("-step", type=int, help="Energy step size in cm^-1", required=True)
    parser.add_argument("-mult", type=float, help="Value that is multiplied by the input energies to convert to wavenumber. Default is 1 eV = 8065.544 cm^-1", default=8065.544)
    parser.add_argument("-input", type=str, help="Relative path to input .csv file", required=True)
    parser.add_argument("-output", type=str, help="Relative path to save output .csv file", required=True)
    parser.add_argument("-header", type=bool, help="Does the .csv have a header?", default=True)

    args = parser.parse_args()
    e_min = args.min
    e_max = args.max
    e_step = args.step
    e_conversion = args.mult
    header = args.header
    input_file = os.getcwd() + "/" + args.input
    output_file = os.getcwd() + "/" + args.output

    energy = []
    bins = []

    with open(input_file, newline="") as fo:
        reader = csv.reader(fo, delimiter=",")
        if header:
            next(reader)
        for row in reader:
            energy.append(float(row[0]))
            bins.append(float(row[1]))

    energy = np.asarray(energy) * e_conversion
    bins = np.asarray(bins)

    cs = CubicSpline(energy, bins)
    new_energy = np.arange(e_min, e_max + e_step, e_step, dtype=int)
    new_bins = cs(new_energy)

    with open(output_file, "w", encoding="UTF8", newline="") as fo:
        writer = csv.writer(fo)
        writer.writerow(["Energy (cm^-1)", "Abundance"])
        for i in range(0, len(new_bins)):
            writer.writerow([new_energy[i], new_bins[i]])


if __name__ == "__main__":
    main()
