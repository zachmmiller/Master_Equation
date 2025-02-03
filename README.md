# Master Equation (MEQ) Usage + Compilation Instructions

Usage
MEQ is configured with one text file that defines all required variables. This file must be named “config.lua” for the program to recognize it and it must be formatted for the Lua programming language. Below is a list of all required variables with a brief description.

Required variables
t_min, t_max, t_step. Simulation time minimum, maximum and step size, respectively. The list of simulated times includes the end points. Given in seconds.

e_min, e_max, e_step. Simulation energy minimum, maximum and step size, respectively. The list of simulated energies includes the end points. Given in wavenumber (cm-1).

temperature. Temperature of the blackbody field. Given in Kelvin.

vib_modes. Array of vibrational normal modes. Given in wavenumber (cm-1). Must be the same length as vib_degen and ir_intens.

vib_degen. Array of vibrational normal mode degeneracies. These describe the degeneracy of each normal mode as an integer value. Must be the same length as vib_modes and ir_intens.

ir_intens. Array of infrared intensities. Given in km/mol. Must be the same length as vib_modes and vib_degen.

TS_vib_modes. Array of transition state normal modes. Given in wavenumber (cm-1). Must be the same length as TS_vib_degen.

TS_vib_degen. Array of transition state normal mode degeneracies. These describe the degeneracy of each normal mode as an integer value. Must be the same length as TS_vib_modes.

e0. Activation energy. Given in wavenumber (cm-1).

initially_boltzmann. Boolean. Dictates whether the initial population is a Boltzmann distribution (set to true) or user defined (set to false). If set to false, a Lua function named initial_population(index) must be defined. See below for more detail. 

initial_population(index). This function is defined by the user and called from C++. The integer parameter is the index of the array in C++, so it starts at 0 and ends at N-1.

number_of_threads. Integer. Sets the number of threads used when populating the transport matrix.

no_RRKM. Boolean. Set to true to prevent the program from computing RRKM dissociation rate constants.
save_modes. Boolean. Set to true to save vibrational modes and other associated, computed values. These values include vibrational energy, infrared intensity, degeneracy, Einstein A coefficient, Einstein B coefficient, Planck Distribution and the Einstein B coefficient multiplied by the Planck Distribution. 

save_initial_condition. Boolean. Set to true to save the initial condition as a .csv file. The file with have two columns. The first is the energy bins and the second is the abundance.

save_boltzmann. Boolean. Set to true to save a Boltzmann distribution (based upon the blackbody field temperature) as a .csv file. The file with have two columns. The first is the energy bins and the second is the abundance.

save_eigenvalues. Boolean. Set to true to save the eigenvalues of the transport matrix as a .csv file. The file with have four columns. The first is the index, the second is the real component, the third is the imaginary component and the last is the magnitude. 

save_time_data. Boolean. Set to true to save the time propagation data as a .csv file. The rows of the file correspond to each time point and the columns correspond to the energy bins. The first row is the values of the energy bins and the first column is the time for each point. 

Loading arrays from .csv files
Loading any file other than “config.lua” is done through the Lua scripting engine. A file titled “csv.lua” is included with the source code of MEQ that details Lua functions for loading .csv files into 2D tables and slicing tables into 1D arrays. The functions are csv.load and slice_csv, respectively. These functions are imported using the dofile(filepath) function in Lua, which behaves essentially like import in Python.

csv.load(filename, delimiter, header). The filename parameter is a string that is the full path to the file. The delimiter parameter is a string that details how the columns of the csv are separated (“,”, “\t”, etc.). The header parameter is a boolean and details whether the first row is a header or not (if the first row is not actual data, set this parameter to true). The data is arranged in memory such that rows in the csv are contiguous blocks (i.e. indexing is done by “some_csv[row][column]”).

csv_col(csv, c). The csv parameter is a Lua table and should be a 2D table (csv.load will return a 2D table even if the .csv file that you loaded only contains one column). The col parameter is an integer that describes what column is returned as a contiguous 1D table, which is suitable for parameters like vib_modes, ir_intens, etc.

csv_row(csv, r). The csv parameter is a Lua table and should be a 2D table (csv.load will return a 2D table even if the .csv file that you loaded only contains one row). The row parameter is an integer that describes what row is returned as a contiguous 1D table, which is suitable for parameters like vib_modes, ir_intens, etc.

Note on loading .csv files: The file must be formatted such that there is a new line character (“\n”) after the last column of the last row. If that new line character is not there, you will get this error: attempt to perform arithmetic on a nil value (local 'e'). To fix this, just open the .csv file in a plain text editor, press return on the last line, save and rerun the MEQ.

Additional symbols for scripting
Lua does not contain a library for file system related tasks out of the box, so things like the current working directory and the directory of the config file are not easily accessible. To get around this, C++ defines these two variables at run time before the config file is run and they can be accessed as config_directory and current_working_directory in the config file. See examples for more detail.

# Compiling
It should be more or less straightforward to compile MEQ on any system if you have a C++ compiler that supports C++20 (GCC 9, Clang 14 or MSVC 16 / VS2019). Once you have the right compiler, all you need is Premake, which is a project file generator for building C/C++ code. Premake will generate a makefile or a visual studio build file by calling “premake5 gmake” or “premake5 vs2022”, respectively in the base directory of the source code. From there, you can move on to compile the code with either build system. To compile with make, just type “make config=release” in the base directory and the C/C++ compiler should do the thing. As of right now, compilation has been tested on Mac, Windows and Linux with LLVM Clang++.



