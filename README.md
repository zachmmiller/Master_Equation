# Master Equation (MEQ) Usage + Compilation Instructions

# Usage
MEQ is configured with one text file that defines all required variables. This file must be named “config.lua” for the program to recognize it and it must be formatted for the Lua programming language. Below is a list of all required variables with a brief description.

## Required variables
Below is a list of all required variables for running MEQ that have to be defined somehow in the config.lua file. Names are case sensitive.  

 - t_min, t_max, t_step. 
    Simulation time minimum, maximum and step size, respectively. The list of simulated times includes the end points. Given in seconds.

 - e_min, e_max, e_step. 
    Simulation energy minimum, maximum and step size, respectively. The list of simulated energies includes the end points. Given in wavenumber (cm-1).

 - temperature. 
    Temperature of the blackbody field. Given in Kelvin.

 - vib_modes. 
    Array of vibrational normal modes. Given in wavenumber (cm-1). Must be the same length as vib_degen and ir_intens.

 - vib_degen. 
    Array of vibrational normal mode degeneracies. These describe the degeneracy of each normal mode as an integer value. Must be the same length as vib_modes and ir_intens.

 - ir_intens. 
    Array of infrared intensities. Given in km/mol. Must be the same length as vib_modes and vib_degen.

 - TS_vib_modes. 
    Array of transition state normal modes. Given in wavenumber (cm-1). Must be the same length as TS_vib_degen.

 - TS_vib_degen. 
    Array of transition state normal mode degeneracies. These describe the degeneracy of each normal mode as an integer value. Must be the same length as TS_vib_modes.

 - e0. 
    Activation energy. Given in wavenumber (cm-1).

 - initially_boltzmann. 
    Boolean. Dictates whether the initial population is a Boltzmann distribution (set to true) or user defined (set to false). If set to false, a Lua function named initial_population(index) must be defined. See below for more detail. 

 - initial_population(index). 
    This function is defined by the user and called from C++. The integer parameter is the index of the array in C++, so it starts at 0 and ends at N-1. Remember that Lua arrays start at 1, so if you are loading data in from a .csv in Lua and passing it to C++ via this function, you will have to compensate by adding 1 to index. 

 - number_of_threads. 
    Integer. Sets the number of threads used when populating the transport matrix.

 - no_RRKM. 
    Boolean.Set to true to prevent the program from computing RRKM dissociation rate constants.

 - save_modes. 
    Boolean. Set to true to save vibrational modes and other associated, computed values. These values include vibrational energy, infrared intensity, degeneracy, Einstein A coefficient, Einstein B coefficient, Planck Distribution and the Einstein B coefficient multiplied by the Planck Distribution. 

 - save_initial_condition. 
    Boolean. Set to true to save the initial condition as a .csv file. The file with have two columns. The first is the energy bins and the second is the abundance.

 - save_boltzmann. 
    Boolean. Set to true to save a Boltzmann distribution (based upon the blackbody field temperature) as a .csv file. The file with have two columns. The first is the energy bins and the second is the abundance.

 - save_eigenvalues. 
    Boolean. Set to true to save the eigenvalues of the transport matrix as a .csv file. The file with have four columns. The first is the index, the second is the real component, the third is the imaginary component and the last is the magnitude. 

 - save_time_data. 
    Boolean. Set to true to save the time propagation data as a .csv file. The rows of the file correspond to each time point and the columns correspond to the energy bins. The first row is the values of the energy bins and the first column is the time for each point. 

## Lua Libraries
Loading .csv files for inputs to MEQ can be done with Lua. A file titled “csv.lua” is included with the source code of MEQ that includes functions for loading .csv files into 2D tables and slicing 2D tables into 1D tables. The functions are csv.load and csv_col/csv_row, respectively. These functions are imported using the dofile(filepath) function, which behaves essentially like import in Python.

 - csv.load(filename, delimiter, header). The filename parameter is a string that is the full path to the file. The delimiter parameter is a string that details how the columns of the csv are separated (“,”, “\t”, etc.). The header parameter is a boolean and details whether the first row is a header or not (if the first row is not actual data, set this parameter to true). The data is arranged in memory such that rows in the csv are contiguous blocks (i.e. indexing is done by “some_csv[row][column]”).

 - csv_col(csv, c). The csv parameter is a Lua table and should be a 2D table (csv.load will return a 2D table even if the .csv file that you loaded only contains one column). The col parameter is an integer that describes what column is returned as a contiguous 1D table, which is suitable for parameters like vib_modes, ir_intens, etc.

 - csv_row(csv, r). The csv parameter is a Lua table and should be a 2D table (csv.load will return a 2D table even if the .csv file that you loaded only contains one row). The row parameter is an integer that describes what row is returned as a contiguous 1D table, which is suitable for parameters like vib_modes, ir_intens, etc.

Note on loading .csv files: The file must be formatted such that there is a new line character (“\n”) after the last column of the last row. If that new line character is not there, you will get this error: "attempt to perform arithmetic on a nil value (local 'e')". To fix this, just open the .csv file in a plain text editor, press return on the last line, save and rerun MEQ.

## Additional symbols for scripting
Lua does not contain a library for file system related tasks out of the box, so things like the current working directory and the directory of the config file are not easily accessible. To get around this, C++ defines these two variables at run time before the config file is run and they can be accessed as config_directory and current_working_directory in the config file. See examples for more detail.

## Python Libraries

There are currently two scripts in the python_libraries directory. The first is importing.py, which provides some functions for plotting the outputs of this program and the second is rebin.py, which is a short script that will rebin a distribution given in eV (or any other unit of energy) to a distribution in wavenumber.

 - importing.py. There are functions three functions for loading saved time data, boltzmann distributions and initial conditions. All functions take one parameter, which is the folder that the saved data is in.
   - get_time_data loads in saved time propagation data and returns three numpy.ndarrays. The first is time (1D array), the second is energy (1D array) and the third is the abundance (2D array, first dimension is time and the second is energy).
   - get_boltzmann loads in a saved boltzmann distribution and returns two numpy.ndarrays. The first is energy (1D array) and the second is abundance (1D array).
   - get_initial_condition loads in a saved initial condition and returns two numpy.ndarrays. The first is energy (1D array) and the second is abundance (1D array).

 - rebin.py. This is a script meant to be run from the command line. It takes 7 arguments and is used to rebin distributions for use with MEQ. The arguments are:
   - "-min". Integer. Required. Sets the energy minumum for binning in wavenumber.
   - "-max". Integer. Required. Sets the energy maximum for binning in wavenunumber.
   - "-step". Integer. Required. Sets the energy step for binning in wavenumber.
   - "-mult". Float. Optional. Default=8065.544 (converts eV to wavenumber). Sets the conversion factor for input energy units to wavenunumber.
   - "-input". String. Required. File path to input data. Must be a two column .csv file. First column is energy and the second is abundance.
   - "-output". String. Required. File path for output. 
   - "-header". Boolean. Optional. Default=True. Set to false if the first row of the input is not a header. Keep as true if the first row is a header.

   How it works: The input distribution is first converted to wavenumber and then the CubicSpline class from scipy.interpolate is used to rebin the distribution in evenly spaced wavenumber bins. Note that renormalization is not performed here, but MEQ will renormalize the initial condition. 

## Examples

There are 10 examples configurations of the program in the examples directory that go over various ways of running the MEQ with Calcium Hexahydrate 2+ as a test system:
 - 80C 
 - 90C
 - 100C
 - 110C
 - 120C
 - 200C
   - Simulates the dissociation rate of Ca(H2O)6 2+ at different temperatures. The initial population is defined by a boltzmann distribution.

 - cooling_1
   - Simulates cooling of Ca(H2O)6 2+ from 200C to 50C by starting the population as a 200C Boltzmann distribution.
    
 - cooling_2
   - Simulates cooling of Ca(H2O)6 2+ from 200C to 50C by starting the population as a 200C depleted Boltzmann distribution from the 200C simulation above (you must run 200C before this one for this one to work). 

 - cooling_3
    - Simulates cooling of Ca(H2O)6 2+ from 200C to 50C by starting the population as a 200C Boltzmann distribution. RRKM rate constants are not computed for this one.

 - boltzmann_comparison
   - Simulates Ca(H2O)6 2+ at 400C without dissociation for comparison to Boltzmann distribution. 


# Compiling
It should be more or less straightforward to compile MEQ on any system if you have a C++ compiler that supports C++20 (GCC 9, Clang 14 or MSVC 16 / VS2019). Once you have a compiler that supports C++20, the program is ready to build. You have a choice of using GNU make (makefiles are provided) or any other build system supported by premake5. Premake5 is a project file generator for building C/C++ code. Premake will generate makefiles or a visual studio build files by calling “premake5 gmake” or “premake5 vs2022”, respectively in the base directory of the source code. To compile with make, just type “make config=release” in the base directory and the C/C++ compiler should do the rest and put an executable in the bin directory. As of right now, compilation has been tested on Mac, Windows and Linux with LLVM Clang++, but GCC or MSVC should work without a problem. To request gcc for example, just change "toolset" in premake5.lua from "clang" to "gcc" and regenerate the build files.


