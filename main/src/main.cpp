#define _USE_MATH_DEFINES

#include <math.h>

#include <algorithm>
#include <condition_variable>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

// #define EIGEN_USE_LAPACKE_STRICT
// #define EIGEN_USE_BLAS

#ifdef RELEASE
#define EIGEN_NO_DEBUG
#endif

#include "../../math/formulas.h"
#include "../../math/master_equation.h"
#include "../../math/states.h"
#include "../../tools/Argument_Parsing.h"
#include "../../tools/Progress_Bar.h"
#include "../../tools/ThreadPool.h"
#include "../../tools/Timer.h"
#include "../../tools/config.h"
#include "../../vendor/Eigen/Dense"
#include "../../vendor/Eigen/Eigenvalues"

namespace fs = std::filesystem;

constexpr int STORAGE_ORDER = Eigen::ColMajor;

int main(int argc, char** argv) {
    std::vector<arg_type> required_args = {{"-f", false}};
    std::vector<arg_type> optional_args = {{"-v", true}};
    std::vector<arg_value_pair> arguments;

    if (Parse_Arguments(argc, argv, &arguments, required_args, optional_args)) {
        std::cout << "Missing required arguments. Exiting." << std::endl;
        return 1;
    }

    bool verbose = false;
    fs::path folder;
    for (auto arg : arguments) {
        if (arg.name == "-f") {
            folder = fs::canonical(arg.value);
        }

        if (arg.name == "-v") {
            verbose = true;
        }
    }

    std::cout << "\nfolder: " << folder << std::endl;

    fs::path config_file;
    fs::directory_iterator files(folder);
    for (const fs::directory_entry& entry : files) {
        if (entry.path().filename().string() == "config.lua") {
            config_file = entry.path();
        }
    }
    std::cout << "config file: " << config_file << std::endl;

    std::string output_directory_str;
    fs::path output_directory;
    double t_min, t_max, t_step, temperature, initial_temperature, initial_integral_abundance;
    int e_min, e_max, e_step, e0, number_of_threads;
    std::vector<int> vib_modes, vib_degen, TS_vib_modes, TS_vib_degen;
    std::vector<double> ir_intens;
    bool initially_boltzmann, no_RRKM, save_initial_condition, save_bins, save_modes, save_boltzmann, save_time_data, save_eigenvalues;
    lua_State* L = luaL_newstate();
    {
        std::cout << "\nLoading configuration..." << std::endl;
        Timer timer("Took");
        luaL_openlibs(L);

        fs::path cwd(fs::canonical("."));
        if (luaL_dostring(L, std::format("current_working_directory = \"{}\"", cwd.string().c_str()).c_str()) != LUA_OK) {
            std::cout << "Error setting current working directory. Message:" << std::endl;
            std::cout << lua_tostring(L, -1) << std::endl;
            lua_close(L);
            return 1;
        }

        if (luaL_dostring(L, std::format("config_directory = \"{}\"", config_file.parent_path().string().c_str()).c_str()) != LUA_OK) {
            std::cout << "Error setting config directory. Message:" << std::endl;
            std::cout << lua_tostring(L, -1) << std::endl;
            lua_close(L);
            return 1;
        }

        if (luaL_dofile(L, config_file.string().c_str()) != LUA_OK) {
            std::cout << "Error running config file. Message:" << std::endl;
            std::cout << lua_tostring(L, -1) << std::endl;
            lua_close(L);
            return 1;
        }

        lua_getglobal(L, "output_directory");
        if (Lua_Load_String(L, &output_directory_str)) {
            std::cout << "Required parameter \"output_directory\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "t_min");
        if (Lua_Load_Number<double>(L, &t_min)) {
            std::cout << "Required parameter \"t_min\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "t_max");
        if (Lua_Load_Number<double>(L, &t_max)) {
            std::cout << "Required parameter \"t_max\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "t_step");
        if (Lua_Load_Number<double>(L, &t_step)) {
            std::cout << "Required parameter \"t_step\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "temperature");
        if (Lua_Load_Number<double>(L, &temperature)) {
            std::cout << "Required parameter \"temperature\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "e_min");
        if (Lua_Load_Number<int>(L, &e_min)) {
            std::cout << "Required parameter \"e_min\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "e_max");
        if (Lua_Load_Number<int>(L, &e_max)) {
            std::cout << "Required parameter \"e_max\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "e_step");
        if (Lua_Load_Number<int>(L, &e_step)) {
            std::cout << "Required parameter \"e_step\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "vib_modes");
        if (Lua_Load_1d_Vector<int>(L, vib_modes)) {
            std::cout << "Required parameter \"vib_modes\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "vib_degen");
        if (Lua_Load_1d_Vector<int>(L, vib_degen)) {
            std::cout << "Required parameter \"vib_degen\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "ir_intens");
        if (Lua_Load_1d_Vector<double>(L, ir_intens)) {
            std::cout << "Required parameter \"ir_intens\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "TS_vib_modes");
        if (Lua_Load_1d_Vector<int>(L, TS_vib_modes)) {
            std::cout << "Required parameter \"TS_vib_modes\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "TS_vib_degen");
        if (Lua_Load_1d_Vector<int>(L, TS_vib_degen)) {
            std::cout << "Required parameter \"TS_vib_degen\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "e0");
        if (Lua_Load_Number<int>(L, &e0)) {
            std::cout << "Required parameter \"e0\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "initially_boltzmann");
        if (Lua_Load_Bool(L, &initially_boltzmann)) {
            std::cout << "Required parameter \"initially_boltzmann\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "initial_temperature");
        if (Lua_Load_Number<double>(L, &initial_temperature)) {
            std::cout << "Required parameter \"initial_temperature\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "initial_integral_abundance");
        if (Lua_Load_Number<double>(L, &initial_integral_abundance)) {
            std::cout << "Required parameter \"initial_integral_abundance\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "no_RRKM");
        if (Lua_Load_Bool(L, &no_RRKM)) {
            std::cout << "Required parameter \"no_RRKM\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "save_modes");
        if (Lua_Load_Bool(L, &save_modes)) {
            std::cout << "Required parameter \"save_modes\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "save_bins");
        if (Lua_Load_Bool(L, &save_bins)) {
            std::cout << "Required parameter \"save_bins\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "save_initial_condition");
        if (Lua_Load_Bool(L, &save_initial_condition)) {
            std::cout << "Required parameter \"save_initial_condition\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "save_boltzmann");
        if (Lua_Load_Bool(L, &save_boltzmann)) {
            std::cout << "Required parameter \"save_boltzmann\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "save_eigenvalues");
        if (Lua_Load_Bool(L, &save_eigenvalues)) {
            std::cout << "Required parameter \"save_eigenvalues\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "save_time_data");
        if (Lua_Load_Bool(L, &save_time_data)) {
            std::cout << "Required parameter \"save_time_data\" not found. Exiting." << std::endl;
            return 1;
        }

        lua_getglobal(L, "number_of_threads");
        if (Lua_Load_Number<int>(L, &number_of_threads)) {
            std::cout << "Required parameter \"number_of_threads\" not found. Exiting." << std::endl;
            return 1;
        }
    }

    output_directory = fs::canonical(output_directory_str);
    if (!fs::is_directory(output_directory)) {
        std::cout << "Output directory is not valid. Exiting." << std::endl;
        return 1;
    }

    if (vib_modes.size() != vib_degen.size() || vib_degen.size() != ir_intens.size()) {
        std::cout << "Parameters vib_modes, vib_degen and ir_intens do not all have the same size. Exiting." << std::endl;
        return 1;
    }

    if (TS_vib_modes.size() != TS_vib_degen.size()) {
        std::cout << "Parameters TS_vib_modes and TS_vib_degen do not have the same size. Exiting." << std::endl;
        return 1;
    }
    {
        int n_modes = 0;
        for (int& i : vib_degen) {
            n_modes += i;
        }
        int n_TS_modes = 0;
        for (int& i : TS_vib_degen) {
            n_TS_modes += i;
        }
        if (n_TS_modes != n_modes - 1 && !no_RRKM) {
            std::cout << "Warning: The number of transition state vibrational modes should\nbe one less than the number of reactant complex vibrational modes.\n" << std::endl;
        }
    }

    int n_t_steps = (int)((t_max - t_min) / t_step) + 1;
    int n_e_bins = ((e_max - e_min) / e_step) + 1;

    Vibrational_Modes modes(vib_modes.data(), vib_degen.data(), ir_intens.data(), temperature, vib_modes.size());

    if (verbose) {
        std::cout << "Blackbody field temperature: " << temperature << " K" << std::endl;
        std::cout << "Initial temperature: " << initial_temperature << " K" << std::endl;

        modes.print();
    }

    if (save_modes) {
        std::cout << "Saving vibrational modes..." << std::endl;
        Timer timer("Took");
        fs::path outpath(output_directory / "modes.csv");
        std::ofstream outfile;
        outfile.open(outpath);
        std::string out_str("Index, mode (cm^-1), degeneracy, IR intensity (km/mol), A (s^-1), B (cm^3 * J^-1 * s^-2), P (J * s * cm^-3), B * P (s^-1)\n");
        outfile << out_str;
        for (int i = 0; i < modes.N; i++) {
            out_str.clear();
            out_str += std::format("{:10d}, {:10d}, {:10d}, {:18.15f}, {:18.15f}, {:18.15e}, {:18.15f}, {:18.15f}",
                                   i,
                                   modes.C[i],
                                   modes.D[i],
                                   modes.I[i],
                                   modes.A[i],
                                   modes.B[i],
                                   modes.P[i],
                                   modes.B[i] * modes.P[i]) +
                       "\n";
            outfile << out_str;
        }
        outfile << out_str;
    }

    // These are all of the linear algebra things we need.

    // Initial population.
    Eigen::Vector<std::complex<double>, -1> N_0 = Eigen::Vector<std::complex<double>, -1>::Zero(n_e_bins);

    // Population at time t
    Eigen::Vector<std::complex<double>, -1> N_t = Eigen::Vector<std::complex<double>, -1>::Zero(n_e_bins);

    // Eigenvalues
    Eigen::Vector<std::complex<double>, -1> e_val = Eigen::Vector<std::complex<double>, -1>::Zero(n_e_bins);

    // Transport matrix
    Eigen::Matrix<double, -1, -1, STORAGE_ORDER> J = Eigen::Matrix<double, -1, -1, STORAGE_ORDER>::Zero(n_e_bins, n_e_bins);

    // Eigenvectors
    Eigen::Matrix<std::complex<double>, -1, -1, STORAGE_ORDER> e_vec = Eigen::Matrix<std::complex<double>, -1, -1, STORAGE_ORDER>::Zero(n_e_bins, n_e_bins);

    // Inverse of eigenvectors
    Eigen::Matrix<std::complex<double>, -1, -1, STORAGE_ORDER> e_vec_inv = Eigen::Matrix<std::complex<double>, -1, -1, STORAGE_ORDER>::Zero(n_e_bins, n_e_bins);

    // Propagator (e^-lambda*t as a diagonal matrix)
    Eigen::Matrix<std::complex<double>, -1, -1, STORAGE_ORDER> lam_mat = Eigen::Matrix<std::complex<double>, -1, -1, STORAGE_ORDER>::Zero(n_e_bins, n_e_bins);

    if (save_bins) {
        std::cout << "Saving bins..." << std::endl;
        Timer timer("Took");
        fs::path outpath(output_directory / "bins.csv");
        std::ofstream outfile;
        outfile.open(outpath);
        std::string out_str;
        for (int i = 0; i < n_e_bins; i++) {
            out_str.clear();
            out_str += std::format("{:8d}\n", e_step * i + e_min);
            outfile << out_str;
        }
    }

    {
        std::cout << "Calculating initial population..." << std::endl;
        Timer timer("Took");

        if (save_boltzmann) {
            Boltzmann(reinterpret_cast<std::complex<double>*>(N_0.data()), n_e_bins, e_min, e_step, vib_modes.data(), vib_degen.data(), vib_modes.size(), temperature);
            fs::path outpath(output_directory / "boltzmann.csv");
            std::ofstream outfile;
            outfile.open(outpath);
            std::string out_str;
            double value;
            for (int i = 0; i < n_e_bins; i++) {
                out_str.clear();
                value = sqrt(N_0.data()[i].real() * N_0.data()[i].real() + N_0.data()[i].imag() * N_0.data()[i].imag());
                out_str += std::format("{:18.15f}\n", value);
                outfile << out_str;
            }
            N_0 *= initial_integral_abundance;
        }

        if (initially_boltzmann) {
            Boltzmann(reinterpret_cast<std::complex<double>*>(N_0.data()), n_e_bins, e_min, e_step, vib_modes.data(), vib_degen.data(), vib_modes.size(), initial_temperature);
            N_0 *= initial_integral_abundance;
        }

        if (!initially_boltzmann) {
            double sum = 0;
            for (int i = 0; i < n_e_bins; i++) {
                lua_getglobal(L, "initial_population");
                lua_pushnumber(L, i);
                if (lua_pcall(L, 1, 1, 0) != 0) {
                    std::cout << "Error calling initial_population. Exiting." << std::endl;
                    return 1;
                }
                if (!lua_isnumber(L, -1)) {
                    std::cout << "Value returned by initial_population is not a number. Exiting." << std::endl;
                    return 1;
                }
                N_0.data()[i] = std::complex<double>(lua_tonumber(L, -1), 0);
                sum += N_0.data()[i].real();
                lua_pop(L, 1);
            }

            // Normalize and scale to initial_integral_abundance
            for (int i = 0; i < n_e_bins; i++) {
                N_0.data()[i] = std::complex<double>(initial_integral_abundance * N_0.data()[i].real() / sum, 0);
            }
        }

        if (save_initial_condition) {
            fs::path outpath(output_directory / "initial_condition.csv");
            std::ofstream outfile;
            outfile.open(outpath);
            std::string out_str;
            double value;
            for (int i = 0; i < n_e_bins; i++) {
                out_str.clear();
                value = sqrt(N_0.data()[i].real() * N_0.data()[i].real() + N_0.data()[i].imag() * N_0.data()[i].imag());
                out_str += std::format("{:18.15f}\n", value);
                outfile << out_str;
            }
        }
    }

    ThreadPool pool;
    pool.Start(number_of_threads, [] {});
    {
        std::cout << "Calculating blackbody absorption / emission rates..." << std::endl;
        Timer timer("Took");
        Progress_Bar bar(n_e_bins, 50, "");
        bar.Start();
        std::mutex bar_mutex;
        std::mutex busy_mutex;
        std::condition_variable busy_cv;
        for (long int i = 0; i < n_e_bins; i++) {
            pool.QueueJob(
                [&modes, e_min, e_step, &J, i, n_e_bins, &bar, &bar_mutex] {
                    J_Column_Blackbody(&modes, e_min, e_step, reinterpret_cast<double*>(J.data()), i, (long int)n_e_bins);
                    {
                        std::unique_lock<std::mutex> lock(bar_mutex);
                        bar.Update();
                    }
                },
                &busy_cv);
        }
        {
            std::unique_lock<std::mutex> lock(busy_mutex);
            busy_cv.wait(lock, [&pool] { return !pool.Busy(); });
        }

        bar.End();
    }
    /*
    // Testing RRKM rate constant with ethyl radical dissociation
    {
        double density_of_states_act = Density_of_States(vib_modes.data() , vib_degen.data(), vib_modes.size(), (long int)14629);
        double sum_of_states_trans = Sum_of_States(TS_vib_modes.data(), TS_vib_degen.data(), TS_vib_modes.size(), (long int)(1588));
        double k_RRKM = sum_of_states_trans / (density_of_states_act * form::h * (1.0 / (form::h * form::c * 100)));
        std::cout << "kRRKM: " << k_RRKM << std::endl;
    }
    */

    if (!no_RRKM) {
        std::cout << "Calculating RRKM dissociation rate constants..." << std::endl;
        Timer timer("Took");
        Progress_Bar bar(n_e_bins, 50, "");
        bar.Start();

        for (long int i = 0; i < n_e_bins; i++) {
            J_Column_RRKM(vib_modes.data(),
                          vib_degen.data(),
                          vib_modes.size(),
                          TS_vib_modes.data(),
                          TS_vib_degen.data(),
                          TS_vib_modes.size(),
                          e0,
                          e_min,
                          e_step,
                          reinterpret_cast<double*>(J.data()),
                          i,
                          (long int)n_e_bins);
            bar.Update();
        }
        bar.End();
    }

    {
        std::cout << "Solving for transport matrix eigenvalues / eigenvectors..." << std::endl;
        Timer timer("Took");
        // Solver for eigenvalues / eigenvectors (could be replaced by LAPACK DGEEV for better performance. This does not seem trivial to do unfortunately.)
        Eigen::EigenSolver<Eigen::Matrix<double, -1, -1, STORAGE_ORDER>> solver(n_e_bins);
        solver.setMaxIterations(10000000);
        solver.compute(J);
        if (solver.info() != Eigen::Success) {
            std::cout << "Solver failed. Exiting." << std::endl;
            return 1;
        }

        e_val = solver.eigenvalues();
        e_vec = solver.eigenvectors();
    }

    {
        std::cout << "Inverting eigenvectors..." << std::endl;
        Timer timer("Took");
        e_vec_inv = e_vec.inverse();
    }

    if (save_eigenvalues) {
        std::cout << "Saving eigenvalues..." << std::endl;
        Timer timer("Took");
        fs::path outpath(output_directory / "eigenvalues.csv");
        std::ofstream outfile;
        outfile.open(outpath);
        std::string out_str("Real component, Imaginary component");
        outfile << out_str;
        for (int i = 0; i < n_e_bins; i++) {
            out_str.clear();
            out_str += std::format("{:18d}, {:18.15f}, {:18.15f}\n", i, e_val.data()[i].real(), e_val.data()[i].imag());
            outfile << out_str;
        }
    }

    N_0 = e_vec_inv * N_0;

    // Uncomment this to remove long term decay / growth of the population that is likely due to floating point rounding in the transport matrix.
    /*
    for (int i = 0; i < n_e_bins; i++) {
        if (abs(e_val.data()[i].real()) < 1e-6 || e_val.data()[i].real() > 0) {
            e_val.data()[i] = std::complex(0, 0);
        }
    }
    */

    if (save_time_data) {
        fs::path outpath(output_directory / "time_data.csv");
        std::ofstream outfile;
        outfile.open(outpath);

        std::cout << "Time propagation..." << std::endl;
        Timer timer("Took");
        double time = 0;
        std::string out_str;
        Progress_Bar bar(n_t_steps, 50, "");
        bar.Start();
        for (int i = 0; i < n_t_steps; i++) {
            time = i * t_step + t_min;
            std::complex<double> ev;
            for (int i = 0; i < n_e_bins; i++) {
                ev = e_val.data()[i];
                std::complex<double> value(exp(ev.real() * time) * cos(ev.imag() * time), exp(ev.real() * time) * sin(ev.imag() * time));
                lam_mat.data()[i * n_e_bins + i] = value;
            }

            N_t = lam_mat * N_0;
            N_t = e_vec * N_t;

            // Below is equivalent to what is going on above, but is way of doing it is much slower.
            // N_t = e_vec * lam_mat * e_vec_inv * N_0;

            if (save_time_data) {
                out_str.clear();
                out_str += std::format("{:18.15f}, ", time);  // format: 18 is the length (in chars) and .15f is the precision
                double value = 0;
                for (int j = 0; j < n_e_bins; j++) {
                    value = sqrt(N_t.data()[j].real() * N_t.data()[j].real() + N_t.data()[j].imag() * N_t.data()[j].imag());
                    out_str += std::format("{:18.15f}", value);
                    if (j < n_e_bins - 1) {
                        out_str += ", ";
                    }
                }
                outfile << out_str << "\n";
            }
            bar.Update();
        }
        bar.End();
    }

    pool.Stop();
    lua_close(L);
    return 0;
}
