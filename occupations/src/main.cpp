//
// Created by Zach Miller on 7/1/25.
//

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

#include <cmath>

#include "../../math/states.h"
#include "../../tools/Argument_Parsing.h"
#include "../../tools/Progress_Bar.h"
#include "../../tools/ThreadPool.h"
#include "../../tools/Timer.h"
#include "../../tools/Timer.h"
#include "../../tools/config.h"

namespace fs = std::filesystem;

int main(int argc, char** argv){
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

    fs::path config_file;
    fs::directory_iterator files(folder);
    for (const fs::directory_entry& entry : files) {
        if (entry.path().filename().string() == "config.lua") {
            config_file = entry.path();
        }
    }

    if (verbose) {
        std::cout << "\nInput folder: " << folder << std::endl;
        std::cout << "Config file: " << config_file << std::endl;
    }

    std::string output_directory_str;
    fs::path output_directory;
    int number_of_threads;
    std::vector<int> energies, vib_modes, vib_degen;
    std::vector<double> ir_intens;
    lua_State* L = luaL_newstate();
    {
        std::cout << "\nLoading configuration..." << std::endl;
        Timer timer("Took");
        luaL_openlibs(L);

        fs::path cwd(fs::canonical("."));
        if (luaL_dostring(L, std::format("current_working_directory = \"{}\"", cwd.generic_string().c_str()).c_str()) != LUA_OK) {
            std::cout << "Error setting current working directory. Message:" << std::endl;
            std::cout << lua_tostring(L, -1) << std::endl;
            lua_close(L);
            return 1;
        }

        if (luaL_dostring(L, std::format("config_directory = \"{}\"", config_file.parent_path().generic_string().c_str()).c_str()) != LUA_OK) {
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

        lua_getglobal(L, "energies");
        if (Lua_Load_1d_Vector<int>(L, energies)) {
            std::cout << "Required parameter \"energies\" not found. Exiting." << std::endl;
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


        lua_getglobal(L, "number_of_threads");
        if (Lua_Load_Number<int>(L, &number_of_threads)) {
            std::cout << "Required parameter \"number_of_threads\" not found. Exiting." << std::endl;
            return 1;
        }
    }

    output_directory = fs::path(output_directory_str);
    if (!fs::is_directory(output_directory)) {
        std::cout << "Output directory does not exist. Creating at " << output_directory.string().c_str() << "\n" << std::endl;
        fs::create_directory(output_directory);
    }
    output_directory = fs::canonical(output_directory);

    Vibrational_Modes modes(vib_modes.data(), vib_degen.data(), ir_intens.data(), 298.15, vib_modes.size());

    if (verbose){
        modes.print();
    }

    ThreadPool pool;
    pool.Start(number_of_threads, [] {});
    {
        std::cout << "Computing occupation probabilities..." << std::endl;
        Timer timer("Took");
        Progress_Bar bar(energies.size(), 50, "");
        bar.Start();
        std::mutex bar_mutex;
        std::mutex busy_mutex;
        std::condition_variable busy_cv;
        for (int energy : energies) {
            pool.QueueJob([&modes, energy, output_directory, &bar, &bar_mutex] {
                    std::string out_file(std::format("{:d}.csv", energy));
                    Occupation occupations(&modes, energy);
                    Occupation_Density(&modes, &occupations, energy);
                    occupations.save(output_directory / out_file, &modes);
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
    pool.Stop();
    return 0;
}