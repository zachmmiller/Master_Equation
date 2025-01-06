//
// Created by Zach Miller on 11/30/23.
//

#include "Argument_Parsing.h"

int Parse_Arguments(int argc, char **argv, std::vector<arg_value_pair> *out, const std::vector<arg_type> &required_args, const std::vector<arg_type> &optional_args) {
    const size_t optional_args_size = optional_args.size();
    const size_t required_args_size = required_args.size();
    bool *required_args_found = new bool[required_args_size];
    std::fill(required_args_found, required_args_found + required_args_size, false);

    std::string arg_name;
    std::string arg_value;
    bool arg_valid = false;
    for (int i = 1; i < argc; i++) {
        arg_valid = false;
        arg_name = std::string(argv[i]);

        for (int j = 0; j < required_args_size; j++) {
            if (required_args[j].name == arg_name) {
                if (required_args[j].is_flag) {
                    required_args_found[j] = true;
                    arg_valid = true;
                    out->push_back({arg_name, ""});
                } else {
                    i++;
                    arg_value = std::string(argv[i]);
                    required_args_found[j] = true;
                    arg_valid = true;
                    out->push_back({arg_name, arg_value});
                }
            }
        }

        if (arg_valid) {
            continue;
        }

        for (int j = 0; j < optional_args_size; j++) {
            if (optional_args[j].name == arg_name) {
                if (optional_args[j].is_flag) {
                    required_args_found[j] = true;
                    arg_valid = true;
                    out->push_back({arg_name, ""});
                } else {
                    i++;
                    arg_value = std::string(argv[i]);
                    required_args_found[j] = true;
                    arg_valid = true;
                    out->push_back({arg_name, arg_value});
                }
            }
        }
        if (!arg_valid) {
            std::cout << arg_name << " is an invalid argument" << std::endl;
            delete[] required_args_found;
            return 1;
        }
    }
    bool required_argument_remaining = false;
    for (int i = 0; i < required_args_size; i++) {
        if (!required_args_found[i]) {
            std::cout << required_args[i].name << " argument not found" << std::endl;
            required_argument_remaining = true;
        }
    }
    delete[] required_args_found;
    if (required_argument_remaining) {
        return 1;
    } else {
        return 0;
    }
}
