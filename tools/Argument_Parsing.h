//
// Created by Zach Miller on 11/30/23.
//

#ifndef CDMS_ANALYSIS_ARGUMENT_PARSING_H
#define CDMS_ANALYSIS_ARGUMENT_PARSING_H

#include <iostream>
#include <string>
#include <vector>

typedef struct {
    std::string name;
    bool is_flag;
} arg_type;

typedef struct {
    std::string name;
    std::string value;
} arg_value_pair;

int Parse_Arguments(int argc, char **argv, std::vector<arg_value_pair> *out, const std::vector<arg_type> &required_args, const std::vector<arg_type> &optional_args);

#endif  // CDMS_ANALYSIS_ARGUMENT_PARSING_H
