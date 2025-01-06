#ifndef EXPORTING_H
#define EXPORTING_H

#include <filesystem>

namespace fs = std::filesystem;

void export_vector(double* vec, size_t N, fs::path file);
void export_matrix(double* mat, size_t dim, fs::path file);

#endif
