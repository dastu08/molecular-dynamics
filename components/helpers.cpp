#include "helpers.h"

#include <fstream>
#include <iostream>

namespace MD {

void array2file(Eigen::ArrayX3d array, std::string filename, std::string head) {
    std::ofstream file(filename);

    if (!file.fail()) {
        // write head line
        file << head << '\n';
        // presicion for printing floating point numbers
        file << std::scientific;
        file.precision(16);
        // print row wise
        for (auto row : array.rowwise()) {
            file << row(0) << ',' << row(1) << ',' << row(2) << '\n';
        }

        std::cout << "[Info] Wrote array to file: " << filename << std::endl;
    } else {
        std::cout << "[Warning] Could not open file: " << filename << std::endl;
    }

    file.close();
}

}  // namespace MD