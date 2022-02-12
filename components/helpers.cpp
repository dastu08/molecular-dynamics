#include "helpers.h"

#include <fstream>
#include <iostream>

namespace MD {

void array2file(Eigen::ArrayXXd array, std::string filename, std::string head) {
    std::ofstream file(filename);

    if (!file.fail()) {
        // write head line
        file << head << '\n';
        // presicion for printing floating point numbers
        file << std::scientific;
        file.precision(16);
        // print row wise
        for (auto row : array.rowwise()) {
            for (int i = 0; i < row.size() - 1; i++) {
                file << row(i) << ',';
            }
            file << row(Eigen::last) << "\n";
        }

        std::cout << "[Info] Wrote array to file: " << filename << std::endl;
    } else {
        std::cout << "[Warning] Could not open file: " << filename << std::endl;
    }

    file.close();
}

void stringList(std::string &string, std::string symbol, uint repetition) {
    // clear string first
    string.empty();
    // append counted symbol
    for (uint i = 0; i < repetition; ++i) {
        string.append(symbol + std::to_string(i) + ',');
    }
    string.erase(string.end()-1);
}

}  // namespace MD