#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <Eigen/Core>
#include <string>

namespace MD {

/* Write the array to a file

Parameters:
    * array:
        The array containing the data which is going to be written to a file
    * filename:
        string of the filename, it can even contain an absolute/relative path
    * head:
        first line of the file, containing header information

Description:
    * Create or override the file from filename.
    * Write the array in it's dimensions (rows, columns) to the file
    with precision 16 decimal digits and scientific format.
*/
void array2file(Eigen::ArrayXXd array, std::string filename, std::string head);

// Compute T = 1/(3N) sum_{i=1}^{N} v_i^2
inline double computeTemperature(Eigen::ArrayX3d& velocities) {
    return velocities.square().sum() / velocities.rows() / 3;
}

/* Repeat the symbol as a comma separated list

Parameters:
    * string:
        reference to the output string
    * symbol:
        the symbol that will be repeated
    * repetition:
        number of repetitions of the symbol

Description:
    * Clear the string.
    * Repeat the symbol appened by a number starting from 0 to repetitions-1.
    * The repetitions are comma seperated but there is no comma at the end.
*/
void stringList(std::string& string, std::string symbol, uint repetition);

}  // namespace MD

#endif  // __HELPERS_H__