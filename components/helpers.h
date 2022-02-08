#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <Eigen/Core>
#include <string>

namespace MD {

void array2file(Eigen::ArrayXXd array, std::string filename, std::string head);

inline double computeTemperature(Eigen::ArrayX3d& velocities) {
    return velocities.square().sum() / velocities.rows() / 3;
}

}  // namespace MD

#endif  // __HELPERS_H__