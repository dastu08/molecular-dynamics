#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <Eigen/Core>
#include <string>

namespace MD {

void array2file(Eigen::ArrayXXd array, std::string filename, std::string head);

}  // namespace MD

#endif  // __HELPERS_H__