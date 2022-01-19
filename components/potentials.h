#ifndef __POTENTIALS_H__
#define __POTENTIALS_H__

#include <Eigen/Core>

namespace MD {

// Compute the LJ energy and forces for N particles in 3 dimensions.
double lennard_jones(Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &forces,
                     uint num_particles);

}  // namespace MD

#endif  // __POTENTIALS_H__