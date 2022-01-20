// define functions that compute the potentials

#include "potentials.h"

namespace MD {

double lennard_jones(const Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &forces,
                     uint num_particles) {
    double energy = 0;
    // reset the forces
    forces = Eigen::ArrayX3d::Zero(num_particles, 3);

    // loop over i-j pairs of particles.
    for (int i = 0; i < num_particles; i++) {
        // since i-j  and j-i are equivalent only start from j = i + 1
        for (int j = i + 1; j < num_particles; j++) {
            // compute the force and energy
        }
    }

    return energy;
}

}  // namespace MD