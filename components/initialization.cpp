#include "initialization.h"

#include <stdlib.h>

#include <iostream>

namespace MD {

uint init_potisionts_3d(Eigen::ArrayX3d &positions,
                        uint num_paricles_cubic,
                        double distance) {
    uint num_particles = num_paricles_cubic * num_paricles_cubic * num_paricles_cubic;
    uint index;
    positions = Eigen::ArrayX3d::Zero(num_particles, 3);

    for (uint ix = 0; ix < num_paricles_cubic; ++ix) {
        for (uint iy = 0; iy < num_paricles_cubic; ++iy) {
            for (uint iz = 0; iz < num_paricles_cubic; ++iz) {
                index = (ix * num_paricles_cubic + iy) * num_paricles_cubic + iz;
                positions(index, 0) = ix * distance;
                positions(index, 1) = iy * distance;
                positions(index, 2) = iz * distance;
            }
        }
    }

    return num_particles;
}

void init_velocities_3d(Eigen::ArrayX3d &velocities,
                        uint num_particles,
                        uint seed,
                        double max_veloctiy) {
    // init random number generator
    srand(seed);
    for (uint i = 0; i < num_particles; ++i) {
        for (uint j = 0; j < 3; j++) {
            velocities(i, j) = ((double)rand() / RAND_MAX - 0.5) * 2 * max_veloctiy;
        }
    }
}

void velocity_drift_removal(Eigen::ArrayX3d &velocities) {
    // center of mass velocity
    Eigen::Vector3d v = velocities.colwise().sum() / velocities.rows();
    // velocities relative to the center of mass
    velocities.rowwise() -= v.transpose().array();
}

void velocity_rescaling(Eigen::ArrayX3d &velocities,
                        uint num_particles,
                        double temp_target) {
    double temp_now = velocities.square().sum() / 3 / num_particles;

    velocities *= sqrt(temp_target / temp_now);

    std::cout << "temp before: " << temp_now
              << "\ntemp after: " << velocities.square().sum() / 3 / num_particles
              << std::endl;
}

}  // namespace MD
