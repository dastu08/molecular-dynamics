#include "initialization.h"

namespace MD {

void init3d(Eigen::ArrayX3d &positions,
            uint num_paricles_cubic,
            double distance) {
    uint num_particles = num_paricles_cubic * num_paricles_cubic * num_paricles_cubic;
    uint index;
    positions = Eigen::ArrayX3d::Zero(num_particles, 3);

    for (uint ix = 0; ix < num_paricles_cubic; ix++) {
        for (uint iy = 0; iy < num_paricles_cubic; iy++) {
            for (uint iz = 0; iz < num_paricles_cubic; iz++) {
                index = (ix * num_paricles_cubic + iy) * num_paricles_cubic + iz;
                positions(index, 0) = ix * distance;
                positions(index, 1) = iy * distance;
                positions(index, 2) = iz * distance;
            }
        }
    }
}

}  // namespace MD
