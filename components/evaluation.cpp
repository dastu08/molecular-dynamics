#include "evaluation.h"

#include <iostream>

namespace MD {

void radial_distribution_hist(Eigen::ArrayX3d& positions,
                              Eigen::VectorXi& r_hist,
                              uint num_bins,
                              double box_length) {
    double delta_r = box_length / num_bins;
    uint num_particles = positions.rows();
    int idx;
    Eigen::Vector3d xij;

    r_hist = Eigen::VectorXi::Zero(num_bins);

    // std::cout << "delta_r = " << delta_r << std::endl;

    // loop over i-j pairs of particles.
    for (uint i = 0; i < num_particles; i++) {
        // since i-j  and j-i are equivalent only start from j = i + 1
        for (uint j = i + 1; j < num_particles; j++) {
            // relative position btwn the particles
            xij = positions.row(i) - positions.row(j);
            // minimum image convention
            // reduce the components to fit inside the box
            for (uint k = 0; k < 3; ++k) {
                xij(k) -= box_length * round(xij(k) / box_length);
            }

            idx = floor(xij.norm() / delta_r);
            if (idx > (int)num_bins) {
                std::cout << "[Error] The computed distance is outside the bin range! Ignoring."
                          << std::endl;
            } else {
                r_hist(idx) += 1;
            }
        }
    }
}

}  // namespace MD