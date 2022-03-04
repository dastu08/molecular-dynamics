// define functions that compute the potentials

#include "potentials.h"

#include <iostream>

namespace MD {

double lennard_jones(const Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &forces,
                     uint num_particles) {
    double energy_total = 0;
    // auxiliary variables
    double r2_inv, r6_inv, fij;
    Eigen::Vector3d xij;

    // reset the forces
    forces = Eigen::ArrayX3d::Zero(num_particles, 3);

    // loop over i-j pairs of particles.
    for (uint i = 0; i < num_particles; i++) {
        // since i-j  and j-i are equivalent only start from j = i + 1
        for (uint j = i + 1; j < num_particles; j++) {
            // relative position btwn the particles
            xij = positions.row(i) - positions.row(j);
            // compute powers of the inverse distance
            r2_inv = 1 / (xij.squaredNorm());
            r6_inv = r2_inv * r2_inv * r2_inv;

            // F_i = - V'(r)/r (x_i - x_j) = -F_j
            // = 48 / r^8 (1/r^6 - 1/2) (x_i - x_j)
            fij = 48 * r2_inv * r6_inv * (r6_inv - 0.5);
            forces.row(i) += fij * xij.array();
            forces.row(j) -= fij * xij.array();

            // add the potential energy of the pair to the total energy
            // only once because it already incorporates the pair of particles
            energy_total += 4 * r6_inv * (r6_inv - 1);
        }
    }

    return energy_total;
}

double lennard_jones(const Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &forces,
                     uint num_particles,
                     double side_length) {
    double energy_total = 0;
    // auxiliary variables
    double r2_inv, r6_inv, fij;
    Eigen::Vector3d xij;

    // reset the forces
    forces = Eigen::ArrayX3d::Zero(num_particles, 3);

    // loop over i-j pairs of particles.
    for (uint i = 0; i < num_particles; i++) {
        // since i-j  and j-i are equivalent only start from j = i + 1
        for (uint j = i + 1; j < num_particles; j++) {
            // relative position btwn the particles
            xij = positions.row(i) - positions.row(j);
            // minimum image convention
            // reduce the components to fit inside the box
            for (uint k = 0; k < 3; ++k) {
                xij(k) -= side_length * round(xij(k) / side_length);
            }

            // compute powers of the inverse distance
            r2_inv = 1 / (xij.squaredNorm());
            r6_inv = r2_inv * r2_inv * r2_inv;

            // F_i = - V'(r)/r (x_i - x_j) = -F_j
            // = 48 / r^8 (1/r^6 - 1/2) (x_i - x_j)
            fij = 48 * r2_inv * r6_inv * (r6_inv - 0.5);
            forces.row(i) += fij * xij.array();
            forces.row(j) -= fij * xij.array();

            // add the potential energy of the pair to the total energy
            // only once because it already incorporates the pair of particles
            energy_total += 4 * r6_inv * (r6_inv - 1);
        }
    }

    return energy_total;
}

}  // namespace MD

namespace MC {

double lennard_jones_rdf(const Eigen::ArrayX3d &positions,
                         uint num_particles,
                         double side_length,
                         Eigen::VectorXi &r_hist,
                         double num_bins) {
    double energy_total = 0;
    // auxiliary variables
    double r2_inv, r6_inv;
    Eigen::Vector3d xij;
    // radial distribution needs this
    double delta_r = side_length / num_bins;
    int idx;
    r_hist = Eigen::VectorXi::Zero(num_bins);

    // loop over i-j pairs of particles.
    for (uint i = 0; i < num_particles; i++) {
        // since i-j  and j-i are equivalent only start from j = i + 1
        for (uint j = i + 1; j < num_particles; j++) {
            // relative position btwn the particles
            xij = positions.row(i) - positions.row(j);
            // minimum image convention
            // reduce the components to fit inside the box
            for (uint k = 0; k < 3; ++k) {
                xij(k) -= side_length * round(xij(k) / side_length);
            }

            // compute powers of the inverse distance
            r2_inv = 1 / (xij.squaredNorm());
            r6_inv = r2_inv * r2_inv * r2_inv;

            // add the potential energy of the pair to the total energy
            // only once because it already incorporates the pair of particles
            energy_total += 4 * r6_inv * (r6_inv - 1);

            // radial distribution calcualtion
            idx = floor(xij.norm() / delta_r);
            // check if the bin index is valid
            if (idx > (int)num_bins) {
                std::cout << "[Error] The computed distance is outside the bin range! Ignoring."
                          << std::endl;
            } else {
                // increase the count of the radial bin
                r_hist(idx) += 1;
            }
        }
    }

    return energy_total;
}

double lennard_jones_single(const Eigen::ArrayX3d &positions,
                            uint num_particles,
                            double side_length,
                            uint particle) {
    double energy_total = 0;
    // auxiliary variables
    double r2_inv, r6_inv;
    Eigen::Vector3d xij;

    // loop over the other particles
    for (uint i = 0; i < num_particles; i++) {
        // exclude self interaction
        if (i != particle) {
            // relative position btwn the particles
            xij = positions.row(particle) - positions.row(i);
            // minimum image convention
            // reduce the components to fit inside the box
            for (uint k = 0; k < 3; ++k) {
                xij(k) -= side_length * round(xij(k) / side_length);
            }

            // compute powers of the inverse distance
            r2_inv = 1 / (xij.squaredNorm());
            r6_inv = r2_inv * r2_inv * r2_inv;

            // add the potential energy of the pair to the total energy
            // only once because it already incorporates the pair of particles
            energy_total += 4 * r6_inv * (r6_inv - 1);
        }
    }

    return energy_total;
}

}  // namespace MC
