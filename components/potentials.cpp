// define functions that compute the potentials

#include "potentials.h"

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
    for (int i = 0; i < num_particles; i++) {
        // since i-j  and j-i are equivalent only start from j = i + 1
        for (int j = i + 1; j < num_particles; j++) {
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
            energy_total += 4 * r6_inv * (r6_inv - 1);
        }
    }

    return energy_total;
}

void lj_test(Eigen::Array<double, Eigen::Dynamic, 5> &data,
             uint num_x_samples) {

    //  constants from Rahman (1964)
    const double sigma = 3.4; // angstrom
    const double x_low = 3 / sigma;
    const double x_high = 10 / sigma;
    const uint num_particles = 2;
    // initialize the positions
    Eigen::ArrayX3d positions = Eigen::ArrayX3d::Zero(num_particles, 3);
    Eigen::ArrayX3d forces;
    double energy;
    uint iter = 0;

    Eigen::VectorXd xs = Eigen::VectorXd::LinSpaced(num_x_samples, x_low, x_high);

    for (double x : xs) {
        // set only x-coordinate of particle 2
        positions(1, 0) = x;

        // compute the LJ energy
        data(iter, 0) = x;
        data(iter, 1) = MD::lennard_jones(positions, forces, num_particles);
        data(iter, {2, 3, 4}) = forces.row(0).transpose();
        ++iter;
    }
}

}  // namespace MD