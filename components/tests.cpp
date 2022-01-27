#include "tests.h"

#include <iostream>

#include "potentials.h"
#include "solvers.h"

namespace MD {

void lj_test(Eigen::Array<double, Eigen::Dynamic, 5> &data,
             uint num_x_samples) {
    //  constants from Rahman (1964)
    const double sigma = 3.4;  // angstrom
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

    std::cout << "[Info] Performed lj_test" << std::endl;
}

void vv_test(Eigen::ArrayXXd &data, double time_step, uint num_t_steps) {
    double sigma = 3.4;
    uint num_particles = 2;
    Eigen::ArrayX3d positions{{0, 0, 0}, {4 / sigma, 0, 0}};
    Eigen::ArrayX3d velocities{{0, 0, 0}, {0, 0, 0}};
    // Eigen::ArrayXXd data_vv = Eigen::ArrayXXd::Zero(num_t_steps, 7);

    MD::velocity_verlet(positions,
                        velocities,
                        MD::lennard_jones,
                        time_step,
                        num_t_steps,
                        num_particles,
                        MD::sample_x2,
                        data);

    std::cout << "[Info] Performed vv_test" << std::endl;
}

}  // namespace MD