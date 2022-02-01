#include "tests.h"

#include <iomanip>
#include <iostream>
#include <sstream>

#include "helpers.h"
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
    Eigen::ArrayX3d positions{{0, 0, 0}, {5 / sigma, 0, 0}};
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

void time_step_test(Eigen::ArrayX2d &data,
                    double time_step_0,
                    double time_step_1,
                    uint num_t_step_samples) {
    double sigma = 3.4;
    double time_step_delta = (time_step_1 - time_step_0) / num_t_step_samples;
    double time_step;
    uint num_particles = 2;
    uint num_t_steps = 2000;
    Eigen::ArrayX3d positions{{0, 0, 0}, {5 / sigma, 0, 0}};
    Eigen::ArrayX3d velocities{{0, 0, 0}, {0, 0, 0}};
    Eigen::ArrayXXd data_vv = Eigen::ArrayXXd::Zero(num_t_steps, 3);
    std::ostringstream oss;

    for (uint i = 0; i < num_t_step_samples; i++) {
        time_step = time_step_0 + i * time_step_delta;
        data(i, 0) = i;
        data(i, 1) = time_step;

        MD::velocity_verlet(positions,
                            velocities,
                            MD::lennard_jones,
                            time_step,
                            num_t_steps,
                            num_particles,
                            MD::sample_energies,
                            data_vv);

        oss << "../data/01/ts_test"
            << std::setfill('0')
            << std::setw(3)
            << i << ".txt";
        MD::array2file(data_vv, oss.str(), "t,epot,ekin");
        oss.clear();
        oss.str("");
    }

    std::cout << "[Info] Performed time_step_test" << std::endl;
}

}  // namespace MD