#include <Eigen/Core>
#include <iostream>

#include "components/helpers.h"
#include "components/potentials.h"
#include "components/solvers.h"

// #define REPORT_1_1
#define REPORT_1_2
// #define FORCE_TEST

// see Rahman for the constant
const double sigma = 3.4;

int main() {
    // Report 1. Dynamics
#ifdef REPORT_1_1
    uint num_samples = 1000;
    Eigen::Array<double, Eigen::Dynamic, 5> data_lj(num_samples, 5);
    MD::lj_test(data_lj, num_samples);
    MD::array2file(data_lj, "../reports/lj_test.txt", "x,energy,fx,fy,fz");
#endif  // REPORT_1_1

#ifdef FORCES_TEST
    // Example for plotting the forces of some particles
    uint num_particles = 5;
    Eigen::ArrayX3d positions{{0, 0, 0},
                              {2, 0, 0},
                              {2, 1.5, 0},
                              {0.8, 2, 0},
                              {0.4, 3, 0}};
    Eigen::ArrayX3d forces;
    double epot = MD::lennard_jones(positions, forces, num_particles);
    MD::array2file((Eigen::ArrayXXd(num_particles, 6) << positions, forces).finished(),
                   "../reports/forces_3_particles.txt", "x,y,z,fx,fy,fz");
#endif  // FORCES_TEST

#ifdef REPORT_1_2
    double time_step = 0.01;
    uint num_t_steps = 800;
    uint num_particles = 2;
    Eigen::ArrayX3d positions{{0, 0, 0}, {4 / sigma, 0, 0}};
    Eigen::ArrayX3d velocities{{0, 0, 0}, {0, 0, 0}};
    Eigen::ArrayXXd data_vv = Eigen::ArrayXXd::Zero(num_t_steps, 7);

    MD::velocity_verlet(positions,
                        velocities,
                        MD::lennard_jones,
                        time_step,
                        num_t_steps,
                        num_particles,
                        MD::sample_x2,
                        data_vv);

    MD::array2file(data_vv, "../reports/vv_test.txt", "t,x1,v1,x2,v2,epot,ekin");

    // std::cout << positions << std::endl;

    // MD::vv_test();
#endif  // REPORT_1_2

    return 0;
}