#include <Eigen/Core>
#include <iostream>

#include "components/helpers.h"
#include "components/potentials.h"
#include "components/solvers.h"
#include "components/tests.h"

// #define REPORT_1_1
// #define REPORT_1_2
#define MANY_PARTICLES
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

#ifdef REPORT_1_2
    double time_step = 0.01;
    uint num_t_steps = 800;
    Eigen::ArrayXXd data_vv = Eigen::ArrayXXd::Zero(num_t_steps, 7);
    MD::vv_test(data_vv, time_step, num_t_steps);
    MD::array2file(data_vv, "../reports/vv_test.txt", "t,x1,v1,x2,v2,epot,ekin");
#endif  // REPORT_1_2

#ifdef MANY_PARTICLES
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

    return 0;
}