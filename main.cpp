#include <Eigen/Core>
#include <iostream>

#include "components/helpers.h"
#include "components/potentials.h"
#include "components/solvers.h"

// #define REPORT_1_1
#define REPORT_1_2

// see Rahman for the constant
const double sigma = 3.4;

int main() {
    // report 1. dynamics
#ifdef REPORT_1_1
    uint num_samples = 400;
    Eigen::Array<double, Eigen::Dynamic, 5> data_lj(num_samples, 5);
    MD::lj_test(data_lj, num_samples);
    MD::array2file(data_lj, "../reports/lj_test.txt", "x,energy,fx,fy,fz");
#endif

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
                        MD::sample_x,
                        data_vv);

    MD::array2file(data_vv, "../reports/vv_test.txt", "t,x1,v1,x2,v2,epot,ekin");

    // std::cout << positions << std::endl;

    // MD::vv_test();
#endif

    return 0;
}