#include <Eigen/Core>
#include <iostream>

#include "components/helpers.h"
#include "components/potentials.h"
#include "components/solvers.h"
#include "components/tests.h"

#define REPORT_1_1
#define REPORT_1_2
#define TIME_STEP_TEST
// #define MANY_PARTICLES
// see Rahman for the constant
const double sigma = 3.4;

void sample_xy5(uint index,
                double time,
                Eigen::ArrayX3d &positions,
                Eigen::ArrayX3d &velocities,
                double e_pot,
                Eigen::ArrayXXd &data) {
    // if ((index % 2) == 1) {
    //     return;
    // }

    data(index, 0) = time;
    data(index, 1) = e_pot;
    data(index, 2) = velocities.square().sum() / 2;
    // x position
    for (int i = 0; i < 5; i++) {
        data(index, 3 + i) = positions(i, 0);
        data(index, 8 + i) = positions(i, 1);
        // data(index, 13 + i) = velocities(i, 0);
        // data(index, 18 + i) = velocities(i, 1);
    }
}

int main() {
    // Report 1. Dynamics
#ifdef REPORT_1_1
    uint num_samples = 1000;
    Eigen::Array<double, Eigen::Dynamic, 5> data_lj(num_samples, 5);
    MD::lj_test(data_lj, num_samples);
    MD::array2file(data_lj, "../data/01/lj_test.txt", "x,energy,fx,fy,fz");
#endif  // REPORT_1_1

#ifdef REPORT_1_2
    double time_step = 0.01; // accurate
    uint num_t_steps = 2000;
    Eigen::ArrayXXd data_vv = Eigen::ArrayXXd::Zero(num_t_steps, 7);
    MD::vv_test(data_vv, time_step, num_t_steps);
    MD::array2file(data_vv, "../data/01/vv_test.txt", "t,x1,v1,x2,v2,epot,ekin");

    time_step = 0.05; // inaccurate
    data_vv = Eigen::ArrayXXd::Zero(num_t_steps, 7);
    MD::vv_test(data_vv, time_step, num_t_steps);
    MD::array2file(data_vv, "../data/01/vv_test2.txt", "t,x1,v1,x2,v2,epot,ekin");
#endif  // REPORT_1_2

#ifdef TIME_STEP_TEST
    uint time_step_samples = 100;
    Eigen::ArrayX2d data_ts = Eigen::ArrayXXd::Zero(time_step_samples, 2);
    MD::time_step_test(data_ts, 0.01, 0.06, time_step_samples);
    MD::array2file(data_ts, "../data/01/ts_test_index.txt", "index,ts");
#endif  // TIME_STEP_TEST

#ifdef MANY_PARTICLES
    // Example for plotting the forces of some particles
    uint num_particles = 5;
    uint num_t_steps = 1280000;
    double time_step = 0.0001;
    Eigen::ArrayX3d positions{{0, 0, 0},
                              {2, 0, 0},
                              {2, 1.5, 0},
                              {0.8, 2, 0},
                              {0.4, 3, 0}};
    Eigen::ArrayX3d forces;
    Eigen::ArrayX3d velocities = Eigen::ArrayX3d::Zero(num_particles, 3);
    Eigen::ArrayXXd data(num_t_steps, 3 + 2 * num_particles);

    // compute the static forces
    MD::lennard_jones(positions, forces, num_particles);
    MD::array2file((Eigen::ArrayXXd(num_particles, 6) << positions, forces).finished(),
                   "../reports/many_particles_forces.txt", "x,y,z,fx,fy,fz");

    // simulate the system
    MD::velocity_verlet(positions, velocities, MD::lennard_jones, time_step,
                        num_t_steps, num_particles, sample_xy5, data);
    MD::array2file(data, "../reports/many_particles_evolution.txt",
                   "t,epot,ekin,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5");  //,vx1,vx2,vx3,vx4,vx5,vy1,vy2,vy3,vy4,vy5");
#endif                                                            // FORCES_TEST

    return 0;
}