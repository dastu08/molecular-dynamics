#include <Eigen/Core>
#include <iostream>

#include "components/helpers.h"
#include "components/initialization.h"
#include "components/potentials.h"
#include "components/solvers.h"
#include "components/tests.h"

// #define REPORT_1
#define REPORT_2
// #define MIC
#define SIM

// see Rahman for the constant
const double sigma = 3.4;                       // angstrom
const double temperature = 0.7916666666666669;  // 95 K in reduced units
const double separation = 1.067324157068294;    // 3.63 angstrom in reduced units
const uint seed = 8028;

int main() {
    // Report 1. Dynamics
#ifdef REPORT_1

    // Potential test
    uint num_samples = 1000;
    Eigen::Array<double, Eigen::Dynamic, 5> data_lj(num_samples, 5);
    MD::lj_test(data_lj, num_samples);
    MD::array2file(data_lj, "../data/01/lj_test.txt", "x,energy,fx,fy,fz");

    // Velocity Verlet test
    double time_step = 0.01;  // accurate
    uint num_t_steps = 2000;
    Eigen::ArrayXXd data_vv = Eigen::ArrayXXd::Zero(num_t_steps, 7);
    MD::vv_test(data_vv, time_step, num_t_steps);
    MD::array2file(data_vv, "../data/01/vv_test.txt", "t,x1,v1,x2,v2,epot,ekin");

    time_step = 0.05;  // inaccurate
    data_vv = Eigen::ArrayXXd::Zero(num_t_steps, 7);
    MD::vv_test(data_vv, time_step, num_t_steps);
    MD::array2file(data_vv, "../data/01/vv_test2.txt", "t,x1,v1,x2,v2,epot,ekin");

    // Scan for the stability limit
    uint time_step_samples = 100;
    Eigen::ArrayX2d data_ts = Eigen::ArrayXXd::Zero(time_step_samples, 2);
    MD::time_step_test(data_ts, 0.01, 0.06, time_step_samples);
    MD::array2file(data_ts, "../data/01/ts_test_index.txt", "index,ts");
#endif  // REPORT_1

// Report 2. Extended system
#ifdef REPORT_2
    Eigen::ArrayX3d positions, positions_wrapped, velocities, forces;
    Eigen::ArrayXXd energies, equilib;
    double box_length, time_step;
    uint n, num_particles, num_t_steps, equilib_steps;

#ifdef MIC
    // Potential test 2 with minimal image convention
    uint num_samples = 1000;
    Eigen::Array<double, Eigen::Dynamic, 5> data_lj(num_samples, 5);
    MD::lj_test2(data_lj, num_samples);
    MD::array2file(data_lj, "../data/02/lj_test2.txt", "x,energy,fx,fy,fz");

    // minimal image convention test
    time_step = 0.01;  // accurate
    num_t_steps = 200;
    Eigen::ArrayXXd data_mic = Eigen::ArrayXXd::Zero(num_t_steps, 7);
    MD::mic_test(data_mic, time_step, num_t_steps, 3);
    MD::array2file(data_mic, "../data/02/mic_test.txt", "t,x1,v1,x2,v2,epot,ekin");
#endif

#ifdef SIM
    // // spacial position init
    n = 8;
    time_step = 0.005;
    num_t_steps = 5000;
    equilib_steps = 500;
    energies = Eigen::ArrayXXd::Zero(num_t_steps, 3);
    equilib = Eigen::ArrayXXd::Zero(equilib_steps, 3);

    num_particles = MD::init_positions_3d(positions, n, separation);
    box_length = n * separation;
    std::cout << "box length: " << box_length << '\n'
              << "separation: " << separation << '\n'
              << "number of particles: " << num_particles
              << std::endl;

    MD::init_velocities_3d(velocities, num_particles, seed, 3);
    MD::velocity_drift_removal(velocities);
    MD::velocity_rescaling(velocities, num_particles, temperature);

    MD::equilibration_phase(positions, velocities, time_step, equilib_steps,
                            num_particles, temperature, box_length,
                            "../data/02/equilibration1.txt");
    MD::equilibration_phase(positions, velocities, time_step, equilib_steps,
                            num_particles, temperature, box_length,
                            "../data/02/equilibration2.txt");
    MD::equilibration_phase(positions, velocities, time_step, equilib_steps,
                            num_particles, temperature, box_length,
                            "../data/02/equilibration3.txt");

    // sampling
    MD::velocity_verlet(positions,
                        velocities,
                        MD::lennard_jones,
                        time_step,
                        num_t_steps,
                        num_particles,
                        MD::sample_energies,
                        energies,
                        box_length);
    std::cout << "sampling of duration " << time_step * num_t_steps
              << ", final temperature: "
              << MD::computeTemperature(velocities) << std::endl;

    positions_wrapped = positions;
    MD::coordinate_wrapping(positions_wrapped, num_particles, box_length);
    MD::array2file(energies, "../data/02/energies.txt", "t,epot,ekin");
    MD::array2file(positions, "../data/02/positions.txt", "x,y,z");
    MD::array2file(positions_wrapped,
                   "../data/02/positions_wrapped.txt", "x,y,z");
    MD::array2file(velocities, "../data/02/velocities.txt", "vx,vy,vz");
#endif
#endif  // REPORT_2

    return 0;
}