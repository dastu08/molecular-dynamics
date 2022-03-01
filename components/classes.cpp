#include "classes.h"

#include <iostream>

#include "helpers.h"
#include "initialization.h"
#include "potentials.h"
#include "samplers.h"
#include "solvers.h"

namespace MD {

Simulation::Simulation(uint n,
                       double time_step,
                       uint num_t_steps,
                       double box_length,
                       uint num_rescales,
                       uint num_rescales_steps,
                       uint num_bins) {
    this->n = n;
    this->time_step = time_step;
    this->num_t_steps = num_t_steps;
    this->box_length = box_length;
    this->num_rescales = num_rescales;
    this->num_rescale_steps = num_rescales_steps;
    this->num_bins = num_bins;

    initialzedFlag = 0;
}

void Simulation::printInfo() {
    std::cout << "[Info] MD Simulation Information\n"
              << "# particles \t\t: " << num_particles << '\n'
              << "time step: \t\t: " << time_step << '\n'
              << "# time steps \t\t: " << num_t_steps << '\n'
              << "separation \t\t: " << separation << '\n'
              << "box length \t\t: " << box_length << '\n'
              << "# velocity rescales \t: " << num_rescales << '\n'
              << "steps btwn rescales \t: " << num_rescale_steps << '\n'
              << "# bins \t\t\t: " << num_bins << '\n'
              << "initializedFlag \t: " << initialzedFlag
              << std::endl;
}

void Simulation::init(double separation, double temperature, uint seed) {
    // set size for data array
    data = Eigen::ArrayXXd::Zero(num_t_steps + num_rescale_steps * num_rescales,
                                 3 + num_bins);

    //  initi positions and velocities
    num_particles = MD::init_positions_3d(positions, n, separation);
    MD::init_velocities_3d(velocities, num_particles, seed, 3);
    MD::velocity_drift_removal(velocities);

    // set flag to ready
    initialzedFlag = 1;
}

void Simulation::run() {
    if (initialzedFlag != 1) {
        std::cout << "[Warning] The simulation was not initialized. "
                  << "Skipping run. (initializedFlag = "
                  << initialzedFlag << ")"
                  << std::endl;
        return;
    }
    // rescaling
    for (uint i = 0; i < num_rescales; ++i) {
        MD::velocity_rescaling(velocities, num_particles, temperature);
        MD::velocity_verlet(positions,
                            velocities,
                            MD::lennard_jones,
                            time_step,
                            num_rescale_steps,
                            num_particles,
                            MD::sample_energies_rdf,
                            data,
                            box_length,
                            num_bins,
                            i * num_rescale_steps);
    }

    // sampling
    MD::velocity_verlet(positions,
                        velocities,
                        MD::lennard_jones,
                        time_step,
                        num_t_steps,
                        num_particles,
                        MD::sample_energies_rdf,
                        data,
                        box_length,
                        num_bins,
                        num_rescale_steps * num_rescales);

    std::cout << "[Info] Finished simulation of duration "
              << time_step * num_t_steps
              << " at final temperature: "
              << MD::computeTemperature(velocities) << std::endl;
}

void Simulation::export2file(std::string filepath) {
    std::string head;
    // create the header string for the bins
    MD::stringList(head, "n", num_bins);

    // save the sampled data (t, epot, ekin, n0, ...) to a file
    MD::array2file(data,
                   filepath + std::to_string(n) + ".txt",
                   "t,epot,ekin," + head);
}

}  // namespace MD

namespace MC {

Simulation::Simulation(uint n,
                       uint num_samples,
                       double box_length,
                       double temperature,
                       uint num_bins) {
    this->n = n;
    this->num_samples = num_samples;
    this->box_length = box_length;
    this->temperature = temperature;
    this->num_bins = num_bins;

    initialzedFlag = 0;
}

void Simulation::printInfo() {
    std::cout << "[Info] MC Simulation Information\n"
              << "# particles \t\t: " << num_particles << '\n'
              << "# samplings \t\t: " << num_samples << '\n'
              << "separation \t\t: " << separation << '\n'
              << "box length \t\t: " << box_length << '\n'
              << "# bins \t\t\t: " << num_bins << '\n'
              << "initializedFlag \t: " << initialzedFlag
              << std::endl;
}

void Simulation::init(double separation, uint seed) {
    // set size for data array
    // data = Eigen::ArrayXXd::Zero(num_samples, 3 + num_bins);
    data = Eigen::ArrayXXd::Zero(num_samples, 4);

    //  initi positions and velocities
    num_particles = MD::init_positions_3d(positions, n, separation);

    this->seed = seed;

    // set flag to ready
    initialzedFlag = 1;
}

void Simulation::run() {
    if (initialzedFlag != 1) {
        std::cout << "[Warning] The simulation was not initialized. "
                  << "Skipping run. (initializedFlag = "
                  << initialzedFlag << ")"
                  << std::endl;
        return;
    }

    MC::metropolis(positions, num_samples, num_particles, seed,
                   nullptr, data, box_length, num_bins);

    MD::array2file(data, "../data/05/hist.txt", "idx,x1,x2,x3");
    // std::cout << "[Info] Finished simulation of duration "
    //           << time_step * num_t_steps
    //           << " at final temperature: "
    //           << MD::computeTemperature(velocities) << std::endl;
}

void Simulation::export2file(std::string filepath) {
    std::string head;
    // create the header string for the bins
    MD::stringList(head, "n", num_bins);

    // save the sampled data (t, epot, ekin, n0, ...) to a file
    MD::array2file(data,
                   filepath + std::to_string(n) + ".txt",
                   "t,epot,ekin," + head);
}

}  // namespace MC
