#include "solvers.h"

#include <iostream>

#include "initialization.h"

namespace MD {

void velocity_verlet(Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double (*force)(const Eigen::ArrayX3d &,
                                     Eigen::ArrayX3d &,
                                     uint),
                     double time_step,
                     uint num_t_steps,
                     uint num_particles,
                     void (*sampler)(uint,
                                     double,
                                     Eigen::ArrayX3d &,
                                     Eigen::ArrayX3d &,
                                     double,
                                     Eigen::ArrayXXd &),
                     Eigen::ArrayXXd &data) {
    Eigen::ArrayX3d forces;
    double e_pot;
    double t = 0;

    std::cout << "[Debug] Velocity verlet algorithm with "
              << num_t_steps << " steps."
              << std::endl;

    if (positions.size() != velocities.size()) {
        std::cout << "[Error] The sizes of positions and velocties don't match.\
            Aborting veloctiy verlet algorithm !"
                  << std::endl;
        return;
    }

    // initial energy/force calculation
    e_pot = force(positions, forces, num_particles);
    if (sampler != nullptr) {
        sampler(0, t, positions, velocities, e_pot, data);
    }
    // run the integration over time
    for (uint i = 1; i < num_t_steps; i++) {
        // assume forces are already for current positons
        positions += velocities * time_step + forces * time_step * time_step / 2;
        // add last force term
        velocities += forces * time_step / 2;
        e_pot = force(positions, forces, num_particles);
        // add current force term
        velocities += forces * time_step / 2;
        t += time_step;
        // call the sampler to save quantities in data
        if (sampler != nullptr) {
            sampler(i, t, positions, velocities, e_pot, data);
        }
    }
}

void velocity_verlet(Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double (*force)(const Eigen::ArrayX3d &,
                                     Eigen::ArrayX3d &,
                                     uint,
                                     double),
                     double time_step,
                     uint num_t_steps,
                     uint num_particles,
                     void (*sampler)(uint,
                                     double,
                                     Eigen::ArrayX3d &,
                                     Eigen::ArrayX3d &,
                                     double,
                                     Eigen::ArrayXXd &,
                                     double),
                     Eigen::ArrayXXd &data,
                     double mic_length) {
    Eigen::ArrayX3d forces;
    double e_pot;
    double t = 0;

    std::cout << "[Debug] Velocity verlet algorithm with "
              << num_t_steps << " steps."
              << std::endl;

    if (positions.size() != velocities.size()) {
        std::cout << "[Error] The sizes of positions and velocties don't match.\
            Aborting veloctiy verlet algorithm !"
                  << std::endl;
        return;
    }

    // initial energy/force calculation
    e_pot = force(positions, forces, num_particles, mic_length);
    if (sampler != nullptr) {
        sampler(0, t, positions, velocities, e_pot, data, mic_length);
    }

    // run the integration over time
    for (uint i = 1; i < num_t_steps; i++) {
        // assume forces are already for current positons
        positions += velocities * time_step + forces * time_step * time_step / 2;
        // add last force term
        velocities += forces * time_step / 2;
        e_pot = force(positions, forces, num_particles, mic_length);
        // add current force term
        velocities += forces * time_step / 2;
        t += time_step;
        // call the sampler to save quantities in data
        if (sampler != nullptr) {
            sampler(i, t, positions, velocities, e_pot, data, mic_length);
        }
    }
}

void sample_x2(uint index,
               double time,
               Eigen::ArrayX3d &positions,
               Eigen::ArrayX3d &velocities,
               double e_pot,
               Eigen::ArrayXXd &data) {
    data(index, 0) = time;
    // particle 1 - x components
    data(index, 1) = positions(0, 0);
    data(index, 2) = velocities(0, 0);
    // particle 2 - x components
    data(index, 3) = positions(1, 0);
    data(index, 4) = velocities(1, 0);
    // energies
    data(index, 5) = e_pot;
    // only need to add the mass later
    data(index, 6) = velocities.square().sum() / 2;
}

void sample_x2_wrapped(uint index,
                       double time,
                       Eigen::ArrayX3d &positions,
                       Eigen::ArrayX3d &velocities,
                       double e_pot,
                       Eigen::ArrayXXd &data,
                       double side_length) {
    Eigen::ArrayX3d pos_wrapped = positions;
    MD::coordinate_wrapping(pos_wrapped, side_length);
    data(index, 0) = time;
    // particle 1 - x components
    data(index, 1) = pos_wrapped(0, 0);
    data(index, 2) = velocities(0, 0);
    // particle 2 - x components
    data(index, 3) = pos_wrapped(1, 0);
    data(index, 4) = velocities(1, 0);
    // energies
    data(index, 5) = e_pot;
    // only need to add the mass later
    data(index, 6) = velocities.square().sum() / 2;
}

void sample_energies(uint index,
                     double time,
                     Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double e_pot,
                     Eigen::ArrayXXd &data) {
    data(index, 0) = time;
    data(index, 1) = e_pot;
    data(index, 2) = velocities.square().sum() / 2;
}

void sample_energies(uint index,
                     double time,
                     Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double e_pot,
                     Eigen::ArrayXXd &data,
                     double side_length) {
    data(index, 0) = time;
    data(index, 1) = e_pot;
    data(index, 2) = velocities.square().sum() / 2;
}

}  // namespace MD
