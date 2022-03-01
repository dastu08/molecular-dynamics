#include "solvers.h"

#include <stdlib.h>

#include <iostream>

#include "helpers.h"

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
                     double box_length) {
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
    e_pot = force(positions, forces, num_particles, box_length);
    if (sampler != nullptr) {
        sampler(0, t, positions, velocities, e_pot, data, box_length);
    }

    // run the integration over time
    for (uint i = 1; i < num_t_steps; i++) {
        // assume forces are already for current positons
        positions += velocities * time_step + forces * time_step * time_step / 2;
        // add last force term
        velocities += forces * time_step / 2;
        e_pot = force(positions, forces, num_particles, box_length);
        // add current force term
        velocities += forces * time_step / 2;
        t += time_step;
        // call the sampler to save quantities in data
        if (sampler != nullptr) {
            sampler(i, t, positions, velocities, e_pot, data, box_length);
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
                                     double,
                                     uint),
                     Eigen::ArrayXXd &data,
                     double box_length,
                     uint num_bins,
                     uint index_offset) {
    Eigen::ArrayX3d forces;
    double e_pot;
    double t = index_offset * time_step;

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
    e_pot = force(positions, forces, num_particles, box_length);
    if (sampler != nullptr) {
        sampler(index_offset, t, positions, velocities, e_pot, data, box_length, num_bins);
    }

    // run the integration over time
    for (uint i = 1; i < num_t_steps; i++) {
        // assume forces are already for current positons
        positions += velocities * time_step + forces * time_step * time_step / 2;
        // add last force term
        velocities += forces * time_step / 2;
        e_pot = force(positions, forces, num_particles, box_length);
        // add current force term
        velocities += forces * time_step / 2;
        t += time_step;
        // call the sampler to save quantities in data
        if (sampler != nullptr) {
            sampler(i + index_offset, t, positions, velocities, e_pot, data, box_length,
                    num_bins);
        }
    }
}

}  // namespace MD

namespace MC {

Eigen::Vector3d move(Eigen::ArrayX3d &positions_old,
                     Eigen::ArrayX3d &postions_new,
                     uint num_particles,
                     double step_size) {
    //   select a random particle
    uint particle = rand() % num_particles;

    Eigen::Vector3d delta;

    delta(0) = step_size * MC::random_double();
    delta(1) = step_size * MC::random_double();
    delta(2) = step_size * MC::random_double();

    return delta;
}

void metropolis(Eigen::ArrayX3d &positions,
                uint num_samples,
                uint num_particles,
                uint seed,
                void (*sampler)(),
                Eigen::ArrayXXd &data,
                double box_length,
                uint num_bins,
                uint index_offset) {
    // init the random number generator
    srand(seed);
    uint particle;
    Eigen::Vector3d delta;

    for (uint i = 0; i < num_samples; ++i) {
        // particle = rand() % num_particles;
        delta = move(positions, positions, num_particles, 0.4);
        data(i, 0) = i;
        data(i, 1) = delta(0);
        data(i, 2) = delta(1);
        data(i, 3) = delta(2);

        // std::cout << "particle : " << particle << std::endl;
    }
}

}  // namespace MC