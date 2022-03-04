#include "solvers.h"

#include <math.h>
#include <stdlib.h>

#include <iostream>

#include "evaluation.h"
#include "helpers.h"
#include "potentials.h"

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

uint move(Eigen::ArrayX3d &positions,
          uint num_particles,
          double step_size) {
    //   select a random particle
    uint particle = rand() % num_particles;
    Eigen::Vector3d delta;

    delta(0) = step_size * MC::random_double_2();
    delta(1) = step_size * MC::random_double_2();
    delta(2) = step_size * MC::random_double_2();

    // move the particle in a random direction
    positions.row(particle) += delta.transpose().array();

    return particle;
}

void metropolis(Eigen::ArrayX3d &positions,
                double (*potential)(const Eigen::ArrayX3d &,
                                    uint,
                                    double,
                                    Eigen::VectorXi &,
                                    double),
                uint num_samples,
                uint num_particles,
                double step_size,
                double beta,
                uint seed,
                void (*sampler)(uint,
                                Eigen::ArrayX3d &,
                                double,
                                Eigen::VectorXi &,
                                Eigen::ArrayXXd &,
                                double,
                                uint),
                Eigen::ArrayXXd &data,
                double box_length,
                uint num_bins,
                uint index_offset) {
    // init the random number generator
    srand(seed);
    double reference, relative_probability, accpetance_rate;
    double energy_new, energy_old, energy_tot;
    uint accpetance_count = 0;
    Eigen::ArrayX3d positions_old = positions;
    Eigen::VectorXi r_hist;
    uint particle;

    std::cout << "[Debug] Metropolis algorithm with "
              << num_samples << " sampling steps."
              << std::endl;

    // initial energy
    energy_tot = potential(positions,
                           num_particles,
                           box_length,
                           r_hist,
                           num_bins);

    if (sampler != nullptr) {
        sampler(index_offset,
                positions,
                energy_tot,
                r_hist,
                data,
                box_length,
                num_bins);
    }

    for (uint i = 0; i < num_samples; ++i) {
        // trial move and compute the new energy
        particle = move(positions, num_particles, step_size);
        energy_old = MC::lennard_jones_single(positions_old, num_particles, box_length, particle);
        energy_new = MC::lennard_jones_single(positions, num_particles, box_length, particle);

        // compute the total initial energy
        // energy_new = potential(positions,
        //                        num_particles,
        //                        box_length,
        //                        r_hist_new,
        //                        num_bins);
        // std::cout << r_hist_new.transpose() << std::endl;

        // accept the new configuration?
        reference = MC::random_double_1();
        relative_probability = exp(-beta * (energy_new - energy_old));
        // acceptance rate: min(1, relative_probability)
        if (relative_probability > reference) {
            // count the number of new accepted configuations
            ++accpetance_count;
            // update the energy/r_hist for the next move
            energy_tot += energy_new - energy_old;
            // energy_old = energy_new;
            // r_hist_old = r_hist_new;
            MD::radial_distribution_hist(positions, r_hist, num_bins, box_length);

            // keep a copy of the configuration
            positions_old = positions;
            // std::cout << "accepted energy: " << energy_new << std::endl;
        } else {
            // keep the old energy/rdf
            // energy_new = energy_old;
            // r_hist_new = r_hist_old;
            // reset the positions
            positions = positions_old;
            // std::cout << "rejected energy: " << energy_new << std::endl;
        }

        if (sampler != nullptr) {
            sampler(i + index_offset,
                    positions,
                    energy_tot,
                    r_hist,
                    data,
                    box_length,
                    num_bins);
        }
    }

    accpetance_rate = (double)accpetance_count / (double)num_samples;
    std::cout << "[Debug] Metropolis acceptance rate: "
              << accpetance_rate << std::endl;
}

}  // namespace MC