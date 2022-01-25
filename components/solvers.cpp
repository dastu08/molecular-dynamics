#include "solvers.h"

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
    // Eigen::VectorXd energies = Eigen::VectorXd::Zero(num_t_steps);
    double e_pot;
    double t = 0;
    // initial energy/force calculation
    e_pot = force(positions, forces, num_particles);
    sampler(0, t, positions, velocities, e_pot, data);

    for (int i = 1; i < num_t_steps; i++) {
        // assume forces are already for current positons
        positions += velocities * time_step + forces * time_step * time_step / 2;
        // add last force term
        velocities += forces * time_step / 2;
        e_pot = force(positions, forces, num_particles);
        // add current force term
        velocities += forces * time_step / 2;
        t += time_step;
        sampler(i, t, positions, velocities, e_pot, data);
    }
}

void sample_x(uint index,
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
    data(index, 6) = velocities.square().sum();
}

}  // namespace MD
