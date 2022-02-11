#ifndef __SOLVERS_H__
#define __SOLVERS_H__

#include <Eigen/Core>

namespace MD {

// Veloctiy Verlet integration of equations of motion.
//
// Paramters:
// - positions: array of particle positions of size Nx3, N is the nubmer of
//  particles
// - velocities: array of particle velocities, same size as positions
// - forces: array of the total forces on each particle, same size as positions
// - force: function with arguments (positions, forces, num_particles) which
//  computes the forces and returns the potential energy
// - time_step: step size in time
// - num_t_steps: number of steps in time, i.e. how often the loop runs
// - num_particles: number of particles, equivalent to the number of rows of the
//  arrays positions and forces
// - sampler: function with parameters
//  (index, time, positions, velocites, e_pot, data) where the relevant
//  quantities are written into data
//  - data: array where the data from the sampler is written to
//
// Description:
//  Starting at time t = 0 integrate equations of motions $mx''(t) = F(t)$ with the
//  veloctiy verlet algorithm:
//  $x(t + h) = x(t) + h*v(t) + h^2 * F(t) / (2m)$
//  $v(t + h) = v(t) + h*(F(t) + F(t + h)) / (2m)$
//  All quantities are used in reduces units. If the positions are in units of
//  sigma and the energy in units of epsilon. Then the force needs to be in
//  units of epsilon/sigma, the velocity will be in units of $\sqrt{epsilon/m}$
//  and the time is in units $sigma \sqrt{m/epsilon}$.
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
                     Eigen::ArrayXXd &data);

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
                     double mic_length);

// Sample the x-components of the first two particles.
//
// Parameters:
// - index: row index of data where the quantities will be written to
// - time: current time during the sampling
// - positions: array of particle positions of size Nx3, N is the nubmer of
//  particles
// - velocities: array of particle velocities, same size as positions
// - e_pot: potential energy of the whole system
// - data: array with at least 6 column where the data will be written to
//
// Description:
//  Save in the index row of data the quantities:
//  time, x1, v1, x2, v2, e_pot, e_kin.
//  e_pot and e_kin are computed using the whole system.
void sample_x2(uint index,
               double time,
               Eigen::ArrayX3d &positions,
               Eigen::ArrayX3d &velocities,
               double e_pot,
               Eigen::ArrayXXd &data);

void sample_x2_wrapped(uint index,
                       double time,
                       Eigen::ArrayX3d &positions,
                       Eigen::ArrayX3d &velocities,
                       double e_pot,
                       Eigen::ArrayXXd &data,
                       double side_length);

// Sample the potential and kinetic energies.
//
// Parameters:
// - index: row index of data where the quantities will be written to
// - time: current time during the sampling
// - positions: array of particle positions of size Nx3, N is the nubmer of
//  particles
// - velocities: array of particle velocities, same size as positions
// - e_pot: potential energy of the whole system
// - data: array with at least 3 column where the data will be written to
//
// Description:
//  Save in the index row of data the quantities: time, e_pot, e_kin.
//  e_pot and e_kin are computed using the whole system.
void sample_energies(uint index,
                     double time,
                     Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double e_pot,
                     Eigen::ArrayXXd &data);
void sample_energies(uint index,
                     double time,
                     Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double e_pot,
                     Eigen::ArrayXXd &data,
                     double side_length);
}  // namespace MD

#endif  // __SOLVERS_H__