#ifndef __SOLVERS_H__
#define __SOLVERS_H__

#include <Eigen/Core>

namespace MD {

/* Veloctiy Verlet integration of equations of motion.

Parameters:
    * positions:
        array of particle positions of size Nx3, N is the nubmer of particles
    * velocities:
        array of particle velocities, same size as positions
    * forces:
        array of the total forces on each particle, same size as positions
    * force:
        function with arguments (positions, forces, num_particles) which
        computes the forces and returns the potential energy
    * time_step:
        step size in time
    * num_t_steps:
        number of steps in time, i.e. how often the loop runs
    * num_particles:
        number of particles, equivalent to the number of rows of the arrays
        positions and forces
    * sampler:
        function with parameters (index, time, positions, velocites, e_pot,
        data) where the relevant quantities are written into data
    * data:
        array where the data from the sampler is written to

Description:
    * Starting at time t = 0 integrate equations of motions mx''(t) = F(t) with
    the veloctiy verlet algorithm:
    $$
        x(t + h) = x(t) + h*v(t) + h^2 * F(t) / (2m)
        v(t + h) = v(t) + h*(F(t) + F(t + h)) / (2m)
    $$
    * All quantities are used in reduces units.
        * [force] = [energy] / [length]
        * [velocity] = sqrt([energy] / [mass])
        * [time] = [length] * sqrt([mass] / [energy])
*/
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

/* Veloctiy Verlet extended with minimum image convention

Overload Parameters:
    * force:
        function with arguments (positions, forces, num_particles, box_length)
        which computes the forces and returns the potential energy
    * sampler:
        function with parameters (index, time, positions, velocites, e_pot,
        data, box_length) where the relevant quantities are written into data
    * box_length:
        side lenght of the box in the minimum image convention

Description:
    * Pass the box_length to the force calulation and the sampler function
    so they can use the minimum image convention.
*/
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
                     double box_length);

/* Veloctiy Verlet extended with radial distribution sampling

Overload Parameters:
    * sampler:
        function with parameters (index, time, positions, velocites, e_pot,
        data, box_length, num_bins) where the relevant quantities are written
        into data
    * numb_bins:
        number of bins used for the radial distribution sampling
    * index_offset:
        offset in the number of time steps so the sampling gets corrected time
        stamps

Description:
    * Set the initial time according to the index_offset.
    * Pass the num_bins to the sampler function so it can sample the radial
    distribution.
*/
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
                     uint index_offset = 0);
}  // namespace MD

namespace MC {

/* Perform a trial move on the positions

Parameters:
    * positions:
        Nx3 array of the positions, has rows as N=num_particles
    * num_particles:
        number of particles of the system, number of rows of positions
    * step_size:
        steps size for the random movement

Description:
    * Pick a random particle.
    * Move the particle in a random direction inside the box
        [-step_size, step_size]^3 around the previous position.
    * Override the particles position in positions array.
    * Return the moved particle index.
*/
uint move(Eigen::ArrayX3d &positions,
          uint num_particles,
          double step_size);

/* Metropolis Monte Carlo algorithm

Parameters:
    * positions:
    * potential:
    * num_samples:
    * num_particles:
    * step_size:
    * beta:
    * seed:
    * sampler:
    * data:
    * box_length:
    * num_bins:
    * index_offset
 */
void metropolis(Eigen::ArrayX3d &positions,
                double (*potential)(const Eigen::ArrayX3d &,
                                    uint,
                                    double),
                double (*potential_single)(const Eigen::ArrayX3d &,
                                           uint,
                                           double,
                                           uint),
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
                void (*radial_hist)(Eigen::ArrayX3d &,
                                    Eigen::VectorXi &,
                                    uint,
                                    double),
                uint index_offset = 0);

}  // namespace MC

#endif  // __SOLVERS_H__