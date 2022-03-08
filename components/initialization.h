#ifndef __INITIALIZATION_H__
#define __INITIALIZATION_H__

#include <Eigen/Core>

namespace MD {

/* Initialize the positions in a grid.

Parameters:
    * positions: 
        array of size n^3x3 where each row is one particles position
    * num_particles_cubic: 
    number of particles per axis, thus the total number of particles is $n^3 =
    N$ which corresponds to the number of rows of the positions array
    * distance: 
        distance between the particles along a coordiante axis
    * box_length: 
        side length of the cubic box

Return:
    Number of particles

Description:
    * Compute the total number of particles.
    * Place the particles inside a qubic lattice of the box $[0, L]^3$.
*/
uint init_positions_3d(Eigen::ArrayX3d &positions,
                       uint num_paricles_cubic,
                       double distance);

/* Initialize the velocities randomly.

Parameters:
    * velocities: 
        array of size Nx3 where N is the number of particles
    * num_particles: 
        number of particles, is the number of rows in the velocities array
    * max_velocity: 
        maximal absolute value for any velocity component
    * seed: 
        seed for the random number generator

Description:
    * Set all components of the velocities array with uniformly distributed 
    double numbers up to the maximum velocity.
*/
void init_velocities_3d(Eigen::ArrayX3d &velocities,
                        uint num_particles,
                        uint seed,
                        double max_veloctiy = 1);

/* Remove a drift velocity.

Parameters:
    * velocities: 
        array of Nx3 where N is the number of particles

Description:
    * Compute the average net velocity.
    * Subtract it from every individual particles velocity.
*/
void velocity_drift_removal(Eigen::ArrayX3d &velocities);

/* Rescale the velocities to match the target temperature.

Parameters:
    * velocites: 
        array of size Nx3 containing the velocities (N is the number of 
        particles)
    * num_particles: 
        number of particles (rows of the veloctities array)
    * temp_target: 
        target temperature for the system

Description:
    * Compute the temperature of the system, aka kinetic energy.
    * Rescale the velocities to match the target temperature.
*/
void velocity_rescaling(Eigen::ArrayX3d &velocities,
                        uint num_particles,
                        double temp_target);

/* Wrap the positions to be in the initial box

Parameters:
    * positions: 
        Nx3 array of the particle positions, N number of particles
    * box_length: 
        side length of the initial box

Description:
    * Wrap the position coordnates such that they fall in the 3D box of the
    specified side length and a vertex at zero.
*/
void coordinate_wrapping(Eigen::ArrayX3d &positions,
                         double box_length);

}  // namespace MD

#endif  // __INITIALIZATION_H__