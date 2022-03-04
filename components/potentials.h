#ifndef __POTENTIALS_H__
#define __POTENTIALS_H__

#include <Eigen/Core>

namespace MD {

/* Compute the LJ energy and forces for N particles in 3 dimensions.

Parameters:
    * positions:
        array of particle positions of size Nx3, N is the number of particles
    * forces:
        array of the total forces on each particle, same size as positions
    * num_particles:
        number of particles, equivalent to the number of rows of the arrays
        positions and forces.

Description:
    * For each particle pair compute the forces from the Lennard-Jones potential
        and sum them all up.
    * Compute the total potential energy of the system.
        V_{LJ}(r) = 4 epsilon ( (sigma/r)^{12} - (sigma/r)^6 )
    * positions in units of sigma.
    * energy in units of epsilon.
    * forces in units of epsilon/sigma.
 */
double lennard_jones(const Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &forces,
                     uint num_particles);

/*Compute the LJ energy and forces for N particles in 3 dimensions.

Overload Parameters:
    * side_length:
        side length of the box for the minimum image convention

Description:
    * Use the minimum image convention for computing the relative distance.
*/
double lennard_jones(const Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &forces,
                     uint num_particles,
                     double side_length);

}  // namespace MD

namespace MC {

/* Compute the LJ energy for N particles in 3 dimensions.

Parameters:
    * positions:
        array of particle positions of size Nx3, N is the number of particles
    * num_particles:
        number of particles, equivalent to the number of rows of the positions
        array
    * side_length:
        side length of the box for the minimum image convention

Description:
    * For each particle pair compute relative distance.
    * Apply the minimum image convention.
    * Compute the LJ energy
        V_{LJ}(r) = 4 epsilon ( (sigma/r)^{12} - (sigma/r)^6 )
    * Add everything up to the total potential energy of the system.
    * positions in units of sigma.
    * energy in units of epsilon.
*/
double lennard_jones(const Eigen::ArrayX3d &positions,
                     uint num_particles,
                     double side_length);

/* Compute the LJ energy of one particle interacting with the others

Parameters:
    * positions:
        array of particle positions of size Nx3, N is the number of particles
    * num_particles:
        number of particles, equivalent to the number of rows of the positions
        array
    * side_length:
        side length of the box for the minimum image convention
    * particle:
        index of the particle which determines the interaction with the others

Description:
    * Compute relative distance to each other particle.
    * Apply the minimum image convention.
    * Compute the LJ energy
        V_{LJ}(r) = 4 epsilon ( (sigma/r)^{12} - (sigma/r)^6 )
    * Add everything up to the total potential energy of the single particle.
    * positions in units of sigma.
    * energy in units of epsilon.
*/
double lennard_jones_single(const Eigen::ArrayX3d &positions,
                            uint num_particles,
                            double side_length,
                            uint particle);

}  // namespace MC

#endif  // __POTENTIALS_H__