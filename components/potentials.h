#ifndef __POTENTIALS_H__
#define __POTENTIALS_H__

#include <Eigen/Core>

namespace MD {

// Compute the LJ energy and forces for N particles in 3 dimensions.
//
// Parameters:
//  - positions : array of particle positions of size Nx3, N is the number of
//  particles
// - forces: array of the total forces on each particle, same size as positions
// - num_particles: number of particles, equivalent to the number of rows of the
//  arrays positions and forces.
//
// Description:
//  For each particle pair compute the forces from the Lennard-Jones potential
//  and sum them all up. Also compute the total potential energy of the system.
double lennard_jones(const Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &forces,
                     uint num_particles);

}  // namespace MD

#endif  // __POTENTIALS_H__