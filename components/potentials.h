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
//  $V_{LJ}(r) = 4 eps ( (sigma/r)^{12} - (sigma/r)^6 )$
//  The implementation has positions in units of sigma and the energy/forces
//  in units of eps.
double lennard_jones(const Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &forces,
                     uint num_particles);

// Test the LJ energy and force calculation.
//
// Parameters:
//  - data: array to hold the (x, energy, force) in the columns
//  - num_x_samples: determines how many steps in x are computed, it will be the
//  number of rows of the array data
// 
// Description:
//  Two particles separated between 3 and 10 angstrom. Compute the forces
//  and the energy for 100 steps in between.
void lj_test(Eigen::Array<double, Eigen::Dynamic, 5> &data,
             uint num_x_samples = 100);

}  // namespace MD

#endif  // __POTENTIALS_H__