#ifndef __INITIALIZATION_H__
#define __INITIALIZATION_H__

#include <Eigen/Core>

namespace MD {

// Initialize the positions in a grid
//
// Parameters:
// - positions: array of size n^3x3 where each row is one particles position
// - num_particles_cubic: number of particles per axis, thus the total number of
//   particles is $n^3 = N$ which corresponds to the number of rows of the
//   positions array
// - distance: distance between the particles along a coordiante axis
// - box_length: side length of the cubic box
//
// Return:
//  Number of particles
//
// Description:
//  Place the particles inside a qubic lattice of the box $[0, L]^3$.
uint init_potisionts_3d(Eigen::ArrayX3d &positions,
                        uint num_paricles_cubic,
                        double distance);

// Initialize the velocities randomly
//
// Parameters:
// - velocities: array of size Nx3 where N is the number of particles
// - num_particles: number of particles, is the number of rows in the velocities
//   array
// - max_velocity: maximal absolute value for any velocity component
// - seed: seed for the random number generator
//
// Description:
//  Set all components of the velocities array with uniformly distributed double
//  numbers up to the maximum velocity.
void init_velocities_3d(Eigen::ArrayX3d &velocities,
                        uint num_particles,
                        double max_veloctiy,
                        uint seed);

// Remove a drift velocity
//
// Parameters:
// - velocities: array of Nx3 where N is the number of particles
//
// Description:
//  Compute the average net velocity and subtract it from every individual
//  particles velocity.
void velocity_drift_removal(Eigen::ArrayX3d &velocities);

}  // namespace MD

#endif  // __INITIALIZATION_H__