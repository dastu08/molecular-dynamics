#ifndef __INITIALIZATION_H__
#define __INITIALIZATION_H__

#include <Eigen/Core>

namespace MD {

// Initialize the positions and velocities
//
// Parameters:
//  - positions : array of size n^3x3 where each row is one particles position
//  - num_particles_cubic: number of particles per axis, thus the total number
//  of particles is $n^3 = N$ which corresponds to the number of rows of the
//  positions array
//  - distance: distance between the particles along a coordiante axis
//  - box_length: side length of the cubic box
//
// Description:
//  Place the particles inside a qubic lattice of the box $[0, L]^3$.
void init3d(Eigen::ArrayX3d &positions,
            uint num_paricles_cubic,
            double distance);

}  // namespace MD

#endif  // __INITIALIZATION_H__