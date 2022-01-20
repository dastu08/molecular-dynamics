#ifndef __INITIALIZATION_H__
#define __INITIALIZATION_H__

#include <Eigen/Core>

namespace MD {

// Initialize the positions and velocities
//
// Parameters:
//  - positions : array of size Nx3 where each row is one particles position
//  - num_particles: number of particles which corresponds to the number of rows
//      of the positions array
// 
// Description:
//  Place the particles inside a qubic lattice of the unit box $[0, 1]^3$.
// void init3d(Eigen::ArrayX3d &positions, uint num_paricles);

}  // namespace MD

#endif  // __INITIALIZATION_H__