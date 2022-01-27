#ifndef __TESTS_H__
#define __TESTS_H__

#include <Eigen/Core>

namespace MD {

// Test the LJ energy and force calculation.
//
// Parameters:
// - data: array to hold the (x, energy, force) in the columns
// - num_x_samples: determines how many steps in x are computed, it will be the
//  number of rows of the array data
//
// Description:
//  Two argon atoms separated between 3 and 10 angstrom. Compute the forces
//  and the energy for 100 steps in between.
void lj_test(Eigen::Array<double, Eigen::Dynamic, 5> &data,
             uint num_x_samples = 100);

// Test the velocity verlet algorithm with 2 bouncing particles.
//
// Parameters:
// - data: array to hold the (t, x1, v1, x2, v2, epot, ekin) in the columns
// - time_step: step size in time
// - num_t_steps: determines how many steps in t are computed, it will be the
//  number of rows of the array data
//
// Description:
//  Test two argon atoms separated by 4 angstrom bouncing in their LJ potential.
void vv_test(Eigen::ArrayXXd &data, double time_step, uint num_t_steps);

}  // namespace MD

#endif  // __TESTS_H__