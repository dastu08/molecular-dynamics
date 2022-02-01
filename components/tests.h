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
//  Test two argon atoms separated by 5 angstrom bouncing in their LJ potential.
void vv_test(Eigen::ArrayXXd &data, double time_step, uint num_t_steps);

// Test the accuracy of the velocity verlet algorithm on the time step.
//
// Parameters:
// - data: array to hold the the columns (index, time_step)
// - time_step_0: start step size in time
// - time_step_1: stop step size in time
// - num_t_step_samples: determines how many steps in time_step are computed, 
//  it will be the number of rows of the array data
//
// Description:
//  Test two argon atoms separated by 5 angstrom bouncing in their LJ potential.
//  Vary the time step and record the list of potential and kintetic energy over
//  time (2000 steps) in a separate file.
void time_step_test(Eigen::ArrayX2d &data,
                    double time_step_0,
                    double time_step_1,
                    uint num_t_step_samples);

}  // namespace MD

#endif  // __TESTS_H__