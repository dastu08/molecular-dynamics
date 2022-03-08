#ifndef __TESTS_H__
#define __TESTS_H__

#include <Eigen/Core>
#include <string>

namespace MD {

/* Test the LJ energy and force calculation.

Parameters:
    * data:
        array to hold the (x, energy, force) in the columns
    * num_x_samples:
        determines how many steps in x are computed, it will be the
        number of rows of the array data

Description:
    * Two argon atoms separated between 3 and 10 angstrom.
    * Compute the forces and the energy for 100 steps in between.
*/
void lj_test(Eigen::Array<double, Eigen::Dynamic, 5> &data,
             uint num_x_samples = 100);

/* Test the velocity verlet algorithm with 2 bouncing particles.

Parameters:
    * data:
        array to hold the (t, x1, v1, x2, v2, epot, ekin) in the columns
    * time_step:
        step size in time
    * num_t_steps:
        determines how many steps in t are computed, it will be the number of
        rows of the array data

Description:
    * Test two argon atoms separated by 5 angstrom bouncing in their LJ
    potential.
*/
void vv_test(Eigen::ArrayXXd &data, double time_step, uint num_t_steps);

/* Test the accuracy of the velocity verlet algorithm on the time step.

Parameters:
    * data:
        array to hold the the columns (index, time_step)
    * time_step_0:
        start step size in time
    * time_step_1:
        stop step size in time
    * num_t_step_samples:
        determines how many steps in time_step are computed, it will be the
        number of rows of the array data

Description:
    * Test two argon atoms separated by 5 angstrom bouncing in their LJ
    potential.
    * Vary the time step and record the list of potential and kintetic energy
    over time (2000 steps) in a separate file.
*/
void time_step_test(Eigen::ArrayX2d &data,
                    double time_step_0,
                    double time_step_1,
                    uint num_t_step_samples);

/* Test the LJ energy and force calculation with minimal image convention.

Parameters:
    * data:
        array to hold the (x, energy, force) in the columns
    * num_x_samples:
        determines how many steps in x are computed, it will be the number of
        rows of the array data

Description:
    * Two argon atoms separated between 3 and 13 angstrom.
    * Compute the forces and the energy for 100 steps in between using the
    minimal image convention with side length 13 angstrom.
*/
void lj_test2(Eigen::Array<double, Eigen::Dynamic, 5> &data,
              uint num_x_samples);

/* Test 2 bouncing particles with the minimal image convention

Parameters:
    * data: 
        array to hold the (t, x1, v1, x2, v2, epot, ekin) in the columns
    * time_step: 
        step size in time
    * num_t_steps: 
        determines how many steps in t are computed, it will be the number of 
        rows of the array data
    * side_length: 
        side length of the minium image convention box

Description:
    * Test two argon atoms separated by 10 angstrom bouncing in their LJ
    potential computed with the minimal image convention.
*/
void mic_test(Eigen::ArrayXXd &data,
              double time_step,
              uint num_t_steps,
              double side_length);

/* Equilibrate the system afer velocity rescaling

Parameters:
    * positions: 
        Nx3 array of the positons
    * velocties: 
        Nx3 array of the velocities
    * time_step: 
        step size in time for velocity verlet
    * num_t_steps: 
        number of time steps to run velocity verlet
    * num_particles: 
        number of particles
    * temperature: 
        target temperature of the system
    * box_length: 
        side length of the box for periodic boundary conditions
    * filename: 
        name of the text file containing the energies sampled during the
        equilibration phase

Description:
    * Rescale the velocities using the target temperature. 
    * Run velocity verlet.
    * Save the sampled energies (t, epot, ekin) in a file.
*/
void equilibration_phase(Eigen::ArrayX3d &positions,
                         Eigen::ArrayX3d &velocities,
                         double time_step,
                         double num_t_steps,
                         uint num_particles,
                         double temperature,
                         double box_length,
                         std::string filename);

}  // namespace MD

#endif  // __TESTS_H__