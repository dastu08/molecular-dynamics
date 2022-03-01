#ifndef __CLASSES_H__
#define __CLASSES_H__

#include <Eigen/Core>
#include <string>

namespace MD {

class Simulation {
   private:
    uint n;
    uint num_particles = 0;
    double time_step;
    uint num_t_steps;
    uint num_bins;
    uint num_rescales;
    uint num_rescale_steps;
    double box_length;
    double temperature;
    double separation = 0;

    Eigen::ArrayX3d positions;
    Eigen::ArrayX3d velocities;
    Eigen::ArrayX3d forces;
    Eigen::ArrayXXd data;

    uint initialzedFlag = 0;

   public:
    /* Set the simulation parameters

    Parameters:
        * n:
            number of particles per dimension
        * time_step:
            step size of the time integration scheme
        * num_t_steps:
            number of time steps for the simulation after rescaling
        * box_length:
            side length of the spacial box of the particles
        * num_rescales:
            number of velocity rescales in the beginning
        * num_rescales_steps:
            number of time steps to run after each velocity rescaling
        * num_bins:
            number of bins for the radial distribution sampling

    Description:
        * Copy the values to the internal variables.
    */
    Simulation(uint n,
               double time_step,
               uint num_t_steps,
               double box_length,
               uint num_rescales,
               uint num_rescales_steps,
               uint num_bins = 200);
    ~Simulation() {}

    /* Print the simulation parameters

    Description:
        * Print the values of the internal variables which were set by the
        constructor.
    */
    void printInfo();

    /* Initialize positions, velocities and data arrays

    Parameters:
        * separation:
            initial particles separation along each spacial dimension
        * temperature:
            target temperature for the velocity rescaling
        * seed:
            seed for the random number generator

    Description:
        * Set the data array to zero but with the correct size.
        * Initialize the positions and update the number of particles.
        * Initialize the velocities.
        * Remove the drift of the velocities.
        * Set initializedFlag to 1.
    */
    void init(double separation, double temperature, uint seed);

    /* Run the time integration

    Description:
        * Check for initializedFlag to be 1.
        * Run the rescales of the velocities.
        * Run the overall sampling phase.
        * Print information after finishing the run.
    */
    void run();

    /* Save sampled data to file

    Parameters:
        * filepath:
            path and file prefix for the file where the data will be written to

    Description:
        * Create the header string for the radial distribuation bins.
        * Save the data array to the file specified by the file path, append
        the value of n to it.
    */
    void export2file(std::string filepath);
};

}  // namespace MD

namespace MC {

class Simulation {
   private:
    uint n;
    uint num_particles = 0;
    uint num_samples;
    uint num_bins;
    uint seed;
    double box_length;
    double temperature;
    double separation = 0;

    Eigen::ArrayX3d positions;
    Eigen::ArrayX3d forces;
    Eigen::ArrayXXd data;

    uint initialzedFlag = 0;

   public:
    /* Set the simulation parameters

    Parameters:
        * n:
            number of particles per dimension
        * num_samples:
            number of samples to generate with Metropolis
        * box_length:
            side length of the spacial box of the particles
        * temperature:
            target temperature for acceptance criterion in Metropolis
        * num_bins:
            number of bins for the radial distribution sampling

    Description:
        * Copy the values to the internal variables.
    */
    Simulation(uint n,
               uint num_samples,
               double box_length,
               double temperature,
               uint num_bins = 200);
    ~Simulation() {}

    /* Print the simulation parameters

    Description:
        * Print the values of the internal variables which were set by the
        constructor.
    */
    void printInfo();

    /* Initialize positions, velocities and data arrays

    Parameters:
        * separation:
            initial particles separation along each spacial dimension
        * seed:
            TODO

    Description:
        * Set the data array to zero but with the correct size.
        * Initialize the positions and update the number of particles.
        * Set initializedFlag to 1.
    */
    void init(double separation, uint seed);

    /* Run the time integration

    Description:
        * Check for initializedFlag to be 1.
        *
        * Print information after finishing the run.
    */
    void run();

    /* Save sampled data to file

    Parameters:
        * filepath:
            path and file prefix for the file where the data will be written to

    Description:
        * Create the header string for the radial distribuation bins.
        * Save the data array to the file specified by the file path, append
        the value of n to it.
    */
    void export2file(std::string filepath);
};

}  // namespace MC

#endif  // __CLASSES_H__