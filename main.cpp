#include <Eigen/Core>
#include <iostream>

#include "components/helpers.h"
#include "components/potentials.h"

int main() {
    // report 1. dynamics

    // initialize the positions
    // Eigen::ArrayX3d positions{{0, 0, 0}, {3.5, 0, 0}};
    // Eigen::ArrayX3d forces;
    // const uint num_particles = 2;
    // double energy;

    // std::cout << "initial positions\n"
    //           << positions << std::endl;

    // MD::array2file(positions, "test.txt", "x,y,z");

    // // compute the LJ energy
    // energy = MD::lennard_jones(positions, forces, num_particles);

    // std::cout << "Total energy: " << energy << std::endl;

    uint num_samples = 200;
    Eigen::Array<double, Eigen::Dynamic, 5> data(num_samples, 5);
    MD::lj_test(data, num_samples);
    MD::array2file(data, "../reports/lj_test.txt", "x,energy,fx,fy,fz");

    return 0;
}