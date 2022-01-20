#include <Eigen/Core>
#include <iostream>

#include "components/helpers.h"
#include "components/potentials.h"

int main() {
    // report 1. dynamics

    // initialize the positions
    Eigen::ArrayX3d positions{{0, 0, 0}, {3.5, 0, 0}};
    Eigen::ArrayX3d forces;
    const uint num_particles = 2;
    double energy;

    std::cout << "initial positions\n"
              << positions << std::endl;

    MD::array2file(positions, "test.txt", "x,y,z");

    // compute the LJ energy
    // energy = MD::lennard_jones(positions, forces, num_particles);

    return 0;
}