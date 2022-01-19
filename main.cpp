#include <Eigen/Core>
#include <iostream>

#include "components/potentials.h"

int main() {
    Eigen::ArrayX3d positions, forces;
    const uint num_particles = 2;
    double energy;

    // report 1. dynamics

    // compute the LJ energy
    energy = MD::lennard_jones(positions, forces, num_particles);

    return 0;
}