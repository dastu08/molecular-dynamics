#include "samplers.h"

#include "initialization.h"

namespace MD {

void sample_x2(uint index,
               double time,
               Eigen::ArrayX3d &positions,
               Eigen::ArrayX3d &velocities,
               double e_pot,
               Eigen::ArrayXXd &data) {
    data(index, 0) = time;
    // particle 1 - x components
    data(index, 1) = positions(0, 0);
    data(index, 2) = velocities(0, 0);
    // particle 2 - x components
    data(index, 3) = positions(1, 0);
    data(index, 4) = velocities(1, 0);
    // energies
    data(index, 5) = e_pot;
    // only need to add the mass later
    data(index, 6) = velocities.square().sum() / 2;
}

void sample_x2_wrapped(uint index,
                       double time,
                       Eigen::ArrayX3d &positions,
                       Eigen::ArrayX3d &velocities,
                       double e_pot,
                       Eigen::ArrayXXd &data,
                       double side_length) {
    Eigen::ArrayX3d pos_wrapped = positions;
    MD::coordinate_wrapping(pos_wrapped, side_length);
    data(index, 0) = time;
    // particle 1 - x components
    data(index, 1) = pos_wrapped(0, 0);
    data(index, 2) = velocities(0, 0);
    // particle 2 - x components
    data(index, 3) = pos_wrapped(1, 0);
    data(index, 4) = velocities(1, 0);
    // energies
    data(index, 5) = e_pot;
    // only need to add the mass later
    data(index, 6) = velocities.square().sum() / 2;
}

void sample_energies(uint index,
                     double time,
                     Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double e_pot,
                     Eigen::ArrayXXd &data) {
    data(index, 0) = time;
    data(index, 1) = e_pot;
    data(index, 2) = velocities.square().sum() / 2;
}

void sample_energies(uint index,
                     double time,
                     Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double e_pot,
                     Eigen::ArrayXXd &data,
                     double side_length) {
    data(index, 0) = time;
    data(index, 1) = e_pot;
    data(index, 2) = velocities.square().sum() / 2;
}

}  // namespace MD
