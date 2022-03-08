#include "samplers.h"

#include <Eigen/Core>

#include "evaluation.h"
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

void sample_energies_rdf(uint index,
                         double time,
                         Eigen::ArrayX3d &positions,
                         Eigen::ArrayX3d &velocities,
                         double e_pot,
                         Eigen::ArrayXXd &data,
                         double side_length,
                         uint num_bins) {
    Eigen::VectorXi r_hist;
    MD::radial_distribution_hist(positions, r_hist, num_bins, side_length);

    data(index, 0) = time;
    data(index, 1) = e_pot;
    data(index, 2) = velocities.square().sum() / 2;
    data(index, Eigen::seqN(3, num_bins)) = r_hist.cast<double>().transpose();
}

}  // namespace MD

namespace MC {

void sample_energies_rdf(uint index,
                         Eigen::ArrayX3d &positions,
                         double e_pot,
                         Eigen::VectorXi &r_hist,
                         Eigen::ArrayXXd &data,
                         double side_length,
                         uint num_bins) {
    data(index, 0) = index;
    data(index, 1) = e_pot;
    data(index, Eigen::seqN(2, num_bins)) = r_hist.cast<double>().transpose();
}

}  // namespace MC
