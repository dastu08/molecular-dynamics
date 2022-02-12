#ifndef __SAMPLERS_H__
#define __SAMPLERS_H__

#include <Eigen/Core>

namespace MD {

// Sample the x-components of the first two particles.
//
// Parameters:
// - index: row index of data where the quantities will be written to
// - time: current time during the sampling
// - positions: array of particle positions of size Nx3, N is the nubmer of
//  particles
// - velocities: array of particle velocities, same size as positions
// - e_pot: potential energy of the whole system
// - data: array with at least 6 column where the data will be written to
//
// Description:
//  Save in the index row of data the quantities:
//  time, x1, v1, x2, v2, e_pot, e_kin.
//  e_pot and e_kin are computed using the whole system.
void sample_x2(uint index,
               double time,
               Eigen::ArrayX3d &positions,
               Eigen::ArrayX3d &velocities,
               double e_pot,
               Eigen::ArrayXXd &data);

void sample_x2_wrapped(uint index,
                       double time,
                       Eigen::ArrayX3d &positions,
                       Eigen::ArrayX3d &velocities,
                       double e_pot,
                       Eigen::ArrayXXd &data,
                       double side_length);

// Sample the potential and kinetic energies.
//
// Parameters:
// - index: row index of data where the quantities will be written to
// - time: current time during the sampling
// - positions: array of particle positions of size Nx3, N is the nubmer of
//  particles
// - velocities: array of particle velocities, same size as positions
// - e_pot: potential energy of the whole system
// - data: array with at least 3 column where the data will be written to
//
// Description:
//  Save in the index row of data the quantities: time, e_pot, e_kin.
//  e_pot and e_kin are computed using the whole system.
void sample_energies(uint index,
                     double time,
                     Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double e_pot,
                     Eigen::ArrayXXd &data);
void sample_energies(uint index,
                     double time,
                     Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double e_pot,
                     Eigen::ArrayXXd &data,
                     double side_length);

}  // namespace MD

#endif  // __SAMPLERS_H__