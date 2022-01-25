#ifndef __SOLVERS_H__
#define __SOLVERS_H__

#include <Eigen/Core>

namespace MD {

void velocity_verlet(Eigen::ArrayX3d &positions,
                     Eigen::ArrayX3d &velocities,
                     double (*force)(const Eigen::ArrayX3d &,
                                     Eigen::ArrayX3d &,
                                     uint),
                     double time_step,
                     uint num_t_steps,
                     uint num_particles,
                     void (*sampler)(uint,
                                     double,
                                     Eigen::ArrayX3d &,
                                     Eigen::ArrayX3d &,
                                     double,
                                     Eigen::ArrayXXd &),
                     Eigen::ArrayXXd &data);

void sample_x(uint index,
              double time,
              Eigen::ArrayX3d &positions,
              Eigen::ArrayX3d &velocities,
              double e_pot,
              Eigen::ArrayXXd &data);

void vv_test();

}  // namespace MD

#endif  // __SOLVERS_H__