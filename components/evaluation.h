#ifndef __EVALUATION_H__
#define __EVALUATION_H__

#include <Eigen/Core>

namespace MD {

/* Compute the histogram for the radial distribution function

Parameters:
    * positions: 
        Nx3 array of the positions
    * r_hist: 
        vector for containing the histogram counts of the radial bins
    * num_bins: 
        number of bins in the radial distribution, length of r_hist
    * box_length: 
        side length of the box, used for minimum image convention

Description:
    * Loop over the pairs of particles. 
    * Compute the distance using the minimum image convention. 
    * Count the distances in bins for a radius from 0 to box_length.
*/
void radial_distribution_hist(Eigen::ArrayX3d& positions,
                              Eigen::VectorXi& r_hist,
                              uint num_bins,
                              double box_length);

}  // namespace MD

#endif  // __EVALUATION_H__