
#pragma once
#include <vector>

/**
 * Time-Weighted Dynamic Time Warping (TWDTW) distance between two
 * multivariate time-series.
 *
 * @param q         Query series as `n x b` (time x bands).
 * @param p         Pattern series as `m x b` (time x bands).
 * @param q_doy     Day-of-year of each query point (length `n`).
 * @param p_doy     Day-of-year of each pattern point (length `m`).
 * @param dist_power `p` value of the Minkowski local distance over bands.
 * @param weight_type 0 = none, 1 = logistic, 2 = linear.
 * @param alpha,beta Logistic time-weight parameters.
 * @param a,b       Linear time-weight parameters.
 *
 * @return The path-length-normalized TWDTW distance.
 */
double twdtw_distance_op(const std::vector<std::vector<double>> &q,
                         const std::vector<std::vector<double>> &p,
                         const std::vector<double> &q_doy,
                         const std::vector<double> &p_doy,
                         double dist_power,
                         int weight_type,
                         double alpha, double beta,
                         double a, double b);
