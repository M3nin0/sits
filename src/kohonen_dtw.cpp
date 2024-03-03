#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "./sits_types.h"

using namespace Rcpp;


/**
 * Minimum of 2 values.
 *
 * @description
 * Auxiliary function to calculate the minimum value of `x` and `y`.
 */
double minval(double x, double y)
{
    // z > nan for z != nan is required by C the standard
    int xnan = std::isnan(x), ynan = std::isnan(y);
    if(xnan || ynan) {
        if(xnan && !ynan) return y;
        if(!xnan && ynan) return x;
        return x;
    }
    return std::min(x,y);
}


/**
 * Calculate the `symmetric2` step pattern.
 *
 * @description
 * This function calculates the `symmetric2` step pattern, which uses a weight
 * of 2 for the diagonal step and 1 for the vertical and horizontal to
 * compensate for the favor of diagonal steps.
 *
 * @note
 * For more information on this step pattern, visit the `IncDTW` package
 * documentation: https://www.rdocumentation.org/packages/IncDTW/versions/1.1.4.4/topics/dtw2vec
 *
 * @reference
 * Leodolter, M., Plant, C., & Brändle, N. (2021). IncDTW: An R Package for
 * Incremental Calculation of Dynamic Time Warping. Journal of Statistical
 * Software, 99(9), 1–23. https://doi.org/10.18637/jss.v099.i09
 *
 * Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping
 * Alignments in R: The dtw Package. Journal of Statistical Software, 31(7),
 * 1–24. https://doi.org/10.18637/jss.v031.i07
 *
 * @return `symmetric2` step pattern value.
 */
double calculate_step_pattern_symmetric2(
    const double gcm10, // vertical
    const double gcm11, // diagonal
    const double gcm01, // horizontal
    const double cm00
) {
    return(cm00 + minval(gcm10, minval(cm00 + gcm11, gcm01)));
}


/**
 * Vector-based Dynamic Time Warping (DTW) distance.
 *
 * @description
 * This function calculates the Dynamic Time Warping (DTW) distance between
 * two sequences using the vector-based algorithm proposed by Leodolter
 * et al. (2021).
 *
 * The complexity of this function, as presented by Leodolter et al. (2021), is
 * equal to O(n).
 *
 * For more information on vector-based DTW, visit:
 * https://doi.org/10.18637/jss.v099.i09
 *
 * @param x A `arma::vec` with time-series values.
 * @param y A `arma::vec` with time-series values.
 *
 * @reference
 * Leodolter, M., Plant, C., & Brändle, N. (2021). IncDTW: An R Package for
 * Incremental Calculation of Dynamic Time Warping. Journal of Statistical
 * Software, 99(9), 1–23. https://doi.org/10.18637/jss.v099.i09
 *
 * @note
 * The implementation of this DTW distance calculation was adapted from the
 * `IncDTW` R package.
 *
 * @return DTW distance.
 */
// [[Rcpp::export]]
double dtw2vec(const arma::vec &x, const arma::vec &y)
{
    int nx = x.size();
    int ny = y.size();

    double *p1 = new double[nx];
    double *p2 = new double[nx];

    double *ptmp;
    double ret;

    // first column
    *p1 = std::abs(x[0] - y[0]);
    for (int i = 1; i < nx; i++)
    {
        p1[i] = std::abs(x[i] - y[0]) + p1[i - 1];
    }

    for (int j = 1; j < ny; j++)
    {
        *p2 = std::abs(x[0] - y[j]) + *(p1);

        for (int i = 1; i < nx; i++)
        {
            *(p2 + i) = calculate_step_pattern_symmetric2(*(p2 + i - 1), *(p1 + i - 1), *(p1 + i), std::abs(x[i] - y[j]));
        }
        ptmp = p1;
        p1 = p2;
        p2 = ptmp;
    }

    ret = *(p1 + nx - 1); // p1[nx-1]

    delete[] p1;
    delete[] p2;

    return (ret);
}


/**
 * Dynamic Time Warping (DTW) distance wrapper.
 *
 * @description
 * This function calculates prepare data from `Kohonen` package and calculate
 * the DTW distance between two 1D time-series.
 *
 * @param p1 A 1D array representing the first time-series.
 * @param p2 A 1D array representing the second time-series.
 * @param np Number of points in arrays `p1` and `p2`.
 * @param nNA Number of NA values in the arrays `p1` and `p2`.
 *
 * @note The function signature was created following the `Kohonen` R package
 *        specifications for custom distance functions.
 *
 * @return The DTW distance between two time-series.
 */
double kohonen_dtw(double *p1, double *p2, int np, int nNA)
{
    arma::vec p1_vec(p1, np, false);
    arma::vec p2_vec(p2, np, false);

    return dtw2vec(p1_vec, p2_vec);
}


// [[Rcpp::export]]
Rcpp::XPtr<DistanceFunctionPtr> dtw()
{
    // Returns a External Pointer, which is used by the `kohonen` package
    // https://cran.r-project.org/doc/manuals/R-exts.html#External-pointers-and-weak-references
    return (Rcpp::XPtr<DistanceFunctionPtr>(new DistanceFunctionPtr(
        &kohonen_dtw)));
}
