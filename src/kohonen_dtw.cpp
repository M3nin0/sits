#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "./sits_types.h"

extern "C"
{
#include "./kohonen_dtw_distance.h"
}

using namespace Rcpp;


/**
 * Create a `symmetric2` step pattern matrix.
 *
 * @description
 * This function calculates the `symmetric2` step pattern matrix
 * as defined in the `dtw` package.
 *
 * @note
 * For more information on this step pattern, visit the `dtw` package
 * documentation: https://www.rdocumentation.org/packages/dtw/versions/1.23-1/topics/stepPattern
 *
 * @reference
 * Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping
 * Alignments in R: The dtw Package. Journal of Statistical Software, 31(7),
 * 1–24. https://doi.org/10.18637/jss.v031.i07
 *
 * @return NumericMatrix of the `symmetric2` step pattern.
 */
// [[Rcpp::export]]
arma::mat step_pattern_factory_symmetric2()
{
    static arma::mat m;

    if (m.n_elem == 0)
    {                     // Populates matrix with `symmetric2` step pattern.
        m.set_size(6, 4);
        m.fill(arma::fill::zeros);

        // Values extracted from `dtw` package.
        m.row(0) = arma::vec({1, 1, 1, -1}).t();
        m.row(1) = arma::vec({1, 0, 0, 2}).t();
        m.row(2) = arma::vec({2, 0, 1, -1}).t();
        m.row(3) = arma::vec({2, 0, 0, 1}).t();
        m.row(4) = arma::vec({3, 1, 0, -1}).t();
        m.row(5) = arma::vec({3, 0, 0, 1}).t();
    }

    return m;
}


/**
 * Create a local cost matrix using `Euclidean` distance.
 *
 * @description
 * This function calculates the cross-distance matrix between two vectors.
 *
 * @reference
 * Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping
 * Alignments in R: The dtw Package. Journal of Statistical Software, 31(7),
 * 1–24. https://doi.org/10.18637/jss.v031.i07
 *
 * @return NumericMatrix with a local cost between two vectors.
 */
// [[Rcpp::export]]
arma::mat local_cost_matrix(arma::mat mat_a, arma::mat mat_b)
{
    // Create column and row vectors for broadcasting
    arma::mat A_mat = arma::repmat(mat_a, 1, mat_a.n_elem);
    arma::mat B_mat = arma::repmat(mat_b.t(), mat_b.n_elem, 1);

    // Compute the absolute differences using matrix operations
    arma::mat local_cost = arma::abs(A_mat - B_mat);

    return local_cost;
}


/**
 * Create `No Window` windowing function
 *
 * @description
 * This function implements the windowing function `NoWindow` presented on the
 * `dtw` package.
 *
 * @reference
 * Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping
 * Alignments in R: The dtw Package. Journal of Statistical Software, 31(7),
 * 1–24. https://doi.org/10.18637/jss.v031.i07
 *
 * @return LogicalMatrix representing the `No Window` function.
 */
// [[Rcpp::export]]
arma::imat windowing_factory_no_window(int n)
{
    static arma::imat m;

    if (m.n_elem == 0)
    {
        m.set_size(n, n);
        m.fill(arma::fill::ones);
    }

    return m;
}


/**
 * Create a empty cost matrix.
 *
 * @description
 * This function creates a cost matrix based on a local cost matrix.
 *
 * @reference
 * Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping
 * Alignments in R: The dtw Package. Journal of Statistical Software, 31(7),
 * 1–24. https://doi.org/10.18637/jss.v031.i07
 *
 * @return NumericMatrix with NA values.
 */
// [[Rcpp::export]]
arma::mat initialize_cost_matrix(arma::mat local_cost)
{
    int n = local_cost.n_rows;

    // Initialize cm with NaN values
    arma::mat cm(n, n);
    cm.fill(NA_REAL);

    // Assign value of local_cost(0, 0) to cm(0, 0)
    cm(0, 0) = local_cost(0, 0);

    return cm;
}


/**
 * Compute the Dynamic Time Warping (DTW) distance between two 2D C++ vectors.
 *
 * @description
 * This function calculates the Dynamic Time Warping (DTW) distance between
 * two sequences that can have a different number of data points but must
 * share the same number of dimensions. An exception is thrown if the dimensions
 * of the input vectors do not match.
 *
 * For more information on DTW, visit:
 * https://en.wikipedia.org/wiki/Dynamic_time_warping
 *
 * @param wm A `LogicalVector` representing the restriction windowing function.
 * @param lm A `NumericMatrix` with the cross-distances between two time-series.
 * @param cm A `NumericMatrix` representing a global cost matrix.
 * @param dir A `NumericMatrix` representing the step pattern applied.
 *
 * @reference
 * Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping
 * Alignments in R: The dtw Package. Journal of Statistical Software, 31(7),
 * 1–24. https://doi.org/10.18637/jss.v031.i07
 *
 * @note
 * The implementation of this DTW distance calculation was adapted from the
 * `dtw` R package.
 *
 * @return The DTW distance.
 */
// [[Rcpp::export]]
double distance_dtw_op(arma::imat wm, arma::mat lm, arma::mat cm, arma::mat dir)
{
    // Get problem size
    int n = lm.n_rows;

    // Get pattern size
    int nsteps = dir.n_rows;

    // Cost matrix (input + output)
    arma::mat cmo = cm;

    // sits: This output is not used, so we disable it's calculation
    // Output 2: smo, INTEGER
    // arma::mat smo(n, n);

    // Dispatch to C
    computeCM(
        n,
        wm.memptr(),
        lm.memptr(),
        &nsteps,
        dir.memptr(),
        cmo.memptr()
    );

    return cmo.at(n - 1, n - 1);
}


/**
 * Dynamic Time Warping (DTW) distance wrapper.
 *
 * @description
 * This function calculates prepare data from `Kohonen` package and calculate
 * the DTW distance between two array of points.
 *
 * @param p1 A 1D array representing the first time-series
 * @param p2 A 1D array representing the second sequence.
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
    arma::vec p1_mat(p1, np, false);
    arma::vec p2_mat(p2, np, false);

    // Local cost matrix
    arma::mat lc_matrix = local_cost_matrix(p1_mat, p2_mat);

    // Step pattern matrix
    arma::mat step_pattern = step_pattern_factory_symmetric2();

    // Window function (in this case `No Window`)
    arma::imat no_window = windowing_factory_no_window(np);

    // Cost matrix
    arma::mat cost_matrix = initialize_cost_matrix(lc_matrix);

    return (distance_dtw_op(
        no_window, lc_matrix, cost_matrix, step_pattern));
}


// [[Rcpp::export]]
Rcpp::XPtr<DistanceFunctionPtr> dtw()
{
    // Returns a External Pointer, which is used by the `kohonen` package
    // https://cran.r-project.org/doc/manuals/R-exts.html#External-pointers-and-weak-references
    return (Rcpp::XPtr<DistanceFunctionPtr>(new DistanceFunctionPtr(
        &kohonen_dtw)));
}
