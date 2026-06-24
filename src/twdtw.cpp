#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "./twdtw.h"

using namespace Rcpp;

/**
 * Compute the p-norm between two time-series.
 *
 * @description
 * The `p-norm`, also known as the `Minkowski space`, is a generalized norm
 * calculation that includes several types of distances based on the value
 * of `p`.
 *
 * Common values of `p` include:
 *
 *  - `p = 1` for the Manhattan (city block) distance;
 *  - `p = 2` for the Euclidean norm (distance).
 *
 * More details about p-norms can be found on Wikipedia:
 * https://en.wikipedia.org/wiki/Norm_(mathematics)#p-norm
 *
 * @param a A `std::vector<double>` with time-series values.
 * @param b A `std::vector<double>` with time-series values.
 * @param p A `double` value of the norm to use, determining the type of
 *          distance calculated.
 *
 * @note
 * Both vectors `a` and `b` must have the same length.
 *
 * @note
 * The implementation of this DTW distance calculation was adapted from the
 * `DTW_cpp` single header library (https://github.com/cjekel/DTW_cpp).
 *
 * @return The `p-norm` value between vectors `a` and `b`.
 */
double dist_p_norm(std::vector<double> a, std::vector<double> b, double p)
{
    double d = 0;

    size_t index;
    size_t a_size = a.size();

    for (index = 0; index < a_size; index++)
    {
        d += std::pow(std::abs(a[index] - b[index]), p);
    }
    return std::pow(d, 1.0 / p);
}

/**
 * Time-weight applied to the local cost of the TWDTW alignment.
 *
 * @description
 * Penalizes the matching of points that are distant in time, using the
 * time-weight functions of Maus et al. (2016, 2019):
 *
 *  - `logistic`: omega = 1 / (1 + exp(-alpha * (g - beta)))
 *  - `linear`:   omega = a * g + b
 *  - `none`:     omega = 0 (plain DTW)
 *
 * where `g` is the (circular) difference in days between two dates.
 *
 * NOTE on the sign of `alpha`: with this (negative-exponent) logistic form,
 * `omega` must increase with `g` so that matching observations far apart in
 * time is penalised. This requires a positive steepness (default alpha = 0.1),
 * matching the `twdtw` package values.
 *
 * @reference
 * Maus, V., Camara, G., Appel, M., & Pebesma, E. (2019). dtwSat:
 * Time-Weighted Dynamic Time Warping for Satellite Image Time Series
 * Analysis in R. Journal of Statistical Software, 88(5), 1-31.
 */
static inline double twdtw_time_weight(double g, int weight_type,
                                       double alpha, double beta,
                                       double a, double b)
{
    // logistic
    if (weight_type == 1) {
        return 1.0 / (1.0 + std::exp(-alpha * (g - beta)));
    }

    // linear
    if (weight_type == 2) {
        return a * g + b;
    }

    // none
    return 0.0;
}

/**
 * Circular difference, in days, between two day-of-year values.
 *
 * @description
 * The elapsed time `g(t1, t2)` of Maus et al. (2016) is the absolute
 * difference in days, taken on a yearly cycle, so the maximum distance
 * between two days is half a year (e.g. Jan 1st and Dec 15th are close).
 */
static inline double twdtw_doy_diff(double d1, double d2)
{
    double diff = std::abs(d1 - d2);

    if (diff > 182.5) {
        diff = 365.0 - diff;
    }

    return diff;
}

double twdtw_distance_op(const std::vector<std::vector<double>> &q,
                         const std::vector<std::vector<double>> &p,
                         const std::vector<double> &q_doy,
                         const std::vector<double> &p_doy,
                         double dist_power,
                         int weight_type,
                         double alpha, double beta,
                         double a, double b)
{
    int n = q.size();
    int m = p.size();

    // Accumulated cost matrix
    std::vector<std::vector<double>> d(n, std::vector<double>(m, 0.0));

    // Local cost: Minkowski distance over bands + time weight
    auto local_cost = [&](int i, int j) {
        // phi
        double phi = dist_p_norm(q[i], p[j], dist_power);

        // g
        double g = twdtw_doy_diff(q_doy[i], p_doy[j]);

        // return!
        return phi + twdtw_time_weight(g, weight_type, alpha, beta, a, b);
    };

    // dtw
    d[0][0] = local_cost(0, 0);

    for (int i = 1; i < n; i++) {
        d[i][0] = d[i - 1][0] + local_cost(i, 0);
    }

    for (int j = 1; j < m; j++) {
        d[0][j] = d[0][j - 1] + local_cost(0, j);
    }

    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            double c = local_cost(i, j);

            // symmetric2 step from https://doi.org/10.18637/jss.v031.i07
            d[i][j] = std::min({
                d[i - 1][j - 1] + 2.0 * c,
                d[i    ][j - 1] +       c,
                d[i - 1][j    ] +       c
            });
        }
    }

    // symmetric2 path-length normalization, so distances are comparable
    // across patterns of different length
    return d[n - 1][m - 1] / static_cast<double>(n + m);
}

/**
 * Compute TWDTW distances between query time-series and class patterns.
 *
 * @description
 * For each query time-series (a row of `values`) and each class pattern,
 * computes the TWDTW distance. The query rows follow the `sits` predictor
 * layout (band-major, time-minor), i.e. the value of band `b` (0-based) at
 * time `t` (0-based) is at column `b * n_times + t`.
 *
 * @param values      Query predictors.
 * @param patterns    List of pattern matrices, each `m_k x n_bands`.
 * @param query_doy   Day-of-year of each query time point (length `n_times`).
 * @param pattern_doy List of day-of-year vectors, one per pattern.
 * @param n_bands     Number of bands.
 * @param dist_power  `p` value of the Minkowski local distance.
 * @param weight_type 0 = none, 1 = logistic, 2 = linear.
 * @param alpha, beta Logistic time-weight parameters.
 * @param a, b        Linear time-weight parameters.
 *
 * @return Matrix `n_samples x n_patterns` of TWDTW distances.
 */
// [[Rcpp::export]]
NumericMatrix C_twdtw_distances(const NumericMatrix &values,
                                const List &patterns,
                                const NumericVector &query_doy,
                                const List &pattern_doy,
                                int n_bands,
                                double dist_power,
                                int weight_type,
                                double alpha, double beta,
                                double a, double b)
{
    int n_samples = values.nrow();
    int n_features = values.ncol();
    int n_times = n_features / n_bands;
    int n_patterns = patterns.size();

    // Pre-extract patterns and their day-of-year as C++ structures
    std::vector<std::vector<std::vector<double>>> pat_mats(n_patterns);
    std::vector<std::vector<double>> pat_doys(n_patterns);

    for (int k = 0; k < n_patterns; k++) {
        NumericMatrix pk = patterns[k]; // m_k x n_bands

        int m = pk.nrow();
        std::vector<std::vector<double>> pm(m, std::vector<double>(n_bands));

        for (int t = 0; t < m; t++) {
            for (int bnd = 0; bnd < n_bands; bnd++) {
                pm[t][bnd] = pk(t, bnd);
            }
        }

        pat_mats[k] = pm;
        pat_doys[k] = as<std::vector<double>>(pattern_doy[k]);
    }

    std::vector<double> q_doy = as<std::vector<double>>(query_doy);

    // Define output
    NumericMatrix out(n_samples, n_patterns);

    // Calculate TWDTW of all samples
    for (int s = 0; s < n_samples; s++) {
        // Rebuild the query as an (n_times x n_bands) matrix
        std::vector<std::vector<double>> q(n_times,
                                           std::vector<double>(n_bands));

        for (int t = 0; t < n_times; t++) {
            for (int bnd = 0; bnd < n_bands; bnd++) {
                q[t][bnd] = values(s, bnd * n_times + t);
            }
        }

        for (int k = 0; k < n_patterns; k++) {
            out(s, k) = twdtw_distance_op(
                q, pat_mats[k], q_doy, pat_doys[k],
                dist_power, weight_type, alpha, beta, a, b);
        }
    }

    // Return!
    return out;
}
