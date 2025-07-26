#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector RECUR(const IntegerMatrix& data, int target_class) {
    int nrow = data.nrow();
    int ncol = data.ncol();

    LogicalVector result(nrow, false);

    for (int i = 0; i < nrow; ++i) {
        std::vector<int> indices;

        // Step 1: Find positions of target_class in the row
        for (int t = 0; t < ncol; ++t) {
            if (data(i, t) == target_class) {
                indices.push_back(t);
            }
        }

        // Step 2: Look for gaps between consecutive appearances
        for (size_t k = 0; k + 1 < indices.size(); ++k) {
            int current = indices[k];
            int next = indices[k + 1];

            if (next - current > 1) {
                bool has_gap = false;
                for (int g = current + 1; g < next; ++g) {
                    if (data(i, g) != target_class) {
                        has_gap = true;
                        break;
                    }
                }

                if (has_gap) {
                    result[i] = true;
                    break;  // stop checking this row
                }
            }
        }
    }

    return result;
}

// [[Rcpp::export]]
LogicalVector CONVERT(const IntegerMatrix& data, int source_class, int target_class) {
    int nrow = data.nrow();
    int ncol = data.ncol();

    LogicalVector result(nrow, false);

    for (int i = 0; i < nrow; ++i) {
        for (int t = 0; t < ncol - 1; ++t) {
            if (data(i, t) == source_class && data(i, t + 1) == target_class) {
                result[i] = true;
                break;
            }
        }
    }

    return result;
}

// [[Rcpp::export]]
LogicalVector EVOLVE(const IntegerMatrix& data, int class_i, int class_j) {
    int nrow = data.nrow();
    int ncol = data.ncol();

    LogicalVector result(nrow, false);

    for (int i = 0; i < nrow; ++i) {
        int first_i = -1;
        int first_j = -1;

        for (int t = 0; t < ncol; ++t) {
            if (first_i != -1 && first_j != -1) {
                break;
            }

            if (data(i, t) == class_i && first_i == -1) {
                first_i = t;
            }

            if (data(i, t) == class_j && first_j == -1) {
                first_j = t;
            }
        }

        if (first_i != -1 && first_j != -1 && first_j > first_i) {
            result[i] = true;
        }
    }

    return result;
}

// [[Rcpp::export]]
LogicalVector KEEPS(const IntegerMatrix& data, int target_class) {
    int nrow = data.nrow();
    int ncol = data.ncol();

    LogicalVector result(nrow, false);

    for (int i = 0; i < nrow; ++i) {
        bool all_match = true;

        for (int j = 0; j < ncol; ++j) {
            if (data(i, j) != target_class) {
                all_match = false;
                break;
            }
        }

        result[i] = all_match;
    }

    return result;
}
