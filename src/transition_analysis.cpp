#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix transition_neighbor_analysis(NumericMatrix data, int reference_class, int neighbor_class) {
    int npixel = data.nrow();
    int nyear = data.ncol();

    if (nyear != 10) {
        stop("Expected exactly 10 years (columns), but got " + std::to_string(nyear));
    }

    // 2015 - 0
    // 2016 - 1
    // 2017 - 2

    // 2018 - 3
    // 2019 - 4
    // 2020 - 5

    // 2021 - 6
    // 2022 - 7
    // 2023 - 8

    // 2024 - 9

    for (int i = 0; i < npixel; i++) {
        for (int j = 0; j + 2 < nyear; j += 3) {

            bool valid_neihbor = data(i, j) == neighbor_class && data(i, j + 2) == neighbor_class;
            bool valid_class = data(i, j + 1) == reference_class;

            if (valid_class && valid_neihbor) {
                data(i, j + 1) = neighbor_class;
            }

            // valid iterations (value of j)
            // 0
            // 3
            // 6
            // 9
        }

        // processing last pixel (out of the rules presented before)
        if (data(i, nyear - 1) == reference_class) {
            data(i, nyear - 1) = neighbor_class;
        }
    }

    return data;
}
