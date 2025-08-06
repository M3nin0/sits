library(sits)

generate_values <- function() {
    matrix(c(
        1, 1, 1, 0, 2, 0, 0, 0, 0, 0,
        2, 2, 12, 2, 2, 0, 2, 12, 2, 0,
        1, 0, 1, 1, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        2, 12, 2, 2, 1, 2, 0, 0, 0, 0,
        2, 1, 2, 2, 1, 2, 2, 1, 2, 1,
        0, 1, NA, 1, 0, 2, 12, 2, 0, 1,
        2, 1, 2, 2, 1, 2, 0, 0, 0, 12,
        2, 1, 2, 2, 1, 2, 2, 12, 2, 12
    ), nrow = 9, byrow = TRUE)
}

# Define data
values <- generate_values()
values_original <- generate_values()

# > values_ref
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    1    1    1    0    2    0    0    0    0     0
# [2,]    2    2    12   2    2    0    2   12    2     0 # <---- (n = 1)
# [3,]    1    0    1    1    1    0    0    0    0     0
# [4,]    0    0    0    0    0    0    0    0    0     0
# [5,]    2   12    2    2    1    2    0    0    0     0 # <---- (n = 1)
# [6,]    2    1    2    2    1    2    2    1    2     1
# [7,]    0    1   NA    1    0    2   12    2    0     1 # <---- (n = 0)
# [8,]    2    1    2    2    1    2    0    0    0    12 # <---- (n = 0)
# [9,]    2    1    2    2    1    2    2   12    2    12 # <---- (n = 1 + change in the last one)

sits:::transition_neighbor_analysis(
    data = values,
    reference_class = 12,  # 12 = "vegetacao_secundaria"
    neighbor_class = 2 # 2 = "Ag_perene"
)

# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    1    1    1    0    2    0    0    0    0     0
# [2,]    2    2   12    2    2    0    2    2    2     0 # <---- (n = 1)
# [3,]    1    0    1    1    1    0    0    0    0     0
# [4,]    0    0    0    0    0    0    0    0    0     0
# [5,]    2    2    2    2    1    2    0    0    0     0 # <---- (n = 1)
# [6,]    2    1    2    2    1    2    2    1    2     1
# [7,]    0    1   NA    1    0    2   12    2    0     1 # <---- (n = 0)
# [8,]    2    1    2    2    1    2    0    0    0    12 # <---- (n = 0)
# [9,]    2    1    2    2    1    2    2    2    2     2 # <---- (n = 1 + change in the last one)

values == values_original
# > values == values_original
# [,1]  [,2] [,3] [,4] [,5] [,6] [,7]  [,8] [,9] [,10]
# [1,] TRUE  TRUE TRUE TRUE TRUE TRUE TRUE  TRUE TRUE  TRUE
# [2,] TRUE  TRUE TRUE TRUE TRUE TRUE TRUE FALSE TRUE  TRUE
# [3,] TRUE  TRUE TRUE TRUE TRUE TRUE TRUE  TRUE TRUE  TRUE
# [4,] TRUE  TRUE TRUE TRUE TRUE TRUE TRUE  TRUE TRUE  TRUE
# [5,] TRUE FALSE TRUE TRUE TRUE TRUE TRUE  TRUE TRUE  TRUE
# [6,] TRUE  TRUE TRUE TRUE TRUE TRUE TRUE  TRUE TRUE  TRUE
# [7,] TRUE  TRUE   NA TRUE TRUE TRUE TRUE  TRUE TRUE  TRUE
# [8,] TRUE  TRUE TRUE TRUE TRUE TRUE TRUE  TRUE TRUE  TRUE
# [9,] TRUE  TRUE TRUE TRUE TRUE TRUE TRUE FALSE TRUE FALSE

which(!values == values_original)
# > which(!values == values_original)
# [1] 14 65 72 90
