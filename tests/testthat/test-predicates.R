
test_that("Starts works as expected", {
    # Define data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, 2, 2,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, NA, 1, 0
    ), nrow = 6, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    2
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    NA    1    0

    # Tests
    # Test starts of class '0' (NA)
    res <- STARTS(values, 0)
    expect_equal(sum(res), 1)

    # Test starts of class '1'
    res <- STARTS(values, 1)
    expect_equal(sum(res), 3)

    # Test starts of class '2'
    res <- STARTS(values, 2)
    expect_equal(sum(res), 1)
})

test_that("Ends works as expected", {
    # Define data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, 2, NA,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 9
    ), nrow = 7, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    NA
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    1    1    0
    # [7,]    0    0    0    0    9

    # Test ends of class '0'
    res <- ENDS(values, 0)
    expect_equal(sum(res), 2)

    # Test ends of class '1'
    res <- ENDS(values, 1)
    expect_equal(sum(res), 2)

    # Test ends of class '2'
    res <- ENDS(values, 2)
    expect_equal(sum(res), 1)

    # Test ends of class '9'
    res <- ENDS(values, 9)
    expect_equal(sum(res), 1)
})

test_that("Edges works as expected", {
    # Define data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, 2, NA,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 9
    ), nrow = 7, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    NA
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    1    1    0
    # [7,]    0    0    0    0    9

    # Test ends of class '0'
    res <- EDGES(values, 0)
    expect_equal(sum(res), 2)

    # Test ends of class '1'
    res <- EDGES(values, 1)
    expect_equal(sum(res), 2)

    # Test ends of class '2' (NA)
    res <- EDGES(values, 2)
    expect_equal(sum(res), 0)

    # Test ends of class '9'
    res <- EDGES(values, 9)
    expect_equal(sum(res), 0)
})

test_that("Persist works as expected", {
    # Define data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, NA, 2,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, 1, 1, NA,
        0, 0, 0, 0, 9
    ), nrow = 7, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    2
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    1    1    NA
    # [7,]    0    0    0    0    9

    # Test persist of class '0'
    res <- PERSIST(values, 0, 1)
    expect_equal(sum(res), 3)

    res <- PERSIST(values, 0 , 5)
    expect_equal(sum(res), 1)

    res <- PERSIST(values, 0, 2)
    expect_equal(sum(res), 0)

    # Test persist of class '1' (NA)
    res <- PERSIST(values, 1, 3)
    expect_equal(sum(res), 2)

    res <- PERSIST(values, 1, 2)
    expect_equal(sum(res), 1)

    # Test persist of class '9'
    res <- PERSIST(values, 9, 1)
    expect_equal(sum(res), 1)

    # Test persist of invalid class
    res <- PERSIST(values, 99, 1)
    expect_equal(sum(res), 0)
})

test_that("Peaks works as expected", {
    # Define data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, 2, 2,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, NA, 0, 9
    ), nrow = 7, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    2
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    1    1    0
    # [7,]    0    0    NA    0    9

    # Test peaks of class '0'
    res <- PEAKS(values, 0)
    expect_true(res[1] == TRUE)
    expect_equal(sum(res), 4)

    # Test peaks of class '1'
    res <- PEAKS(values, 1)
    expect_equal(sum(res), 0)

    # Test peaks of class '9' (NA)
    res <- PEAKS(values, 9)
    expect_equal(sum(res), 0)
})

test_that("Recur works as expected", {
    # Define data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, 2, 2,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, NA, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 9
    ), nrow = 7, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    2
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    NA   0    1    1
    # [6,]    0    1    1    1    0
    # [7,]    0    0    0    0    9

    # Test recur of class '0'
    res <- RECUR(values, 0)
    expect_equal(sum(res), 1)

    # Test recur of class '1' (NA)
    res <- RECUR(values, 1)
    expect_equal(which(res), 3)
    expect_equal(sum(res), 1)

    # Test recur of class '2'
    res <- RECUR(values, 2)
    expect_equal(sum(res), 0)
})

test_that("Convert works as expected", {
    # Define data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, 2, 2,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 9,
        0, NA, 4, 1, 19
    ), nrow = 8, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    2
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    1    1    0
    # [7,]    0    0    0    0    9
    # [8,]    0    NA    4    1   19

    # Tests
    # Test convert of class '0'
    res <- CONVERT(values, 0, 2)
    expect_equal(sum(res), 1)

    res <- CONVERT(values, 0, 1)
    expect_equal(sum(res), 3)

    res <- CONVERT(values, 0, 9)
    expect_equal(sum(res), 1)

    res <- CONVERT(values, 0, 0)
    expect_equal(sum(res), 2)
    expect_equal(which(res), c(4, 7))

    res <- CONVERT(values, 0, 9)
    expect_equal(sum(res), 1)

    res <- CONVERT(values, 0, 19)
    expect_equal(sum(res), 0)

    # Test convert of class '1'
    res <- CONVERT(values, 1, 0)
    expect_equal(sum(res), 4)

    res <- CONVERT(values, 1, 19)
    expect_equal(sum(res), 0)

    res <- CONVERT(values, 1, 1)
    expect_equal(sum(res), 4)
})

test_that("Evolve works as expected", {
    # Define data
    values <- matrix(c(
        1, NA, 1, 0, 2,
        2, 2, 2, 2, 2,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 9,
        0, 4, 4, 1, 19
    ), nrow = 8, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    NA    1   0    2
    # [2,]    2    2    2    2    2
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    1    1    0
    # [7,]    0    0    0    0    9
    # [8,]    0    4    4    1   19

    # Tests
    # Test evolve of class '0' (NA)
    res <- EVOLVE(values, 0, 2)
    expect_equal(sum(res), 0)

    res <- EVOLVE(values, 0, 9)
    expect_equal(sum(res), 1)

    res <- EVOLVE(values, 0, 19)
    expect_equal(sum(res), 1)

    res <- EVOLVE(values, 9, 99)
    expect_equal(sum(res), 0)

    # Test evolve of class '1' (NA)
    res <- EVOLVE(values, 1, 2)
    expect_equal(sum(res), 0)

    res <- EVOLVE(values, 1, 19)
    expect_equal(sum(res), 1)

    # Test evolve of invalid class
    res <- EVOLVE(values, -99, 123)
    expect_equal(sum(res), 0)
})

test_that("Keeps works as expected", {
    # Define data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, 2, 2,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 9,
        5, NA, 5, 5, 5,
        5, 5, 5, 5, 5
    ), nrow = 9, byrow = TRUE)

    # Visualize data
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    2
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    1    1    0
    # [7,]    0    0    0    0    9
    # [8,]    5    NA    5    5    5
    # [9,]    5    5    5    5    5

    # Test keeps of class '0'
    res <- KEEPS(values, 0)
    expect_equal(sum(res), 1)
    expect_equal(which(res), 4)

    # Test keeps of class '1'
    res <- KEEPS(values, 1)
    expect_equal(sum(res), 0)

    # Test keeps of class '5' (NA)
    res <- KEEPS(values, 5)
    expect_equal(sum(res), 1)
    expect_equal(which(res), 9)
})

test_that("Expressions conversion as expected", {
    # Auxiliary data
    labels = c(
        "0" = "Class 0",
        "1" = "Class 1",
        "2" = "Class 2",
        "4" = "Class 4",
        "5" = "Class 5",
        "9" = "Class 9",
        "19" = "Class 19"
    )

    # Starts
    expr <- substitute(Starts("Class 0"))
    expected <- substitute(STARTS(values, 0L))

    expr <- .transitions_expand_expr(expr, labels)
    expect_true(identical(expr, expected))

    # Ends
    expr <- substitute(Ends("Class 9"))
    expected <- substitute(ENDS(values, 9L))

    expr <- .transitions_expand_expr(expr, labels)
    expect_true(identical(expr, expected))

    # Persist
    expr <- substitute(Persist("Class 5", 2))
    expected <- substitute(PERSIST(values, 5L, 2))

    expr <- .transitions_expand_expr(expr, labels)
    expect_true(identical(expr, expected))

    # Peaks
    expr <- substitute(Peaks("Class 2"))
    expected <- substitute(PEAKS(values, 2L))

    expr <- .transitions_expand_expr(expr, labels)
    expect_true(identical(expr, expected))

    # Recur
    expr <- substitute(Recur("Class 4"))
    expected <- substitute(RECUR(values, 4L))

    expr <- .transitions_expand_expr(expr, labels)
    expect_true(identical(expr, expected))

    # Convert
    expr <- substitute(Convert("Class 4", "Class 5"))
    expected <- substitute(CONVERT(values, 4L, 5L))

    expr <- .transitions_expand_expr(expr, labels)
    expect_true(identical(expr, expected))

    # Evolve
    expr <- substitute(Evolve("Class 9", "Class 19"))
    expected <- substitute(EVOLVE(values, 9L, 19L))

    expr <- .transitions_expand_expr(expr, labels)
    expect_true(identical(expr, expected))

    # Keeps
    expr <- substitute(Keeps("Class 2"))
    expected <- substitute(KEEPS(values, 2L))

    expr <- .transitions_expand_expr(expr, labels)
    expect_true(identical(expr, expected))
})

test_that("Expressions logics works as expected", {
    # Auxiliary data
    labels = c(
        "0" = "Class 0",
        "1" = "Class 1",
        "2" = "Class 2",
        "4" = "Class 4",
        "5" = "Class 5",
        "9" = "Class 9",
        "19" = "Class 19"
    )

    # Conversion 1 - Unique expression
    expr1 <- substitute(Starts("Class 0"))
    expected <- substitute(STARTS(values, 0L))

    expr1 <- .transitions_expand_expr(expr1, labels)
    expect_true(identical(expr1, expected))

    # Conversion 2 - not
    expr2 <- substitute(!Starts("Class 0"))
    expected2 <- substitute(!STARTS(values, 0L))

    expr2 <- .transitions_expand_expr(expr2, labels)
    expect_true(identical(expr2, expected2))

    # Conversion 3 - and
    expr3 <- substitute(Starts("Class 0") & Starts("Class 1"))
    expected3 <- substitute(STARTS(values, 0L) & STARTS(values, 1L))

    expr3 <- .transitions_expand_expr(expr3, labels)
    expect_true(identical(expr3, expected3))

    # Conversion 4 - or
    expr4 <- substitute(Starts("Class 0") | Starts("Class 1"))
    expected4 <- substitute(STARTS(values, 0L) | STARTS(values, 1L))

    expr4 <- .transitions_expand_expr(expr4, labels)
    expect_true(identical(expr4, expected4))

    # Conversion 5 - and + or
    expr5 <- substitute(Starts("Class 0") | Starts("Class 1") & Ends("Class 19"))
    expected5 <- substitute(STARTS(values, 0L) | STARTS(values, 1L) & ENDS(values, 19L))

    expr5 <- .transitions_expand_expr(expr5, labels)
    expect_true(identical(expr5, expected5))

    # Conversion 6 - and + or + not
    expr6 <- substitute(
        Starts("Class 0") |
            Starts("Class 1") &
            Ends("Class 19") &
            !Peaks("Class 9")
    )
    expected6 <- substitute(
        STARTS(values, 0L) |
            STARTS(values, 1L) &
            ENDS(values, 19L) &
            !PEAKS(values, 9L)
    )

    expr6 <- .transitions_expand_expr(expr6, labels)
    expect_true(identical(expr6, expected6))
})

test_that("Expressions call works as expected", {
    # Auxiliary data
    values <- matrix(c(
        1, 1, 1, 0, 2,
        2, 2, 2, 2, 2,
        1, 0, 1, 1, 1,
        0, 0, 0, 0, 0,
        1, 1, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 9,
        0, 4, 4, 1, 19,
        5, 5, NA, 5, 5,
        5, 5, 5, 5, 5
    ), nrow = 10, byrow = TRUE)

    labels = c(
        "0" = "Class 0",
        "1" = "Class 1",
        "2" = "Class 2",
        "4" = "Class 4",
        "5" = "Class 5",
        "9" = "Class 9",
        "19" = "Class 19"
    )

    # Data visualization
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    1    0    2
    # [2,]    2    2    2    2    2
    # [3,]    1    0    1    1    1
    # [4,]    0    0    0    0    0
    # [5,]    1    1    0    1    1
    # [6,]    0    1    1    1    0
    # [7,]    0    0    0    0    9
    # [8,]    0    4    4    1   19
    # [9,]    5    5    NA    5    5
    # [10,]   5    5    5    5    5

    # Define evaluation environment
    env <- list2env(list(
        values = values
    ))

    # Expression 1 - Starts
    expr <- substitute(Starts("Class 0"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 4)

    expr <- substitute(Starts("Class 0") | Starts("Class 5"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 5)

    # Expression 2 - Ends
    expr <- substitute(Ends("Class 19"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 1)

    # Expression 3 - Persist
    expr <- substitute(Persist("Class 4", 2))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 1)

    expr <- substitute(Persist("Class 4", 2) | Persist("Class 1", 2))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 2)

    # Expression 4 - Peaks
    expr <- substitute(Peaks("Class 9"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 1)

    expr <- substitute(Peaks("Class 1") & Peaks("Class 19"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 1)

    # Expression 5 - Recur
    expr <- substitute(Recur("Class 1"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 2)

    expr <- substitute(Recur("Class 5"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 0)

    # Expression 6 - Convert
    expr <- substitute(Convert("Class 0", "Class 1"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 3)

    # Using more () to test transformation as well
    expr <- substitute(((Convert("Class 0", "Class 1")) | (Convert("Class 1", "Class 19"))))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 4)

    # Expression 7 - Evolve
    expr <- substitute(Evolve("Class 0", "Class 19"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 1)

    expr <- substitute(Evolve("Class 4", "Class 19") | Evolve("Class 1", "Class 2"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 2)

    expr <- substitute(Evolve("Class 4", "Class 19") & Evolve("Class 1", "Class 2"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 0)

    # Expression 7 - Keeps
    expr <- substitute(Keeps("Class 0"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 1)

    expr <- substitute(Keeps("Class 0") | Keeps("Class 5"))
    expr <- .transitions_expand_expr(expr, labels)
    res <- eval(expr, envir = env)

    expect_equal(sum(res), 2)
})
