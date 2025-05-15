## check_2column_matrix

using("tinysnapshot")

opts <- options(cli.width = 80)

expect_error(phutil:::check_2column_matrix(0))
m <- as.matrix(noisy_circle_ripserr)
expect_error(phutil:::check_2column_matrix(m))
m <- m[, 2:3]
m1 <- m
m1[1, 1] <- "a"
expect_error(phutil:::check_2column_matrix(m1))
m1 <- m
m1[1, 1] <- NA
expect_error(phutil:::check_2column_matrix(m1))
m1 <- m
m1[1, 1] <- 1
m1[1, 2] <- 0
expect_false(phutil:::check_2column_matrix(m1))
expect_snapshot_print(
  phutil:::check_2column_matrix(m1),
  label = "print-check-2column-matrix-warn"
)
expect_snapshot_print(
  phutil:::check_2column_matrix(m1, warn = FALSE),
  label = "print-check-2column-matrix-nowarn"
)
m1 <- m
m1[1, 1] <- -1
expect_true(phutil:::check_2column_matrix(m1))

options(opts)
