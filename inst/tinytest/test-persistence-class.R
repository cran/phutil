using("tinysnapshot")

opts <- options(cli.width = 80)

m <- as.matrix(noisy_circle_ripserr)

mat <- as_persistence(m)
expect_equal(mat$metadata$data, "?")
expect_equal(mat$metadata$engine, "?")
expect_equal(mat$metadata$filtration, "?")
expect_equal(mat$metadata$call, "?")
expect_equal(mat$metadata$parameters, list())
expect_snapshot_print(mat, label = "print-mat-persistence")

params <- list(maxdimension = 1L, maxscale = 1.6322)
mat <- as_persistence(m, engine = "ripserr::vietoris_rips", parameters = params)
expect_equal(mat$metadata$engine, "ripserr::vietoris_rips")
expect_equal(mat$metadata$parameters, params)
expect_snapshot_print(mat, label = "print-mat-persistence-engine")

mat <- as_persistence(m, filtration = "vietoris-rips")
expect_equal(mat$metadata$filtration, "vietoris-rips")
expect_snapshot_print(mat, label = "print-mat-persistence-filtration")

more_df <- noisy_circle_ripserr
more_df$morevar <- 1L
expect_inherits(as_persistence(more_df), "persistence")

df_2cols <- noisy_circle_ripserr
df_2cols$dimension <- NULL
ps <- as_persistence(df_2cols)
expect_inherits(ps, "persistence")

wrong_df <- noisy_circle_ripserr
colnames(wrong_df) <- c("x", "y", "z")
expect_error(as_persistence(wrong_df))

ripserr <- as_persistence(noisy_circle_ripserr)
expect_equal(sapply(ripserr$pairs, ncol), c(2L, 2L))
expect_equal(ripserr$metadata$ordered_pairs, TRUE)
expect_equal(as.character(ripserr$metadata$data), "?")
expect_equal(ripserr$metadata$engine, "ripserr::<vietoris_rips/cubical>")
expect_equal(ripserr$metadata$filtration, "Vietoris-Rips/cubical")
expect_equal(as.character(ripserr$metadata$call), "?")
expect_equal(length(ripserr$metadata$parameters), 0L)
expect_snapshot_print(ripserr, label = "print-ripserr-persistence")

tda_rips <- as_persistence(getElement(noisy_circle_tda_rips, "diagram"))
expect_equal(sapply(tda_rips$pairs, ncol), c(2L, 2L))

tda_rips <- as_persistence(noisy_circle_tda_rips)
expect_equal(sapply(tda_rips$pairs, ncol), c(2L, 2L))
expect_equal(as.character(tda_rips$metadata$data), "noisy_circle_points")
expect_equal(tda_rips$metadata$engine, "TDA::ripsDiag")
expect_equal(tda_rips$metadata$filtration, "Vietoris-Rips")
expect_equal(
  as.character(tda_rips$metadata$call),
  c("TDA::ripsDiag", "noisy_circle_points", "1", "1.6322")
)
expect_equal(length(tda_rips$metadata$parameters), 2L)
expect_equal(names(tda_rips$metadata$parameters), c("maxdimension", "maxscale"))
expect_equal(tda_rips$metadata$parameters$maxdimension, 1L)
expect_equal(tda_rips$metadata$parameters$maxscale, 1.6322)
expect_snapshot_print(
  as_persistence(tda_rips),
  label = "print-tda-rips-persistence"
)

expect_equal(ncol(get_pairs(tda_rips, dimension = 1)), 2L)
expect_error(get_pairs(m, dimension = 1))
expect_equal(nrow(get_pairs(tda_rips, dimension = 2)), 0L)
expect_equal(ncol(get_pairs(tda_rips, dimension = 2)), 2L)

expect_equal(ncol(as.data.frame(tda_rips)), 3L)
row_names <- c(paste("0", 1:100, sep = "-"), paste("1", 1:2, sep = "-"))
expect_equal(
  rownames(as.data.frame(tda_rips, optional = FALSE)),
  row_names
)
row_names <- paste0("row", 1:102)
expect_equal(
  rownames(as.data.frame(tda_rips, row.names = row_names)),
  row_names
)

d <- dist(cbind(x = c(0, 3, 0), y = c(0, 0, 4)))
hc <- hclust(d, method = "complete")
expect_error(as_persistence(hc, birth = 6))
hclust <- as_persistence(hc)
expect_equal(hclust$pairs[[1]][, 1], rep(0, 3L))
expect_snapshot_print(hclust, label = "print-hclust-persistence")
# negative distances
d2 <- d
d2[2] <- -d2[2]
hc2 <- hclust(d2, method = "single")
hclust2 <- as_persistence(hc2)
expect_equal(hclust2$pairs[[1]][, 1], rep(-Inf, 3L))

# Test that persistence is correctly structured when some degrees are missing
x <- data.frame(
  dimension = c(0L, 3L, 0L, 2L),
  birth = c(0.0, 0.5, 0.2, 0.7),
  death = c(0.5, 1.0, 0.3, 1.2),
  morevar = c(1L, 2L, 3L, 4L)
)
xp <- as_persistence(x)
expect_length(xp$pairs, 4L)
expect_equal(sapply(xp$pairs, nrow), c(2L, 0L, 1L, 1L))
expect_equal(sapply(xp$pairs, ncol), rep(2L, 4L))

# Test that persistence complains when some "dimension" values are not whole
# numbers
y <- data.frame(
  dimension = c(0L, 3L, NA, 2L),
  birth = c(0.0, 0.5, 0.2, 0.7),
  death = c(0.5, 1.0, 0.3, 1.2),
  morevar = c(1L, 2L, 3L, 4L)
)
expect_message(
  as_persistence(y),
  pattern = "Negative, infinite, and missing dimensions will be omitted."
)

# Test that as_persistence() errors out if provided with a matrix with
# less than 2 columns
x <- matrix(1:6, ncol = 1)
expect_error(
  as_persistence(x),
  pattern = "The matrix must have at least 2 columns."
)

# Test that as_persistence() correctly processes unnamed matrices
x <- matrix(1:6, ncol = 3)
expect_inherits(as_persistence(x), "persistence")

# Test that as.matrix.persistence() works

## Test that as.matrix.persistence() return 0-row matrix if no pairs are
## present
x <- data.frame(
  dimension = integer(),
  birth = numeric(),
  death = numeric()
)
xp <- as_persistence(x)
expect_snapshot_print(xp, label = "print-empty-persistence")
expect_inherits(as.matrix(xp), "matrix")
expect_equal(nrow(as.matrix(xp)), 0L)
expect_equal(ncol(as.matrix(xp)), 3L)
expect_equal(colnames(as.matrix(xp)), c("dimension", "birth", "death"))

## Test that as.matrix.persistence() returns a matrix with the correct number of
## rows and columns
x <- as_persistence(list(
  matrix(c(0, 1, 2, 3), ncol = 2),
  matrix(c(0, 1, 2, 3), ncol = 2)
))
xm <- as.matrix(x)
expect_inherits(xm, "matrix")
expect_equal(xm[, "dimension"], c(0L, 0L, 1L, 1L))
expect_equal(xm[, "birth"], c(0, 1, 0, 1))
expect_equal(xm[, "death"], c(2, 3, 2, 3))

options(opts)
