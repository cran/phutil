mat <- matrix(c(0, 1, 2, 1, 0.5, 2, 2, 1, 3), ncol = 3, byrow = TRUE)
colnames(mat) <- c("dimension", "birth", "death")
pd <- as_persistence(mat)
df <- data.frame(
  dimension = c(1, 1, 1, 1, 1),
  birth = c(1, .6, 1.3, 2.1, 1.99),
  death = c(2, 1.2, .7, 4.2, 1.99)
)


#as_persistence.default correctly processes a basic matrix"
expect_inherits(pd, "persistence")
expect_equal(length(pd$pairs), 3) # 3 unique degrees (0,1,2)
expect_equal(
  unname(pd$pairs[[1]]),
  matrix(c(1, 2), ncol = 2),
  ignore_attr = TRUE
) # Degree 0
expect_equal(
  unname(pd$pairs[[2]]),
  matrix(c(0.5, 2), ncol = 2),
  ignore_attr = TRUE
) # Degree 1
expect_equal(
  unname(pd$pairs[[3]]),
  matrix(c(1, 3), ncol = 2),
  ignore_attr = TRUE
) # Degree 2

#correctly processes output from `____$diag` from the {TDA} package
coords <- df[, c("birth", "death")]

df3d <- t(data.frame(
  w = c(1.1, -.2, 1.07),
  x = c(1.3, 2.1, 1.99),
  y = c(0.7, 4.2, 1.99),
  z = c(0.2, 0.4, 0.35)
))

pd <- TDA::alphaComplexDiag(df3d, maxdimension = 1)$diagram
pd[, c(2, 3)] <- sqrt(pd[, c(2, 3)])
pd_p <- as_persistence(pd)

expect_inherits(pd_p, "persistence")
expect_equal(pd_p$metadata$parameters$maxdimension, 1)
expect_equal(length(pd_p$pairs[[1]][1, ]), 2) #correct Format


pd2 <- TDA::alphaShapeDiag(df3d, maxdimension = 2)$diagram
pd2[, c(2, 3)] <- sqrt(pd2[, c(2, 3)])
pd_p2 <- as_persistence(pd2)

expect_inherits(pd_p2, "persistence")
expect_equal(pd_p$metadata$parameters$maxdimension, 1)
expect_equal(length(pd_p2$pairs[[1]][1, ]), 2) #correct Format


FltRips <- TDA::ripsFiltration(
  X = df,
  maxdimension = 1,
  maxscale = 1.5,
  dist = "euclidean",
  library = "Dionysus",
  printProgress = FALSE
)
DiagFltRips <- TDA::filtrationDiag(
  filtration = FltRips,
  maxdimension = 1,
  library = "Dionysus",
  location = TRUE,
  printProgress = FALSE
)
pd3 <- DiagFltRips$diagram
pd3[, c(2, 3)] <- sqrt(pd3[, c(2, 3)])
pd_p3 <- as_persistence(pd3)

expect_inherits(pd_p3, "persistence")
expect_equal(pd_p3$metadata$parameters$maxdimension, 1)
expect_equal(length(pd_p3$pairs[[1]][1, ]), 2) #correct Format


Diag1 <- TDA::gridDiag(
  coords,
  TDA::distFct,
  lim = cbind(c(-1, 1), c(-1, 1)),
  maxdimension = 1,
  by = 0.05,
  sublevel = TRUE,
  printProgress = FALSE
)
pd4 <- Diag1$diagram
pd4[, c(2, 3)] <- sqrt(pd4[, c(2, 3)])
pd_p4 <- as_persistence(pd4)

expect_inherits(pd_p4, "persistence")
expect_equal(pd_p4$metadata$parameters$maxdimension, 1)
expect_equal(length(pd_p4$pairs[[1]][1, ]), 2) #correct Format


pd5 <- TDA::ripsDiag(df, maxdimension = 1, maxscale = 10)$diagram
pd5[, c(2, 3)] <- sqrt(pd5[, c(2, 3)])
pd_p5 <- as_persistence(pd5)

expect_inherits(pd_p5, "persistence")
expect_equal(length(pd_p5$pairs), 1)
expect_equal(length(pd_p5$pairs[[1]][1, ]), 2) #correct Format


#as_persistence.PHom correctly processes output from `ripserr::cubical` from {ripserr}
if (requireNamespace("ripserr", quietly = TRUE)) {
  pd <- ripserr::cubical(volcano)
  pd_p <- as_persistence(pd)

  rip_ver <- as.character(utils::packageVersion("ripserr"))

  if (utils::compareVersion(rip_ver, "0.3.0") >= 0) {
    #ripserr ≥ 0.3.0
    expect_inherits(pd_p, "persistence")
    expect_equal(pd_p$metadata$engine, "ripserr::<vietoris_rips/cubical>")
    expect_equal(pd_p$metadata$filtration, "Vietoris-Rips/cubical")
    expect_equal(length(pd_p$pairs[[1]][1, ]), 2)
  } else {
    #ripserr < 0.3.0
    expect_inherits(pd_p, "persistence")
    expect_equal(pd_p$metadata$engine, "ripserr::<vietoris_rips/cubical>")
    expect_equal(pd_p$metadata$filtration, "Vietoris-Rips/cubical")
    expect_equal(length(pd_p$pairs[[1]][1, ]), 2)
  }
}

#as_persistence.list handles a list of matrices correctly"
pd <- as_persistence(list(
  matrix(c(1, 2), ncol = 2),
  matrix(c(0.5, 2), ncol = 2)
))

expect_inherits(pd, "persistence")
expect_equal(length(pd$pairs), 2) #Two degrees (0 and 1)
expect_equal(pd$pairs[[1]], matrix(c(1, 2), ncol = 2))
expect_equal(pd$pairs[[2]], matrix(c(0.5, 2), ncol = 2))


#as_persistence.persistence returns the object unchanged
mat <- matrix(c(0, 1, 2, 1, 0.5, 2), ncol = 3, byrow = TRUE)
colnames(mat) <- c("dimension", "birth", "death")
pd <- as_persistence(mat)
pd2 <- as_persistence(pd)

expect_identical(pd, pd2)


#as.data.frame.persistence creates correct format
mat <- matrix(c(0, 1, 2, 1, 0.5, 2), ncol = 3, byrow = TRUE)
colnames(mat) <- c("dimension", "birth", "death")
pd <- as_persistence(mat)
df <- as.data.frame(pd)

expect_inherits(df, "data.frame")
expect_equal(colnames(df), c("dimension", "birth", "death"))
expect_equal(as.numeric(df[1, ]), c(0, 1, 2))


#`get_pairs` correctly grabs pairs from persistence object

pd_p <- as_persistence(pd)
expect_equal(get_pairs(pd_p, 1), pd_p$pairs[[2]])
expect_true(all(is.na(get_pairs(pd_p, 3)))) #nonexistent dimension returns empty matrix
expect_error(get_pairs(pd)) #verifies failure given incorrect class
