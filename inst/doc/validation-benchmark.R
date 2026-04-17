## -----------------------------------------------------------------------------
#| label: setup
#| include: false
library(phutil)
init_par <- par()


## ----define small PDs---------------------------------------------------------
X <- rbind(
  c(1, 3),
  c(3, 5)
)
Y <- rbind(
  c(3, 4)
)


## -----------------------------------------------------------------------------
#| label: fig-plot-small
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
#| fig-cap: "Overlaid persistence diagrams $X$ (circles) and $Y$ (diamond) with dashed segments connecting optimally matched pairs."
oldpar <- par(mar = c(4, 4, 1, 1) + .1)
plot(
  NA_real_,
  xlim = c(0, 6), ylim = c(0, 6), asp = 1, xlab = "birth", ylab = "death"
)
abline(a = 0, b = 1)
points(X, pch = 1)
points(Y, pch = 5)
segments(X[, 1], X[, 2], c(2, Y[, 1]), c(2, Y[, 2]), lty = 2)
par(mar = init_par$mar)


## ----validate small PDs with Hera---------------------------------------------
wasserstein_distance(X, Y, p = 1)
wasserstein_distance(X, Y, p = 2)
bottleneck_distance(X, Y)


## ----validate small PDs with Dionysus-----------------------------------------
TDA::wasserstein(cbind(0, X), cbind(0, Y), p = 1, dimension = 0)
sqrt(TDA::wasserstein(cbind(0, X), cbind(0, Y), p = 2, dimension = 0))
TDA::bottleneck(cbind(0, X), cbind(0, Y), dimension = 0)


## ----validate small PD vs empty-----------------------------------------------
# empty PD
E <- matrix(NA_real_, nrow = 0, ncol = 2)
# with dimension column
E_ <- cbind(matrix(NA_real_, nrow = 0, ncol = 1), E)
# distance from empty using phutil/Hera
wasserstein_distance(E, X, p = 1)
wasserstein_distance(E, X, p = 2)
bottleneck_distance(E, X)
# distance from empty using TDA/Dionysus
TDA::wasserstein(E_, cbind(0, X), p = 1, dimension = 0)
sqrt(TDA::wasserstein(E_, cbind(0, X), p = 2, dimension = 0))
TDA::bottleneck(E_, cbind(0, X), dimension = 0)


## ----compute large PDs fake, eval=FALSE---------------------------------------
# set.seed(28415)
# n <- 24
# PDs1 <- lapply(seq(n), function(i) {
#   S1 <- tdaunif::sample_trefoil(n = 120, sd = .05)
#   as_persistence(TDA::ripsDiag(S1, maxdimension = 2, maxscale = 6))
# })
# PDs2 <- lapply(seq(n), function(i) {
#   S2 <- cbind(tdaunif::sample_arch_spiral(n = 120, arms = 2), 0)
#   S2 <- tdaunif::add_noise(S2, sd = .05)
#   as_persistence(TDA::ripsDiag(S2, maxdimension = 2, maxscale = 6))
# })


## ----compute large PDs true, echo=FALSE---------------------------------------
n <- 24
PDs1 <- trefoils
PDs2 <- arch_spirals


## -----------------------------------------------------------------------------
#| label: benchmark phutil and TDA
#| warning: false
PDs1_ <- lapply(lapply(PDs1, as.data.frame), as.matrix)
PDs2_ <- lapply(lapply(PDs2, as.data.frame), as.matrix)
# iterate over homological degrees and Wasserstein powers
bm_all <- list()
PDs_i <- seq_along(PDs1)
for (dimension in seq(0, 2)) {
  # compute
  bm_1 <- do.call(rbind, lapply(seq_along(PDs1), function(i) {
    as.data.frame(microbenchmark::microbenchmark(
      TDA = TDA::wasserstein(
        PDs1_[[i]], PDs2_[[i]], dimension = dimension, p = 1
      ),
      phutil = wasserstein_distance(
        PDs1[[i]],  PDs2[[i]],  dimension = dimension, p = 1
      ),
      times = 1, unit = "ns"
    ))
  }))
  bm_2 <- do.call(rbind, lapply(seq_along(PDs1), function(i) {
    as.data.frame(microbenchmark::microbenchmark(
      TDA = sqrt(TDA::wasserstein(
        PDs1_[[i]], PDs2_[[i]], dimension = dimension, p = 2
      )),
      phutil = wasserstein_distance(
        PDs1[[i]],  PDs2[[i]],  dimension = dimension, p = 2
      ),
      times = 1, unit = "ns"
    ))
  }))
  bm_inf <- do.call(rbind, lapply(seq_along(PDs1), function(i) {
    as.data.frame(microbenchmark::microbenchmark(
      TDA = TDA::bottleneck(
        PDs1_[[i]], PDs2_[[i]], dimension = dimension
      ),
      phutil = bottleneck_distance(
        PDs1[[i]],  PDs2[[i]],  dimension = dimension
      ),
      times = 1, unit = "ns"
    ))
  }))
  # annotate and combine
  bm_1$power <- 1; bm_2$power <- 2; bm_inf$power <- Inf
  bm_res <- rbind(bm_1, bm_2, bm_inf)
  bm_res$degree <- dimension
  bm_all <- c(bm_all, list(bm_res))
}
bm_all <- do.call(rbind, bm_all)


## -----------------------------------------------------------------------------
#| label: fig-benchmark-large
#| fig-width: 8
#| fig-height: 4
#| fig-align: 'center'
#| fig-retina: 2
#| fig-cap: "Benchmark comparison of Dionysus via {TDA} and Hera via {phutil} on
#| large persistence diagrams: Jitter plots of runtime distributions
#| (time measured in seconds)."
bm_all <- transform(
  bm_all,
  expr = factor(as.character(expr), levels = c("TDA", "phutil")),
  time = unlist(time) * 10e-9
)
bm_all <- subset(bm_all, select = c(expr, degree, power, time))
xrans <- lapply(seq(0, 2), function(d) range(subset(bm_all, degree == d, time)))
par(mfcol = c(3, 3), mar = c(2, 2, 2, 2) + .1)
for (d in seq(0, 2)) for (p in c(1, 2, Inf)) {
  bm_d_p <- subset(bm_all, degree == d & power == p)
  plot(
    x = bm_d_p$time, xlim = xrans[[d + 1]],
    y = jitter(as.integer(bm_d_p$expr)), yaxt = "n",
    pch = 19
  )
  axis(2, at = c(1, 2), labels = levels(bm_d_p$expr))
  if (p == 1) axis(
    3, at = mean(xrans[[d+1]]),
    tick = FALSE, labels = paste("degree: ", d), padj = 0
  )
  if (d == 2) axis(
    4, at = 1.5,
    tick = FALSE, labels = paste("power: ", p), padj = 0
  )
}
par(mfcol = init_par$mfcol)

