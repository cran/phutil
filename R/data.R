#' Toy Data: Noisy circle
#'
#' A simulated data set consisting of 100 points sampled from a circle with
#' additive Gaussian noise using a standard deviation of 0.05.
#'
#' The point cloud stored in `noisy_circle_points` has been generated using the
#' [**tdaunif**](https://tdaverse.github.io/tdaunif/) package using the
#' [`tdaunif::sample_circle()`](https://tdaverse.github.io/tdaunif/reference/circles.html)
#' function. Specifically, the following parameters were used: `n = 100`, `sd =
#' 0.05` and a seed of 1234.
#'
#' The persistence diagram stored in `noisy_circle_ripserr` has been computed
#' using the [**ripserr**](https://tdaverse.github.io/ripserr/) package with the
#' [`ripserr::vietoris_rips()`](https://tdaverse.github.io/ripserr/reference/vietoris_rips.html)
#' function. Specifically, the following parameters were used: `max_dim = 1L`.
#'
#' The persistence diagram stored in `noisy_circle_tda_rips` has been computed
#' using the **TDA** package with the
#' [`TDA::ripsDiag()`](https://www.rdocumentation.org/packages/TDA/versions/1.9.1/topics/ripsDiag)
#' function. Specifically, the following parameters were used: `maxdimension =
#' 1L` and `maxscale = 1.6322`.
#'
#' @format
#' ## `noisy_circle_points`
#' A matrix with 100 rows and 2 columns listing the coordinates of the points.
#'
#' ## `noisy_circle_ripserr`
#' An object of class 'PHom' as returned by the
#' [`ripserr::vietoris_rips()`](https://tdaverse.github.io/ripserr/reference/vietoris_rips.html)
#' function, which is a data frame with 3 variables:
#'
#' - `dimension`: the dimension/degree of the feature,
#' - `birth`: the birth value of the feature,
#' - `death`: the death value of the feature.
#'
#' ## `noisy_circle_tda_rips`
#' A list of length 1 containing an object of class 'diagram' as returned by the
#' [`TDA::ripsDiag()$diagram`](https://www.rdocumentation.org/packages/TDA/versions/1.9.1/topics/ripsDiag)
#' function, which is a matrix with 3 columns:
#'
#' - `dimension`: the dimension/degree of the feature,
#' - `birth`: the birth value of the feature,
#' - `death`: the death value of the feature.
#'
#' @source <https://tdaverse.github.io/tdaunif/reference/circles.html>,
#'   <https://tdaverse.github.io/ripserr/reference/vietoris_rips.html>,
#'   <https://www.rdocumentation.org/packages/TDA/versions/1.9.1/topics/ripsDiag>
#' @name noisy_circle

#' @rdname noisy_circle
"noisy_circle_points"

#' @rdname noisy_circle
"noisy_circle_ripserr"

#' @rdname noisy_circle
"noisy_circle_tda_rips"

#' Toy Data: A sample of persistence diagrams
#'
#' A collection of 100 samples of size 100 on the sphere from which a
#' persistence diagram is computed using the [`TDA::ripsDiag()`] function with
#' parameters `maxdimension = 1L` and `maxscale = 1.6322`. Each diagram has been
#' generated using the [`tdaunif::sample_2sphere()`] function with the following
#' parameters: `n = 100` and `sd = 0.05`. The seed was fixed to 1234.
#'
#' @format
#' A list of length 100, where each element is an object of class 'persistence'.
"persistence_sample"

#' A sample of persistence diagrams from the trefoil
#'
#' A collection of 24 samples of size 120 on the trefoil from which a
#' persistence diagram is computed using the [`TDA::ripsDiag()`] function with
#' `maxdimension = 2` and `maxscale = 6`. Each diagram has been generated using the
#' [`tdaunif::sample_trefoil()`] function with the following parameters: `n =
#' 120` and `sd = 0.05`. The seed was fixed to 28415.
#'
#' @format
#' A list of length 24, where each element is an object of class 'persistence'.
"trefoils"

#' A sample of persistence diagrams from the arch spiral
#'
#' A collection of 24 samples of size 120 on the arch spiral from which a
#' persistence diagram is computed using the [`TDA::ripsDiag()`] function with
#' `maxdimension = 2` and `maxscale = 6`. Each diagram has been generated using
#' the [`tdaunif::sample_arch_spiral()`] function with the following parameters:
#' `n = 120`, `arms = 2`` and `sd = 0.05`. The seed was fixed to 28415.
#'
#' @format
#' A list of length 24, where each element is an object of class 'persistence'.
"arch_spirals"
