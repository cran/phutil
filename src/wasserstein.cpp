#include "hera/wasserstein.h"
#include "diagram_parser.h"

#include <cmath>
#include <limits>
#include <string>

double wassersteinDist(PairVector& diagramA,
                       PairVector& diagramB,
                       const double wasserstein_power = 1.0,
                       const double delta = 0.01,
                       const double internal_p = hera::get_infinity<double>(),
                       const double initial_epsilon = 0.0,
                       const double epsilon_common_ratio = 5.0,
                       const int max_bids_per_round = 1,
                       const int max_num_phases = std::numeric_limits<int>::max(),
                       const bool tolerate_max_iter_exceeded = false,
                       const bool return_matching = false,
                       const bool match_inf_points = true,
                       const bool print_relative_tolerance = false,
                       const bool verbose = false)
{
  using Real = double;
  hera::AuctionParams<Real> params;
  params.wasserstein_power = wasserstein_power;
  params.delta = delta;
  params.internal_p = internal_p;
  params.initial_epsilon = initial_epsilon;
  params.epsilon_common_ratio = epsilon_common_ratio;
  params.max_bids_per_round = max_bids_per_round;
  params.max_num_phases = max_num_phases;
  params.tolerate_max_iter_exceeded = tolerate_max_iter_exceeded;
  params.return_matching = return_matching;
  params.match_inf_points = match_inf_points;

  if (params.wasserstein_power < 1.0)
  {
    std::string msg = "Wasserstein_degree was \"" +
      std::to_string(params.wasserstein_power) +
      "\", must be a number >= 1.0. Cannot proceed.";
    cpp11::stop(msg.c_str());
  }

  if (params.wasserstein_power == 1.0)
  {
    hera::remove_duplicates<Real>(diagramA, diagramB);
  }

  if (params.delta <= 0.0)
  {
    std::string msg = "relative error was \"" +
      std::to_string(params.delta) +
      "\", must be a number > 0.0. Cannot proceed.";
    cpp11::stop(msg.c_str());
  }

  auto res = hera::wasserstein_cost_detailed(diagramA, diagramB, params);

  return res.distance;
}

[[cpp11::register]]
double wassersteinDistance(const cpp11::doubles_matrix<>& x,
                           const cpp11::doubles_matrix<>& y,
                           const double delta = 0.01,
                           const double wasserstein_power = 1.0)
{
  PairVector diagramA, diagramB;
  parseMatrix(x, diagramA);
  parseMatrix(y, diagramB);
  return wassersteinDist(diagramA, diagramB, wasserstein_power, delta);
}

[[cpp11::register]]
cpp11::doubles wassersteinPairwiseDistances(const cpp11::list& x,
                                            const double delta = 0.01,
                                            const double wasserstein_power = 1.0,
                                            const unsigned int ncores = 1)
{
  unsigned int N = x.size();
  unsigned int K = N * (N - 1) / 2;
  cpp11::writable::doubles result(K);
  std::vector<PairVector> pairs(N);

  for (int n = 0;n < N;++n)
  {
    auto diagram = cpp11::as_cpp<cpp11::doubles_matrix<>>(x[n]);
    parseMatrix(diagram, pairs[n]);
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores)
#endif
  for (int k = 0;k < K;++k)
  {
    unsigned int i = N - 2 - std::floor(std::sqrt(-8 * k + 4 * N * (N - 1) - 7) / 2.0 - 0.5);
    unsigned int j = k + i + 1 - N * (N - 1) / 2 + (N - i) * ((N - i) - 1) / 2;
    result[k] = wassersteinDist(pairs[i], pairs[j], wasserstein_power, delta);
  }

  return result;
}
