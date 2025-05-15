#include "hera/bottleneck.h"
#include "diagram_parser.h"

#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

double bottleneckDist(const PairVector& diagramA,
                      const PairVector& diagramB,
                      const double delta = 0.01)
{
  hera::bt::MatchingEdge<double> e;

  if (delta > 0.0)
  {
    return hera::bottleneckDistApprox(diagramA, diagramB, delta, e, true);
  }

  if (delta == 0.0)
  {
    int decPrecision { 0 };
    return hera::bottleneckDistExact(diagramA, diagramB, decPrecision);
  }

  std::string msg = "relative error was \"" +
    std::to_string(delta) +
    "\", must be a number >= 0.0. Cannot proceed.";
  cpp11::stop(msg.c_str());
}

[[cpp11::register]]
double bottleneckDistance(const cpp11::doubles_matrix<>& x,
                          const cpp11::doubles_matrix<>& y,
                          const double delta = 0.01)
{
  PairVector diagramA, diagramB;
  parseMatrix(x, diagramA);
  parseMatrix(y, diagramB);
  return bottleneckDist(diagramA, diagramB, delta);
}

[[cpp11::register]]
cpp11::doubles bottleneckPairwiseDistances(const cpp11::list& x,
                                           const double delta = 0.01,
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
    result[k] = bottleneckDist(pairs[i], pairs[j], delta);
  }

  return result;
}
