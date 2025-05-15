#include <cpp11.hpp>
#include <utility>
#include <vector>

using PairVector = std::vector<std::pair<double,double>>;

void parseMatrix(const cpp11::doubles_matrix<>& matrix, PairVector& result);
