#include "diagram_parser.h"

void parseMatrix(const cpp11::doubles_matrix<>& matrix, PairVector& result) {
  unsigned int numPairs = matrix.nrow();
  result.reserve(numPairs);
  for (unsigned int i = 0;i < numPairs; ++i) {
    result.emplace_back(matrix(i, 0), matrix(i, 1));
  }
}
