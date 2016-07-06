#pragma once
#include <complex>
#include <vector>
#include <functional>

namespace cheb {
  using CPLX = std::complex<double>;
  using CPoints = std::vector<CPLX>;
  using Points = std::vector<double>;
  Points  chebpoints(size_t n);
  const double PI = 3.141592653589793;
  using  Approximant = std::function<double(double)>;
  Points evaluate(const Approximant& f, size_t n);
  const Points& chebcoeffs(size_t degree, size_t density);
}
