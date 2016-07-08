#pragma once
#include <complex>
#include <vector>
#include <functional>

namespace cheb {
  using CPLX = std::complex<double>;
  using CPoints = std::vector<CPLX>;
  using Points = std::vector<double>;
  const double PI = 3.141592653589793;

  const Points& chebpoints(size_t n);
  const Points& chebcoeffs(size_t degree, size_t density);

  using  Approximant = std::function<double(double)>;
  Points evaluate(const Approximant& f, size_t n);
  Points approximate(const Approximant& f);
  
  double evaluateT(const size_t degree, const double x);
}
