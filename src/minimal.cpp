#include <cheb/Cheb.h>

#include <iostream>
#include <map>
#include <utility>
namespace cheb {
  // TODO(JCK): cache these, return reference to cache
  Points chebpoints(size_t n) {
    if (n==0) {
      return Points{1};
    }
    Points out;
    out.reserve(n+1);
    for(size_t i=0; i <= n; i++) {
      CPLX arg(0,PI * i / (1.0 * n));
      CPLX z = std::exp(arg);
      out.push_back(z.real());
    }
    return out;
  }

  static std::map<std::pair<size_t, size_t>, Points> coeffCache;
  const Points& chebcoeffs(size_t degree, size_t density) {
    const auto key = std::make_pair(degree,density);
    if (coeffCache.count(key)){
      return coeffCache.at(key);
    }
    Points out;
    out.reserve(density+1);
    for(size_t i=0; i <= density; i++) {
      // x : i/(1.0 * density)
      const double x = 1 - 2.0*(density - i) /( density);
      const double t = degree * std::acos(x);
      out.push_back(std::cos(t));
    }
    coeffCache[key] = out;
    
    return coeffCache.at(key);
    
  }
  Points evaluate(const Approximant& f, size_t n) {
    const Points& ptx = chebpoints(n);
    Points out;
    out.reserve(ptx.size());
    for(const auto& pt : ptx) {
      out.push_back(f(pt));
    }
    return out;
  }
}
