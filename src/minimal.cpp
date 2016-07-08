#include <cheb/Cheb.h>

#include <iostream>
#include <map>
#include <utility>
#include <cassert>
namespace cheb {
  // TODO(JCK): cache these, return reference to cache
  static std::map<size_t, Points> pointsCache;
  
  const Points& chebpoints(size_t n) {
    if (pointsCache.count(n)) {
      return pointsCache.at(n);
    }
    if (n==0) {
      pointsCache[0] = Points{};
      return pointsCache.at(0);
    }
    Points out;
    out.reserve(n);
    // for(size_t i=0; i <= n; i++) {
    //   CPLX arg(0,PI * i / (1.0 * n));
    //   CPLX z = std::exp(arg);
    //   out.push_back(z.real());
    // }
    for(size_t i=0; i < n; i++) {
      
      out.push_back(std::cos((2*i+1.0) / (2.0 * n)*cheb::PI));
    }
    pointsCache[n] = out;
    return pointsCache.at(n);
  }
  
  double evaluateT(const size_t degree, const double x) {
    if (degree == 0){ return 1;}
    if (degree == 1){ return x;}
    double a = 1;
    double b = x;
    double c;
    for(size_t i=2; i <= degree; i++) {
      c = 2*x*b - a;
      a = b;
      b = c;
    }
    return c;
  }

  static std::map<std::pair<size_t, size_t>, Points> coeffCache;

  // This method evaluates the specified chebpoly at the specified
  // chebpoints
  const Points& chebcoeffs(size_t degree, size_t density) {
    const auto key = std::make_pair(degree,density);
    if (coeffCache.count(key)){
      return coeffCache.at(key);
    }
    Points out;
    const Points& nodes = cheb::chebpoints(density);
    out.reserve(nodes.size());
    for(const auto& pt: nodes) {
      out.push_back(cheb::evaluateT(degree, pt));
    }
    coeffCache[key] = out;
    return coeffCache.at(key); 
  }
  
  Points evaluate(const Approximant& f, size_t degree) {
    const Points& ptx = chebpoints(degree);
    Points out;
    out.reserve(ptx.size());
    for(const auto& pt : ptx) {
      out.push_back(f(pt));
    }
    return out;
  }

  Points approximate(const Approximant& f) {
    size_t start{10};
    Points fx = evaluate(f, start);
    Points out;
    
    for(size_t i=0; i <= start; i++) {
      // take dot product of fx and chebcoeffs
      const Points& localPoints = chebcoeffs(i, start);
      double sum{0.0};
      for(size_t j=0; j < fx.size(); j++) {
	std::cout << i << " " << j << " "
		  << fx.at(j) << " " << localPoints.at(j) << " " << sum << std::endl;
	sum += fx.at(j) * localPoints.at(j);
      }
      std::cout << "SUM : " << sum << std::endl;
      out.push_back(sum);
    }
    return out;
  }
}
