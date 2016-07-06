#include <gtest/gtest.h>
#include <cheb/Cheb.h>

TEST(points, basic) {
  for(size_t i=0; i < 100; i++) {
    const auto ptx = cheb::chebpoints(i);
    EXPECT_EQ(ptx.size(), i+1);
    for(const auto& pt : ptx) {
      EXPECT_GE(pt, -1.0);
      EXPECT_LE(pt,  1.0);
    }
  }
}

TEST(points, zero) {
  const auto ptx = cheb::chebpoints(0);
  EXPECT_EQ(ptx.size(), 1);
  EXPECT_EQ(ptx.at(0), 1);
}

TEST(coeffs, quadratic) {
  const cheb::Points target = cheb::chebpoints(10);
  const cheb::Points coeffs = cheb::chebcoeffs(2, 10);
  EXPECT_EQ(target.size(), coeffs.size());
  for(size_t i=0; i < target.size(); i++) {
    const double x = i * 0.2 - 1;
    const double t2x = 2*x*x-1;
    EXPECT_NEAR(t2x, coeffs.at(i), 1e-10);
    std::cout << i << " " << x << " " << t2x << " " << coeffs.at(i) << std::endl;
  }
}
