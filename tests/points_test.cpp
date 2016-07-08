#include <gtest/gtest.h>
#include <cheb/Cheb.h>

TEST(points, basic) {
  for(size_t i=0; i < 100; i++) {
    const auto ptx = cheb::chebpoints(i);
    EXPECT_EQ(ptx.size(), i);
    for(const auto& pt : ptx) {
      EXPECT_GE(pt, -1.0);
      EXPECT_LE(pt,  1.0);
    }
  }
}

TEST(points, zero) {
  const auto ptx = cheb::chebpoints(1);
  EXPECT_EQ(ptx.size(), 1);
  EXPECT_NEAR(ptx.at(0), 0, 1e-10);
}

TEST(coeffs, quadratic) {
  const cheb::Points target = cheb::chebpoints(10);
  const cheb::Points coeffs = cheb::chebcoeffs(2, 10);
  EXPECT_EQ(target.size(), coeffs.size());
  for(size_t i=0; i < target.size(); i++) {
    const double x = i * 0.2 - 1;
    const double t2x = 2*x*x-1;
    EXPECT_NEAR(t2x, coeffs.at(i), 1e-10);
  }
}


TEST(approx, constant) {
  auto fx = [](double x){ return 1.0;};
  const cheb::Points out = cheb::approximate(fx);
  EXPECT_EQ(out.size(), 17);
  const cheb::Points coeffs = cheb::chebpoints(16);
  for(size_t i=0; i < out.size(); i++) {
    std::cout << i << " " << coeffs.at(i) << " " << out.at(i) << std::endl;
  }
}


TEST(points, print) {
  for(size_t i=0; i <= 20; i++) {
    const auto ptx = cheb::chebpoints(i);
    for(const auto& pt: ptx) {
      std::cout << pt << " ";
    }
    std::cout << std::endl;
  }
}


TEST(eval, t2) {
  for(double i=-1.0; i < 1.0; i += 0.1) {
    double observed = cheb::evaluateT(2,i);
    double expected = 2*i*i-1;
    EXPECT_NEAR(observed, expected, 1e-10);
  }
}


TEST(eval, t3) {
  for(double i=-1.0; i < 1.0; i += 0.1) {
    double observed = cheb::evaluateT(3,i);
    double expected = 4*i*i*i-3*i;
    EXPECT_NEAR(observed, expected, 1e-10);
  }
}


TEST(weights, t2t10) {
  const auto ptx = cheb::chebcoeffs(2,10);
  double sum{0.0};
  std::accumulate(ptx.begin(), ptx.end(), sum);
  EXPECT_NEAR(sum, 0, 1e-10);
  EXPECT_NEAR(ptx.at(0), .951, 1e-2);
}


TEST(eval, lowdeg) {
  auto t0 = [](double x){ return 1.0;};
  auto t1 = [](double x){ return x;};
  auto t2 = [](double x){ return 2*x*x-1;};
  auto fx = [](double x){ return x*x-2*x+1;};
  auto P0 = cheb::evaluate(t0, 1);
  EXPECT_EQ(P0.size(), 1);
  for(const auto ptx : P0) {
    EXPECT_NEAR(ptx, 1, 1e-10);
  }
  P0 = cheb::evaluate(t0, 2);
  EXPECT_EQ(P0.size(), 2);
  for(const auto ptx : P0) {
    EXPECT_NEAR(ptx, 1, 1e-10);
  }
  P0 = cheb::evaluate(t0, 30);
  EXPECT_EQ(P0.size(), 30);
  for(const auto ptx : P0) {
    EXPECT_NEAR(ptx, 1, 1e-10);
  }
  P0 = cheb::evaluate(t0, 50);
  EXPECT_EQ(P0.size(), 50);
  for(const auto ptx : P0) {
    EXPECT_NEAR(ptx, 1, 1e-10);
  }
  P0 = cheb::evaluate(t0, 200);
  EXPECT_EQ(P0.size(), 200);
  for(const auto ptx : P0) {
    EXPECT_NEAR(ptx, 1, 1e-10);
  }

  for(size_t density = 3; density < 10; density++) {
    const auto& p1 = cheb::evaluate(t1, density);
    EXPECT_EQ(p1.size(), density);
    double p1sum{0.0};
    std::accumulate(p1.begin(), p1.end(), p1sum);
    EXPECT_NEAR(p1sum, 0, 1e-10);
    for(size_t i=0; i <= density/2; i++) {
      EXPECT_NEAR(p1.at(i), -p1.at(density - 1 -i), 1e-10);
    }
    
    const auto& p2 = cheb::evaluate(t2, density);
    EXPECT_EQ(p2.size(), density);
    double p2sum{0.0};
    std::accumulate(p2.begin(), p2.end(), p2sum);
    EXPECT_NEAR(p2sum, 0, 1e-10);
    for(size_t i=0; i <= density/2; i++) {
      EXPECT_NEAR(p2.at(i), p2.at(density - 1 -i), 1e-10);
    }


  }
}
