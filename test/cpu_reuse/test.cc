#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "curve.hh"
#include "store/curve_data.hh"
#include "store/store.hh"
#include "store/store_data.hh"

namespace {

void require(bool condition, const std::string & message)
{
  // Unlike assert(), keep checks enabled in Release builds.
  if (!condition) throw std::runtime_error(message);
}

// Independent, small-field oracle: F_p[u]/(u^2-d), with d a nonsquare.
// No logarithm/reduction tables from the implementation under test are used.
std::tuple<unsigned int, unsigned int> direct_count(
    const Curve & curve, const FqElementTable & field, unsigned int exponent)
{
  const unsigned int p = curve.prime();
  unsigned int d = 2;
  std::vector<bool> prime_squares(p, false);
  for (unsigned int x = 0; x < p; ++x) prime_squares[x*x % p] = true;
  while (prime_squares[d]) ++d;
  const unsigned int q = exponent == 1 ? p : p*p;
  const auto multiply = [p, d](unsigned int x, unsigned int y) {
    const unsigned int a = x % p, b = x / p, c = y % p, e = y / p;
    return (a*c + d*b*e) % p + p*((a*e + b*c) % p);
  };
  std::vector<bool> squares(q, false);
  for (unsigned int y = 0; y < q; ++y) squares[multiply(y, y)] = true;
  std::vector<unsigned int> coefficients;
  for (auto e : curve.rhs_coeff_exponents()) coefficients.push_back(field.at_nmod(e));
  unsigned int unramified = 0, ramified = 0;
  for (unsigned int x = 0; x < q; ++x) {
    unsigned int value = 0;
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
      value = multiply(value, x);
      value = (value % p + *it) % p + p*(value / p);
    }
    if (value == 0) ++ramified;
    else if (squares[value]) unramified += 2;
  }
  if (curve.degree() % 2) ++ramified;
  else if (squares[coefficients.back()]) unramified += 2;
  return std::make_tuple(unramified, ramified);
}

void compare(Curve & baseline, Curve & reused, ReductionTable & fp,
             ReductionTable & fp2, CpuEvaluationCache & cache1,
             CpuEvaluationCache & cache2)
{
  baseline.count(fp);
  baseline.count(fp2);
  reused.count(cache1);
  reused.count(cache2);
  std::ostringstream context;
  context << baseline;
  require(baseline.number_of_points() == reused.number_of_points(),
          "Point-count mismatch: " + context.str());
  require(baseline.hasse_weil_offsets() == reused.hasse_weil_offsets(),
          "Hasse-Weil mismatch: " + context.str());
}

void test_coefficient_changes(unsigned int p)
{
  auto field = std::make_shared<FqElementTable>(p, 1);
  ReductionTable fp(p, 1), fp2(p, 2);
  for (unsigned int degree : {5u, 6u}) {
    for (unsigned int varying : {0u, 2u, degree}) {
      CpuEvaluationCache cache1(fp, varying), cache2(fp2, varying);
      std::vector<unsigned int> coefficients(degree + 1, p - 1);
      coefficients[0] = 0;
      coefficients[1] = 0;
      coefficients[degree - 1] = 0;
      coefficients[degree] = 0;
      // Exercise zero coefficients, leading coefficient removal/fallback,
      // cache hits, and changes to a supposedly fixed coefficient.
      for (unsigned int phase = 0; phase < 2; ++phase) {
        coefficients[varying == 0 ? 1 : 0] = phase;
        for (unsigned int a = 0; a < p; ++a) {
          coefficients[varying] = a;
          Curve baseline(field, coefficients), reused(field, coefficients);
          compare(baseline, reused, fp, fp2, cache1, cache2);
          if (p <= 7 && baseline.has_squarefree_rhs()) {
            for (unsigned int exponent : {1u, 2u})
              require(reused.number_of_points().at(exponent) ==
                      direct_count(reused, *field, exponent),
                      "Independent arithmetic mismatch at p=" + std::to_string(p)
                      + ", exponent=" + std::to_string(exponent));
          }
        }
      }
      require(cache1.rebuild_count() == 2 && cache2.rebuild_count() == 2,
              "Fixed coefficient change must rebuild each field cache");
      require(cache1.reuse_count() > 0 && cache2.reuse_count() > 0,
              "Expected actual reuse, not only equivalent recomputation");
      // Same cache, changed polynomial degree: signature must invalidate it.
      coefficients[varying] = 0;
      coefficients.push_back(0);
      Curve baseline(field, coefficients), reused(field, coefficients);
      compare(baseline, reused, fp, fp2, cache1, cache2);
      require(cache1.rebuild_count() == 3 && cache2.rebuild_count() == 3,
              "Degree change must rebuild each field cache");
    }
  }
}

std::string test_enumeration(unsigned int p, unsigned int package_size,
                             bool marked, std::size_t & total_reuses)
{
  auto field = std::make_shared<FqElementTable>(p, 1);
  ReductionTable fp(p, 1), fp2(p, 2);
  using ResultStore = Store<HyCu::CurveData::ExplicitRamificationHasseWeil,
                            HyCu::StoreData::Count>;
  ResultStore baseline_store, reused_store;
  std::size_t tested = 0;
  for (CurveIterator blocks(*field, 2, marked, package_size);
       !blocks.is_end(); blocks.step()) {
    const auto bounds = blocks.as_block();
    const auto varying = CpuEvaluationCache::select_coefficient(bounds);
    CpuEvaluationCache cache1(fp, varying), cache2(fp2, varying);
    for (BlockIterator it(bounds); !it.is_end(); it.step()) {
      Curve baseline(field, it.as_position()), reused(field, it.as_position());
      if (!baseline.has_squarefree_rhs()) continue;
      compare(baseline, reused, fp, fp2, cache1, cache2);
      require(baseline.ramification_type() == reused.ramification_type(),
              "Ramification mismatch");
      baseline_store.register_curve(baseline);
      reused_store.register_curve(reused);
      ++tested;
    }
    total_reuses += cache1.reuse_count() + cache2.reuse_count();
  }
  std::ostringstream baseline_output, reused_output;
  baseline_store.insert(baseline_output);
  reused_store.insert(reused_output);
  require(tested > 0, "Enumeration must exercise squarefree curves");
  require(baseline_output.str() == reused_output.str(),
          "Aggregated output mismatch at p=" + std::to_string(p));
  return reused_output.str();
}

} // namespace

int main()
{
  try {
    require(CpuEvaluationCache::select_coefficient({{0,1}, {0,3}, {1,4}}) == 2,
            "Coefficient selection must break ties by highest index");
    require(CpuEvaluationCache::select_coefficient({{0,5}, {0,2}}) == 0,
            "Coefficient selection must prefer widest range");
    require(CpuEvaluationCache::select_coefficient({{0,1}, {2,3}}) == 2,
            "Singleton blocks must select fallback");
    for (unsigned int p : {3u,5u,7u,11u,13u,17u,19u,23u,29u,31u,37u,41u,43u,47u,53u,59u}) {
      test_coefficient_changes(p);
      std::cout << "Coefficient tests passed: p=" << p << std::endl;
    }
    std::size_t total_reuses = 0;
    for (unsigned int p : {3u,5u,7u}) {
      for (bool marked : {false, true}) {
        const auto reference = test_enumeration(p, 1, marked, total_reuses);
        for (unsigned int package : {7u,30u,2500u})
          require(reference == test_enumeration(p, package, marked, total_reuses),
                  "Result changed with package size at p=" + std::to_string(p));
      }
      std::cout << "Full enumeration passed: p=" << p << std::endl;
    }
    require(total_reuses > 0, "Enumeration must exercise actual cache reuse");
    std::cout << "All CPU reuse checks passed. Reuses in enumeration: "
              << total_reuses << std::endl;
    return 0;
  }
  catch (const std::exception & error) {
    std::cerr << "FAIL: " << error.what() << std::endl;
    return 1;
  }
}
