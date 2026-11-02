#ifndef HYCU_CPU_EVALUATION_CACHE_HH
#define HYCU_CPU_EVALUATION_CACHE_HH

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

class Curve;
class ReductionTable;

// One cache per worker, block and counting field. The table must outlive it
// and must not be rebuilt while the cache is in use.
class CpuEvaluationCache
{
  public:
    CpuEvaluationCache(ReductionTable & table, std::size_t varying_coefficient)
      : table_(table), varying_coefficient_(varying_coefficient) {}

    // Match BlockIterator's fastest-changing dimension: widest range, then
    // highest coefficient index. Return bounds.size() for a singleton block.
    static std::size_t select_coefficient(
        const std::vector<std::tuple<unsigned int, unsigned int>> & bounds)
    {
      std::size_t selected = bounds.size();
      unsigned int widest = 1;
      for (std::size_t ix = 0; ix < bounds.size(); ++ix) {
        const unsigned int lower = std::get<0>(bounds[ix]);
        const unsigned int upper = std::get<1>(bounds[ix]);
        if (upper > lower && upper - lower > 1 && upper - lower >= widest) {
          widest = upper - lower;
          selected = ix;
        }
      }
      return selected;
    }

    std::size_t rebuild_count() const { return rebuild_count_; }
    std::size_t reuse_count() const { return reuse_count_; }

  private:
    friend class Curve;
    ReductionTable & table_;
    const std::size_t varying_coefficient_;
    std::vector<unsigned int> fixed_coefficients_;
    std::vector<std::uint32_t> fixed_values_;
    std::vector<std::uint32_t> monomial_exponents_;
    std::size_t rebuild_count_ = 0;
    std::size_t reuse_count_ = 0;
};

#endif
