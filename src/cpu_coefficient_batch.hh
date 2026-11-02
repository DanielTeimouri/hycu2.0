#ifndef HYCU_CPU_COEFFICIENT_BATCH_HH
#define HYCU_CPU_COEFFICIENT_BATCH_HH

#include <cstdint>
#include <memory>
#include <vector>
#include "cpu_evaluation_cache.hh"
#include "ntt.hh"

class ReductionTable;
class Curve;

struct CpuBatchStatistics {
  std::uint64_t families = 0, lookups = 0, fallbacks = 0;
  std::uint64_t histogram_ns = 0, convolution_ns = 0;
};

// Coordinates a+b*g use the generator from ReductionTable.
// Share field plans only when the Zech tables match.
class CpuBatchField {
public:
  static std::shared_ptr<const CpuBatchField> get(const ReductionTable & table);
private:
  friend class CpuCoefficientBatch;
  friend class Curve;
  explicit CpuBatchField(const ReductionTable & table);
  unsigned int p_, q_, exponent_;
  HyCu::NttPlan ntt_;
  std::vector<std::int32_t> increments_;
  std::vector<unsigned int> coordinates_;
  std::vector<std::vector<std::uint32_t>> kernels_;
};

// Per-worker workspace; family results are reset at each block boundary.
class CpuCoefficientBatch {
public:
  CpuCoefficientBatch(std::shared_ptr<ReductionTable> table,
                      CpuBatchStatistics & statistics);
  void begin_block(std::size_t varying);
private:
  friend class Curve;
  void prepare(const std::vector<unsigned int> & coefficients);
  std::shared_ptr<ReductionTable> table_;
  std::shared_ptr<const CpuBatchField> field_;
  CpuBatchStatistics & statistics_;
  std::size_t varying_ = 0;
  std::unique_ptr<CpuEvaluationCache> fallback_;
  std::vector<unsigned int> signature_, monomials_, roots_, unramified_;
  std::vector<int> histogram_;
  std::vector<std::uint32_t> work_, spectrum_;
};

#endif
