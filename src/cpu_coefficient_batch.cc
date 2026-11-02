#include "cpu_coefficient_batch.hh"
#include "reduction_table.hh"

#include <algorithm>
#include <chrono>
#include <mutex>
#include <stdexcept>

namespace {
unsigned int add_exponents(unsigned int a, unsigned int b, unsigned int zero,
                           const std::vector<std::int32_t> & increments) {
  if (a == zero) return b;
  if (b == zero) return a;
  if (a > b) std::swap(a, b);
  const unsigned int z = increments[b - a];
  if (z == zero) return zero;
  const auto s = a + z;
  return s >= zero ? s - zero : s;
}
using Clock = std::chrono::steady_clock;
std::uint64_t elapsed(Clock::time_point start) {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
}
}

std::shared_ptr<const CpuBatchField> CpuBatchField::get(const ReductionTable & t) {
  // Return no batch plan for unsupported fields or OpenCL counting.
  // Require 2*q < mod to recover signed character sums uniquely.
  if (t.is_opencl_enabled() || t.prime < 3 || t.prime > 1009 ||
      (t.prime_exponent != 1 && t.prime_exponent != 2) ||
      std::uint64_t(t.prime_power) * 2 >= HyCu::NttPlan::mod) return {};
  static std::mutex mutex;
  static std::vector<std::weak_ptr<const CpuBatchField>> plans;
  std::lock_guard<std::mutex> lock(mutex);
  for (auto it = plans.begin(); it != plans.end();) {
    auto plan = it->lock();
    if (!plan) { it = plans.erase(it); continue; }
    if (plan->p_ == t.prime && plan->exponent_ == t.prime_exponent &&
        plan->increments_ == *t.incrementation_table) return plan;
    ++it;
  }
  auto plan = std::shared_ptr<const CpuBatchField>(new CpuBatchField(t));
  plans.push_back(plan);
  return plan;
}

CpuBatchField::CpuBatchField(const ReductionTable & t)
  : p_(t.prime), q_(t.prime_power), exponent_(t.prime_exponent),
    ntt_(2 * p_ - 1), increments_(*t.incrementation_table), coordinates_(q_, q_) {
  const auto zero = q_ - 1;
  // Recover the prime subfield by addition in the existing Zech representation.
  std::vector<unsigned int> base(p_, zero);
  for (unsigned int a = 1; a < p_; ++a)
    base[a] = add_exponents(base[a - 1], 0, zero, increments_);
  const unsigned int rows = q_ / p_;
  kernels_.resize(rows, std::vector<std::uint32_t>(ntt_.size()));
  for (unsigned int b = 0; b < rows; ++b) {
    const unsigned int bg = b ? (base[b] + 1) % zero : zero;
    auto & kernel = kernels_[b];
    for (unsigned int a = 0; a < p_; ++a) {
      const auto e = add_exponents(base[a], bg, zero, increments_);
      if (coordinates_[e] != q_) throw std::runtime_error("Batch field basis collision");
      coordinates_[e] = b * p_ + a;
      kernel[a] = e == zero ? 0 : (e & 1) ? HyCu::NttPlan::mod - 1 : 1;
    }
#ifdef HYCU_CPU_BATCH_REFERENCE_NTT
    ntt_.transform(kernel);
#else
    ntt_.forward_spectrum(kernel);
#endif
  }
}

CpuCoefficientBatch::CpuCoefficientBatch(std::shared_ptr<ReductionTable> table,
                                       CpuBatchStatistics & statistics)
  : table_(std::move(table)), field_(CpuBatchField::get(*table_)), statistics_(statistics) {
  if (field_) {
    histogram_.resize(field_->q_);
    roots_.resize(field_->p_);
    unramified_.resize(field_->p_);
    work_.resize(field_->ntt_.size());
    spectrum_.resize(field_->ntt_.size());
  }
}

void CpuCoefficientBatch::begin_block(std::size_t varying) {
  if (field_ && (monomials_.empty() || varying_ != varying)) {
    monomials_.resize(field_->q_ - 1);
    for (unsigned int ix = 0; ix < monomials_.size(); ++ix)
      monomials_[ix] = std::uint64_t(varying) * (ix + 1) % monomials_.size();
  }
  varying_ = varying;
  signature_.clear();
  fallback_.reset(new CpuEvaluationCache(*table_, varying));
}

void CpuCoefficientBatch::prepare(const std::vector<unsigned int> & coefficients) {
  // Reuse family results only while all fixed coefficients stay unchanged.
  bool same = signature_.size() == coefficients.size();
  if (same) for (std::size_t i = 0; i < coefficients.size(); ++i)
    if (i != varying_ && signature_[i] != coefficients[i]) { same = false; break; }
  ++statistics_.lookups;
  if (same) return;
  signature_ = coefficients;
  // Zero the varying coefficient to evaluate H(x).
  signature_[varying_] = table_->prime_power_pred;
  ++statistics_.families;
  const auto start = Clock::now();
  const unsigned int zero = field_->q_ - 1, p = field_->p_;
  const auto & reductions = *table_->exponent_reduction_table;
  const auto & increments = *table_->incrementation_table;
  std::fill(histogram_.begin(), histogram_.end(), 0);
  std::fill(roots_.begin(), roots_.end(), 0);
  for (unsigned int x = 1; x <= zero; ++x) {
    unsigned int value = signature_[0], power = 0;
    for (std::size_t k = 1; k < signature_.size(); ++k) {
      power = reductions[power + x];
      if (signature_[k] != zero)
        value = add_exponents(value, reductions[signature_[k] + power], zero, increments);
    }
    const auto monomial = monomials_[x - 1];
    // t = H(x)/x^k; chi(f_c(x)) = chi(x^k) * chi(t+c).
    const auto quotient = value == zero ? zero :
        (value >= monomial ? value - monomial : value + zero - monomial);
    const auto coordinate = field_->coordinates_[quotient];
    histogram_[coordinate] += (monomial & 1) ? -1 : 1;
    // Count roots separately because signed histogram weights can cancel.
    if (coordinate < p) ++roots_[coordinate ? p - coordinate : 0];
  }
  statistics_.histogram_ns += elapsed(start);
  const auto convolution_start = Clock::now();
  std::fill(spectrum_.begin(), spectrum_.end(), 0);
  for (unsigned int b = 0; b < field_->q_ / p; ++b) {
    std::fill(work_.begin(), work_.end(), 0);
    bool nonzero = false;
    // Reverse indices to turn cyclic correlation into convolution.
    for (unsigned int a = 0; a < p; ++a) {
      const int h = histogram_[b*p + a];
      nonzero |= h != 0;
      work_[a ? p-a : 0] = h < 0 ? HyCu::NttPlan::mod - unsigned(-h) : unsigned(h);
    }
    if (!nonzero) continue;
#ifdef HYCU_CPU_BATCH_REFERENCE_NTT
    field_->ntt_.transform(work_);
#else
    field_->ntt_.forward_spectrum(work_);
#endif
    const auto & kernel = field_->kernels_[b];
    for (std::size_t i = 0; i < work_.size(); ++i)
      spectrum_[i] = HyCu::NttPlan::add(spectrum_[i], HyCu::NttPlan::mul(work_[i], kernel[i]));
  }
  // Linearity lets all row products share one inverse transform.
#ifdef HYCU_CPU_BATCH_REFERENCE_NTT
  field_->ntt_.transform(spectrum_, true);
#else
  field_->ntt_.inverse_spectrum(spectrum_);
#endif
  for (unsigned int c = 0; c < p; ++c) {
    // Fold indices modulo p and decode the signed character sum.
    const auto residue = HyCu::NttPlan::add(spectrum_[c], spectrum_[c+p]);
    const std::int64_t sum = residue > HyCu::NttPlan::mod/2
        ? std::int64_t(residue) - HyCu::NttPlan::mod : residue;
    // Subtract roots to keep only points with nonzero y.
    const std::int64_t count = std::int64_t(zero) - roots_[c] + sum;
    if (count < 0 || count > 2 * std::int64_t(zero) || (count & 1))
      throw std::runtime_error("Invalid batch character sum");
    unramified_[c] = count;
  }
  statistics_.convolution_ns += elapsed(convolution_start);
}
