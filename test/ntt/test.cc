#include "ntt.hh"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

namespace {
using Plan = HyCu::NttPlan;
using Values = std::vector<std::uint32_t>;
void require(bool ok, const char * message) {
  if (!ok) throw std::runtime_error(message);
}
std::size_t reversed(std::size_t value, std::size_t n) {
  std::size_t result = 0;
  for (; n > 1; n >>= 1, value >>= 1) result = 2*result + (value&1);
  return result;
}
void transforms() {
  std::mt19937 random(20260918);
  for (std::size_t n = 1; n <= 4096; n *= 2) {
    Plan plan(n);
    for (unsigned int trial = 0; trial < 36; ++trial) {
      Values input(n);
      for (std::size_t i = 0; i < n; ++i) {
        if (trial == 0) input[i] = 0;
        else if (trial == 1) input[i] = Plan::mod-1;
        else if (trial == 2) input[i] = (i == n-1 ? 1 : 0);
        else if (trial == 3) input[i] = i&1 ? Plan::mod-1 : 1;
        else input[i] = random()%Plan::mod;
      }
      auto reference = input, fast = input;
      plan.transform(reference);
      plan.forward_spectrum(fast);
      for (std::size_t i = 0; i < n; ++i)
        require(fast[reversed(i,n)] == reference[i], "Spectrum differs from reference");
      if (n <= 32) {
        // Direct O(n^2) DFT, independent of butterfly traversal/permutation.
        const auto root = Plan::power(3,(Plan::mod-1)/n);
        for (std::size_t k = 0; k < n; ++k) {
          std::uint32_t expected = 0;
          for (std::size_t j = 0; j < n; ++j)
            expected = Plan::add(expected,Plan::mul(input[j],Plan::power(root,(j*k)%n)));
          require(fast[reversed(k,n)] == expected,"Direct DFT mismatch");
        }
      }
      plan.inverse_spectrum(fast);
      require(fast == input,"Roundtrip mismatch");
      // Arbitrary frequency input, not just outputs from our forward path.
      reference = input;
      for (std::size_t i = 0; i < n; ++i) fast[reversed(i,n)] = input[i];
      plan.transform(reference,true);
      plan.inverse_spectrum(fast);
      require(fast == reference,"Inverse spectrum mismatch");
    }
  }
  bool rejected = false;
  try { Plan invalid(0); } catch (const std::invalid_argument &) { rejected = true; }
  require(rejected,"Zero transform length accepted");
  Plan plan(8);
  Values short_buffer(7);
  rejected = false;
  try { plan.forward_spectrum(short_buffer); } catch (const std::invalid_argument &) { rejected = true; }
  require(rejected,"Wrong forward buffer length accepted");
  rejected = false;
  try { plan.inverse_spectrum(short_buffer); } catch (const std::invalid_argument &) { rejected = true; }
  require(rejected,"Wrong inverse buffer length accepted");
}
void correlations() {
  for (unsigned int p : {3u,43u,47u,59u,101u,127u,131u,163u,199u,257u,509u,1009u}) {
    Plan plan(2*p-1);
    Values a(plan.size()), b(plan.size());
    std::vector<int> h(p), k(p);
    for (unsigned int i = 0; i < p; ++i) {
      h[i] = int((i*13+7)%19)-9;
      k[i] = int((i*11+3)%3)-1;
      a[i ? p-i : 0] = h[i]<0 ? Plan::mod-unsigned(-h[i]) : unsigned(h[i]);
      b[i] = k[i]<0 ? Plan::mod-1 : unsigned(k[i]);
    }
    plan.forward_spectrum(a); plan.forward_spectrum(b);
    for (std::size_t i = 0; i < a.size(); ++i) a[i] = Plan::mul(a[i],b[i]);
    plan.inverse_spectrum(a);
    for (unsigned int c = 0; c < p; ++c) {
      std::int64_t expected = 0;
      for (unsigned int i = 0; i < p; ++i) expected += h[i]*k[(i+c)%p];
      const auto residue = Plan::add(a[c],a[c+p]);
      const auto actual = residue>Plan::mod/2 ? std::int64_t(residue)-Plan::mod : residue;
      require(actual == expected,"Direct signed correlation mismatch");
    }
  }
}

struct Family {
  Plan plan;
  std::vector<Values> inputs, old_kernels, fast_kernels;
  Values work, spectrum;
  explicit Family(unsigned int p) : plan(2*p-1), work(plan.size()), spectrum(plan.size()) {
    std::mt19937 random(p);
    for (unsigned int row = 0; row < p; ++row) {
      Values input(plan.size()), kernel(plan.size());
      for (unsigned int i = 0; i < p; ++i) {
        const int h = int(random()%9)-4, k = int(random()%3)-1;
        input[i ? p-i : 0] = h<0 ? Plan::mod-unsigned(-h) : unsigned(h);
        kernel[i] = k<0 ? Plan::mod-1 : unsigned(k);
      }
      inputs.push_back(input);
      auto fast_kernel = kernel;
      plan.transform(kernel); plan.forward_spectrum(fast_kernel);
      old_kernels.push_back(std::move(kernel));
      fast_kernels.push_back(std::move(fast_kernel));
    }
  }
  template<bool Fast> void calculate() {
    std::fill(spectrum.begin(),spectrum.end(),0);
    const auto & kernels = Fast ? fast_kernels : old_kernels;
    for (std::size_t row = 0; row < inputs.size(); ++row) {
      work = inputs[row];
      if (Fast) plan.forward_spectrum(work); else plan.transform(work);
      for (std::size_t i = 0; i < work.size(); ++i)
        spectrum[i] = Plan::add(spectrum[i],Plan::mul(work[i],kernels[row][i]));
    }
    if (Fast) plan.inverse_spectrum(spectrum); else plan.transform(spectrum,true);
  }
  template<bool Fast> double measure(unsigned int iterations) {
    const auto start = std::chrono::steady_clock::now();
    for (unsigned int i = 0; i < iterations; ++i) calculate<Fast>();
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/iterations;
  }
};
void benchmark() {
  for (unsigned int p : {43u,59u,101u,127u,131u,163u,199u,509u,1009u}) {
    Family family(p);
    family.calculate<false>();
    const auto expected = family.spectrum;
    family.calculate<true>();
    require(expected == family.spectrum,"Family convolution differs");
    const unsigned int iterations = std::max(3u, unsigned(1000000/(p*family.plan.size())));
    // Warm both paths; alternate sample order to reduce ordering bias.
    family.measure<false>(iterations); family.measure<true>(iterations);
    for (unsigned int sample = 0; sample < 7; ++sample) {
      double old_time, fast_time;
      if (sample&1) { fast_time=family.measure<true>(iterations); old_time=family.measure<false>(iterations); }
      else { old_time=family.measure<false>(iterations); fast_time=family.measure<true>(iterations); }
      require(expected == family.spectrum,"Benchmark result differs");
      std::cout << std::setprecision(10) << "{\"prime\":" << p
        << ",\"ntt_length\":" << family.plan.size() << ",\"sample\":" << sample+1
        << ",\"iterations\":" << iterations << ",\"reference_seconds_per_family\":" << old_time
        << ",\"optimized_seconds_per_family\":" << fast_time << "}" << std::endl;
    }
  }
}
}
int main(int argc, char ** argv) {
  try {
    transforms(); correlations();
    std::cerr << "PASS: reference/direct DFT, inverse, roundtrips, buffer guards, signed correlation through p=1009\n";
    if (argc == 2 && std::string(argv[1]) == "--benchmark") benchmark();
    else if (argc != 1) throw std::invalid_argument("usage: hycu-ntt-test [--benchmark]");
    return 0;
  } catch (const std::exception & e) { std::cerr << "FAIL: " << e.what() << std::endl; return 1; }
}
