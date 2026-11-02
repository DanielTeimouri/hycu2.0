#ifndef HYCU_NTT_HH
#define HYCU_NTT_HH

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace HyCu {

// Exact radix-2 transform. All operands are canonical residues modulo mod.
// Plans are immutable; transform buffers belong to the caller.
class NttPlan {
public:
  static constexpr std::uint32_t mod = 998244353;

  static std::uint32_t add(std::uint32_t a, std::uint32_t b) {
    const auto s = a + b;
    return s >= mod ? s - mod : s;
  }
  static std::uint32_t mul(std::uint32_t a, std::uint32_t b) {
    return std::uint64_t(a) * b % mod;
  }
  static std::uint32_t power(std::uint32_t a, std::uint32_t b) {
    std::uint32_t r = 1;
    for (; b; b >>= 1, a = mul(a, a)) if (b & 1) r = mul(r, a);
    return r;
  }

  explicit NttPlan(std::size_t minimum) : size_(1) {
    if (!minimum || minimum > (1u << 23))
      throw std::invalid_argument("Unsupported NTT length");
    while (size_ < minimum) size_ <<= 1;
    reverse_.resize(size_);
    for (std::size_t i = 1; i < size_; ++i)
      reverse_[i] = (reverse_[i >> 1] >> 1) | ((i & 1) ? size_ >> 1 : 0);
    roots_.resize(size_);
    root_quotients_.resize(size_);
    inverse_roots_.resize(size_);
    inverse_quotients_.resize(size_);
    for (std::size_t half = 1; half < size_; half <<= 1) {
      const auto w = power(3, (mod - 1) / (2 * half));
      const auto inverse_w = power(w, mod - 2);
      roots_[half] = 1;
      inverse_roots_[half] = 1;
      for (std::size_t j = 1; j < half; ++j) {
        roots_[half + j] = mul(roots_[half + j - 1], w);
        inverse_roots_[half + j] = mul(inverse_roots_[half + j - 1], inverse_w);
      }
      for (std::size_t j = 0; j < half; ++j) {
        root_quotients_[half + j] = (std::uint64_t(roots_[half + j]) << 32) / mod;
        inverse_quotients_[half + j] = (std::uint64_t(inverse_roots_[half + j]) << 32) / mod;
      }
    }
    inverse_size_ = power(size_, mod - 2);
  }

  std::size_t size() const { return size_; }
  // Reference transform with natural-order input and output.
  void transform(std::vector<std::uint32_t> & a, bool inverse = false) const {
    if (a.size() != size_) throw std::invalid_argument("NTT buffer length");
    for (std::size_t i = 0; i < size_; ++i)
      if (i < reverse_[i]) std::swap(a[i], a[reverse_[i]]);
    for (std::size_t half = 1; half < size_; half <<= 1)
      for (std::size_t start = 0; start < size_; start += 2 * half)
        for (std::size_t j = 0; j < half; ++j) {
          const auto u = a[start + j];
          const auto v = mul(a[start + j + half], roots_[half + j]);
          a[start + j] = add(u, v);
          a[start + j + half] = u >= v ? u - v : u + mod - v;
        }
    if (inverse) {
      std::reverse(a.begin() + 1, a.end());
      for (auto & v : a) v = mul(v, inverse_size_);
    }
  }

  // Forward DIF transform: natural-order input, bit-reversed spectrum.
  // Kernels and products use the same spectrum order.
  void forward_spectrum(std::vector<std::uint32_t> & a) const {
    if (a.size() != size_) throw std::invalid_argument("NTT buffer length");
    for (std::size_t half = size_ >> 1; half; half >>= 1)
      for (std::size_t start = 0; start < size_; start += 2 * half)
        for (std::size_t j = 0; j < half; ++j) {
          const auto u = a[start + j], v = a[start + j + half];
          a[start + j] = add(u, v);
          const auto difference = u >= v ? u - v : u + mod - v;
          a[start + j + half] = mul_fixed(difference, roots_[half + j], root_quotients_[half + j]);
        }
  }

  // DIT with inverse roots consumes the bit-reversed spectrum directly and
  // returns natural-order coefficients. Only final normalization remains.
  void inverse_spectrum(std::vector<std::uint32_t> & a) const {
    if (a.size() != size_) throw std::invalid_argument("NTT buffer length");
    for (std::size_t half = 1; half < size_; half <<= 1)
      for (std::size_t start = 0; start < size_; start += 2 * half)
        for (std::size_t j = 0; j < half; ++j) {
          const auto u = a[start + j];
          const auto v = mul_fixed(a[start + j + half], inverse_roots_[half + j], inverse_quotients_[half + j]);
          a[start + j] = add(u, v);
          a[start + j + half] = u >= v ? u - v : u + mod - v;
        }
    for (auto & v : a) v = mul(v, inverse_size_);
  }

private:
  // For a,w < mod and quotient=floor(w*2^32/mod), q estimates floor(a*w/mod)
  // with error at most one. Thus a*w-q*mod lies in [0,2*mod), below 2^32.
  // Unsigned low-word arithmetic is exact for that residue; one subtraction
  // gives the canonical representative. No floating point or lazy residues.
  static std::uint32_t mul_fixed(std::uint32_t a, std::uint32_t w, std::uint32_t quotient) {
    const auto q = std::uint32_t((std::uint64_t(a) * quotient) >> 32);
    const std::uint32_t remainder = a * w - q * mod;
    return remainder >= mod ? remainder - mod : remainder;
  }
  std::size_t size_;
  std::uint32_t inverse_size_;
  std::vector<std::size_t> reverse_;
  std::vector<std::uint32_t> roots_;
  std::vector<std::uint32_t> root_quotients_, inverse_roots_, inverse_quotients_;
};

} // namespace HyCu
#endif
