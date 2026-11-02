/*============================================================================

    This file is part of HyCu.

    A polynomial over an FqField, with coefficients stored from low to high.

===============================================================================*/

#ifndef HYCU_FLINT_FQ_POLYNOMIAL_HH
#define HYCU_FLINT_FQ_POLYNOMIAL_HH

#include <cstddef>
#include <initializer_list>
#include <vector>

#include "flint_fq.hh"


namespace hycu
{

class FqPolynomial
{
  public:
    typedef FqElement Coefficient;

    // Construct the zero polynomial over field.
    explicit FqPolynomial(const FqField & field);

    // Coefficients are ordered c_0, c_1, ..., c_d.
    FqPolynomial(const FqField & field,
                 const std::vector<Coefficient> & coefficients);
    FqPolynomial(const FqField & field,
                 std::initializer_list<Coefficient> coefficients);

    FqPolynomial(const FqPolynomial & other) = default;
    FqPolynomial(FqPolynomial && other) = default;
    FqPolynomial & operator=(FqPolynomial other);

    void swap(FqPolynomial & other);

    const FqField & field() const;
    bool is_zero() const;
    // The zero polynomial has no degree and causes a domain_error.
    std::size_t degree() const;

    std::size_t coefficient_count() const;
    const Coefficient & coefficient(std::size_t index) const;
    const std::vector<Coefficient> & coefficients() const;

    // Read-only contiguous storage of coefficient objects. The pointer is
    // invalidated when this polynomial is assigned to or destroyed.
    const Coefficient * coefficient_data() const;

    // Evaluate c_0 + c_1*x + ... + c_d*x^d using Horner's method.
    Coefficient evaluate(const FqElement & point) const;

  private:
    void assign_coefficients(const std::vector<Coefficient> & coefficients);
    void remove_trailing_zeroes();

    FqField field_;
    std::vector<Coefficient> coefficients_;
};

} // namespace hycu

#endif
