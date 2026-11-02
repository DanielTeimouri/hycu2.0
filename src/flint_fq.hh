/*============================================================================

    This file is part of HyCu.

    A small C++ wrapper around FLINT's fq_nmod API.

===============================================================================*/

#ifndef HYCU_FLINT_FQ_HH
#define HYCU_FLINT_FQ_HH

#include <cstddef>
#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

#include <flint/fq_nmod.h>


namespace hycu
{

namespace detail
{
class FqContext;
}


class FqElement;


// The field F_p[X] / C(X). Copies of FqField share the same FLINT context.
class FqField
{
  public:
    typedef ulong Coefficient;

    // Let FLINT select a monic irreducible polynomial of the given degree.
    FqField(Coefficient prime, std::size_t degree,
            const std::string & variable = "T");

    // Coefficients of C(X), from constant to leading coefficient.
    FqField(Coefficient prime,
            const std::vector<Coefficient> & modulus_coefficients,
            const std::string & variable = "T");

    Coefficient prime() const;
    std::size_t degree() const;
    std::vector<Coefficient> modulus_coefficients() const;

    FqElement zero() const;
    FqElement one() const;
    // The residue class of X. It need not generate the multiplicative group.
    FqElement x() const;
    FqElement element(Coefficient constant) const;
    FqElement element(const std::vector<Coefficient> & coefficients) const;
    FqElement element(std::initializer_list<Coefficient> coefficients) const;

    bool is_compatible_with(const FqField & other) const;

    // Escape hatch for FLINT algorithms not yet covered by this wrapper.
    const fq_nmod_ctx_struct * native_handle() const;

  private:
    std::shared_ptr<detail::FqContext> context_;

    friend class FqElement;
};


class FqElement
{
  public:
    typedef FqField::Coefficient Coefficient;

    FqElement(const FqElement & other);
    FqElement(FqElement && other);
    ~FqElement();

    FqElement & operator=(const FqElement & other);
    FqElement & operator=(FqElement && other);

    FqElement & operator+=(const FqElement & other);
    FqElement & operator-=(const FqElement & other);
    FqElement & operator*=(const FqElement & other);
    FqElement & operator/=(const FqElement & other);

    FqElement operator-() const;
    FqElement inverse() const;
    FqElement pow(Coefficient exponent) const;

    bool is_zero() const;
    bool is_one() const;
    Coefficient coefficient(std::size_t index) const;
    std::vector<Coefficient> coefficients() const;

    // A read-only view of FLINT's dense, normalised nmod_poly storage.
    // Only coefficient_count() entries are meaningful. The pointer is
    // invalidated by any non-const operation on this element.
    const Coefficient * coefficient_data() const;
    std::size_t coefficient_count() const;
    std::size_t coefficient_capacity() const;

    fq_nmod_struct * native_handle();
    const fq_nmod_struct * native_handle() const;
    const fq_nmod_ctx_struct * native_context() const;

    bool is_compatible_with(const FqElement & other) const;

  private:
    explicit FqElement(const FqField & field);
    FqElement(const FqField & field, Coefficient constant);
    FqElement(const FqField & field,
              const std::vector<Coefficient> & coefficients);

    void require_compatible(const FqElement & other) const;

    std::shared_ptr<detail::FqContext> context_;
    fq_nmod_t value_;

    friend class FqField;
    friend bool operator==(const FqElement & lhs, const FqElement & rhs);
};


FqElement operator+(FqElement lhs, const FqElement & rhs);
FqElement operator-(FqElement lhs, const FqElement & rhs);
FqElement operator*(FqElement lhs, const FqElement & rhs);
FqElement operator/(FqElement lhs, const FqElement & rhs);

bool operator==(const FqElement & lhs, const FqElement & rhs);
bool operator!=(const FqElement & lhs, const FqElement & rhs);

} // namespace hycu

#endif
