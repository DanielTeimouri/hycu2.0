/*============================================================================

    This file is part of HyCu.

===============================================================================*/

#include "flint_fq.hh"

#include <limits>
#include <stdexcept>

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>


namespace
{

void
validate_field_parameters(
    ulong prime,
    std::size_t degree,
    const std::string & variable
    )
{
  if ( !n_is_prime(prime) )
    throw std::invalid_argument("FqField: characteristic must be prime");
  if ( degree == 0 )
    throw std::invalid_argument("FqField: extension degree must be positive");
  if ( degree > static_cast<std::size_t>(std::numeric_limits<slong>::max()) )
    throw std::invalid_argument("FqField: extension degree is too large for FLINT");
  if ( variable.empty() )
    throw std::invalid_argument("FqField: variable name must not be empty");
}


class NModPolynomial
{
  public:
    explicit NModPolynomial(ulong modulus)
    {
      nmod_poly_init(value_, modulus);
    }

    ~NModPolynomial()
    {
      nmod_poly_clear(value_);
    }

    nmod_poly_struct * get()
    {
      return value_;
    }

  private:
    NModPolynomial(const NModPolynomial &) = delete;
    NModPolynomial & operator=(const NModPolynomial &) = delete;

    nmod_poly_t value_;
};

} // namespace


namespace hycu
{

namespace detail
{

class FqContext
{
  public:
    FqContext(ulong prime, std::size_t degree, const std::string & variable)
    {
      validate_field_parameters(prime, degree, variable);
      fq_nmod_ctx_init_ui(value, prime, static_cast<slong>(degree),
                          variable.c_str());
    }

    FqContext(ulong prime,
              const std::vector<ulong> & modulus_coefficients,
              const std::string & variable)
    {
      if ( modulus_coefficients.size() < 2 )
        throw std::invalid_argument(
            "FqField: modulus must have positive degree");

      const std::size_t degree = modulus_coefficients.size() - 1;
      validate_field_parameters(prime, degree, variable);

      if ( modulus_coefficients.back() % prime != 1 )
        throw std::invalid_argument("FqField: modulus must be monic");

      NModPolynomial modulus(prime);
      for ( std::size_t ix = 0; ix < modulus_coefficients.size(); ++ix )
        nmod_poly_set_coeff_ui(modulus.get(), static_cast<slong>(ix),
                               modulus_coefficients[ix]);

      if ( !nmod_poly_is_irreducible(modulus.get()) )
        throw std::invalid_argument("FqField: modulus must be irreducible");

      fq_nmod_ctx_init_modulus(value, modulus.get(), variable.c_str());
    }

    ~FqContext()
    {
      fq_nmod_ctx_clear(value);
    }

    bool is_compatible_with(const FqContext & other) const
    {
      return this == &other
          || (fq_nmod_ctx_prime(value) == fq_nmod_ctx_prime(other.value)
              && fq_nmod_ctx_degree(value) == fq_nmod_ctx_degree(other.value)
              && nmod_poly_equal(fq_nmod_ctx_modulus(value),
                                 fq_nmod_ctx_modulus(other.value)));
    }

    fq_nmod_ctx_t value;

  private:
    FqContext(const FqContext &) = delete;
    FqContext & operator=(const FqContext &) = delete;
};

} // namespace detail


FqField::
FqField(
    Coefficient prime,
    std::size_t degree,
    const std::string & variable
    ) :
  context_(std::make_shared<detail::FqContext>(prime, degree, variable))
{}


FqField::
FqField(
    Coefficient prime,
    const std::vector<Coefficient> & modulus_coefficients,
    const std::string & variable
    ) :
  context_(std::make_shared<detail::FqContext>(
      prime, modulus_coefficients, variable))
{}


FqField::Coefficient
FqField::
prime()
  const
{
  return fq_nmod_ctx_prime(this->context_->value);
}


std::size_t
FqField::
degree()
  const
{
  return static_cast<std::size_t>(fq_nmod_ctx_degree(this->context_->value));
}


std::vector<FqField::Coefficient>
FqField::
modulus_coefficients()
  const
{
  const nmod_poly_struct * modulus =
      fq_nmod_ctx_modulus(this->context_->value);
  const std::size_t size = static_cast<std::size_t>(nmod_poly_length(modulus));

  std::vector<Coefficient> coefficients(size);
  for ( std::size_t ix = 0; ix < size; ++ix )
    coefficients[ix] = nmod_poly_get_coeff_ui(modulus, static_cast<slong>(ix));
  return coefficients;
}


FqElement
FqField::
zero()
  const
{
  return FqElement(*this);
}


FqElement
FqField::
one()
  const
{
  return FqElement(*this, 1);
}


FqElement
FqField::
x()
  const
{
  FqElement result(*this);
  fq_nmod_gen(result.value_, this->context_->value);
  return result;
}


FqElement
FqField::
element(
    Coefficient constant
    )
  const
{
  return FqElement(*this, constant);
}


FqElement
FqField::
element(
    const std::vector<Coefficient> & coefficients
    )
  const
{
  return FqElement(*this, coefficients);
}


FqElement
FqField::
element(
    std::initializer_list<Coefficient> coefficients
    )
  const
{
  return this->element(std::vector<Coefficient>(coefficients));
}


bool
FqField::
is_compatible_with(
    const FqField & other
    )
  const
{
  return this->context_->is_compatible_with(*other.context_);
}


const fq_nmod_ctx_struct *
FqField::
native_handle()
  const
{
  return this->context_->value;
}


FqElement::
FqElement(
    const FqField & field
    ) :
  context_(field.context_)
{
  fq_nmod_init(this->value_, this->context_->value);
}


FqElement::
FqElement(
    const FqField & field,
    Coefficient constant
    ) :
  context_(field.context_)
{
  fq_nmod_init(this->value_, this->context_->value);
  fq_nmod_set_ui(this->value_, constant, this->context_->value);
}


FqElement::
FqElement(
    const FqField & field,
    const std::vector<Coefficient> & coefficients
    ) :
  context_(field.context_)
{
  if ( coefficients.size()
       > static_cast<std::size_t>(std::numeric_limits<slong>::max()) )
    throw std::invalid_argument("FqElement: polynomial is too large for FLINT");

  fq_nmod_init(this->value_, this->context_->value);
  for ( std::size_t ix = 0; ix < coefficients.size(); ++ix )
    nmod_poly_set_coeff_ui(this->value_, static_cast<slong>(ix), coefficients[ix]);
  fq_nmod_reduce(this->value_, this->context_->value);
}


FqElement::
FqElement(
    const FqElement & other
    ) :
  context_(other.context_)
{
  fq_nmod_init(this->value_, this->context_->value);
  fq_nmod_set(this->value_, other.value_, this->context_->value);
}


FqElement::
FqElement(
    FqElement && other
    ) :
  context_(other.context_)
{
  fq_nmod_init(this->value_, this->context_->value);
  fq_nmod_swap(this->value_, other.value_, this->context_->value);
}


FqElement::
~FqElement()
{
  fq_nmod_clear(this->value_, this->context_->value);
}


FqElement &
FqElement::
operator=(
    const FqElement & other
    )
{
  if ( this != &other ) {
    this->require_compatible(other);
    fq_nmod_set(this->value_, other.value_, this->context_->value);
  }
  return *this;
}


FqElement &
FqElement::
operator=(
    FqElement && other
    )
{
  if ( this != &other ) {
    this->require_compatible(other);
    fq_nmod_swap(this->value_, other.value_, this->context_->value);
  }
  return *this;
}


void
FqElement::
require_compatible(
    const FqElement & other
    )
  const
{
  if ( !this->is_compatible_with(other) )
    throw std::invalid_argument("FqElement: operands belong to different fields");
}


FqElement &
FqElement::
operator+=(
    const FqElement & other
    )
{
  this->require_compatible(other);
  fq_nmod_add(this->value_, this->value_, other.value_, this->context_->value);
  return *this;
}


FqElement &
FqElement::
operator-=(
    const FqElement & other
    )
{
  this->require_compatible(other);
  fq_nmod_sub(this->value_, this->value_, other.value_, this->context_->value);
  return *this;
}


FqElement &
FqElement::
operator*=(
    const FqElement & other
    )
{
  this->require_compatible(other);
  fq_nmod_mul(this->value_, this->value_, other.value_, this->context_->value);
  return *this;
}


FqElement &
FqElement::
operator/=(
    const FqElement & other
    )
{
  this->require_compatible(other);
  if ( other.is_zero() )
    throw std::domain_error("FqElement: division by zero");
  fq_nmod_div(this->value_, this->value_, other.value_, this->context_->value);
  return *this;
}


FqElement
FqElement::
operator-()
  const
{
  FqElement result(*this);
  fq_nmod_neg(result.value_, result.value_, result.context_->value);
  return result;
}


FqElement
FqElement::
inverse()
  const
{
  if ( this->is_zero() )
    throw std::domain_error("FqElement: zero has no inverse");
  FqElement result(*this);
  fq_nmod_inv(result.value_, result.value_, result.context_->value);
  return result;
}


FqElement
FqElement::
pow(
    Coefficient exponent
    )
  const
{
  FqElement result(*this);
  fq_nmod_pow_ui(result.value_, result.value_, exponent, result.context_->value);
  return result;
}


bool
FqElement::
is_zero()
  const
{
  return fq_nmod_is_zero(this->value_, this->context_->value) != 0;
}


bool
FqElement::
is_one()
  const
{
  return fq_nmod_is_one(this->value_, this->context_->value) != 0;
}


FqElement::Coefficient
FqElement::
coefficient(
    std::size_t index
    )
  const
{
  if ( index >= static_cast<std::size_t>(
          fq_nmod_ctx_degree(this->context_->value)) )
    return 0;
  return nmod_poly_get_coeff_ui(this->value_, static_cast<slong>(index));
}


std::vector<FqElement::Coefficient>
FqElement::
coefficients()
  const
{
  const std::size_t size =
      static_cast<std::size_t>(fq_nmod_ctx_degree(this->context_->value));
  std::vector<Coefficient> result(size);
  for ( std::size_t ix = 0; ix < size; ++ix )
    result[ix] = this->coefficient(ix);
  return result;
}


const FqElement::Coefficient *
FqElement::
coefficient_data()
  const
{
  return this->value_->coeffs;
}


std::size_t
FqElement::
coefficient_count()
  const
{
  return static_cast<std::size_t>(this->value_->length);
}


std::size_t
FqElement::
coefficient_capacity()
  const
{
  return static_cast<std::size_t>(this->value_->alloc);
}


fq_nmod_struct *
FqElement::
native_handle()
{
  return this->value_;
}


const fq_nmod_struct *
FqElement::
native_handle()
  const
{
  return this->value_;
}


const fq_nmod_ctx_struct *
FqElement::
native_context()
  const
{
  return this->context_->value;
}


bool
FqElement::
is_compatible_with(
    const FqElement & other
    )
  const
{
  return this->context_->is_compatible_with(*other.context_);
}


FqElement
operator+(
    FqElement lhs,
    const FqElement & rhs
    )
{
  lhs += rhs;
  return lhs;
}


FqElement
operator-(
    FqElement lhs,
    const FqElement & rhs
    )
{
  lhs -= rhs;
  return lhs;
}


FqElement
operator*(
    FqElement lhs,
    const FqElement & rhs
    )
{
  lhs *= rhs;
  return lhs;
}


FqElement
operator/(
    FqElement lhs,
    const FqElement & rhs
    )
{
  lhs /= rhs;
  return lhs;
}


bool
operator==(
    const FqElement & lhs,
    const FqElement & rhs
    )
{
  return lhs.is_compatible_with(rhs)
      && fq_nmod_equal(lhs.value_, rhs.value_, lhs.context_->value) != 0;
}


bool
operator!=(
    const FqElement & lhs,
    const FqElement & rhs
    )
{
  return !(lhs == rhs);
}

} // namespace hycu
