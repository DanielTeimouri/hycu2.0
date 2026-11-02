/*============================================================================

    This file is part of HyCu.

===============================================================================*/

#include "flint_fq_polynomial.hh"

#include <stdexcept>
#include <utility>


namespace hycu
{

FqPolynomial::
FqPolynomial(
    const FqField & field
    ) :
  field_(field),
  coefficients_(1, field.zero())
{}


FqPolynomial::
FqPolynomial(
    const FqField & field,
    const std::vector<Coefficient> & coefficients
    ) :
  field_(field)
{
  this->assign_coefficients(coefficients);
}


FqPolynomial::
FqPolynomial(
    const FqField & field,
    std::initializer_list<Coefficient> coefficients
    ) :
  field_(field)
{
  this->assign_coefficients(std::vector<Coefficient>(coefficients));
}


FqPolynomial &
FqPolynomial::
operator=(
    FqPolynomial other
    )
{
  this->swap(other);
  return *this;
}


void
FqPolynomial::
swap(
    FqPolynomial & other
    )
{
  using std::swap;
  swap(this->field_, other.field_);
  this->coefficients_.swap(other.coefficients_);
}


void
FqPolynomial::
assign_coefficients(
    const std::vector<Coefficient> & coefficients
    )
{
  const FqElement zero = this->field_.zero();

  if ( coefficients.empty() ) {
    this->coefficients_.push_back(zero);
    return;
  }

  this->coefficients_.reserve(coefficients.size());
  for ( std::size_t ix = 0; ix < coefficients.size(); ++ix ) {
    if ( !zero.is_compatible_with(coefficients[ix]) )
      throw std::invalid_argument(
          "FqPolynomial: coefficient belongs to a different field");

    FqElement coefficient(zero);
    coefficient = coefficients[ix];
    this->coefficients_.push_back(std::move(coefficient));
  }

  this->remove_trailing_zeroes();
}


void
FqPolynomial::
remove_trailing_zeroes()
{
  while ( this->coefficients_.size() > 1
          && this->coefficients_.back().is_zero() )
    this->coefficients_.pop_back();
}


const FqField &
FqPolynomial::
field()
  const
{
  return this->field_;
}


bool
FqPolynomial::
is_zero()
  const
{
  return this->coefficients_.size() == 1
      && this->coefficients_.front().is_zero();
}


std::size_t
FqPolynomial::
degree()
  const
{
  if ( this->is_zero() )
    throw std::domain_error("FqPolynomial: zero polynomial has no degree");
  return this->coefficients_.size() - 1;
}


std::size_t
FqPolynomial::
coefficient_count()
  const
{
  return this->coefficients_.size();
}


const FqPolynomial::Coefficient &
FqPolynomial::
coefficient(
    std::size_t index
    )
  const
{
  return this->coefficients_.at(index);
}


const std::vector<FqPolynomial::Coefficient> &
FqPolynomial::
coefficients()
  const
{
  return this->coefficients_;
}


const FqPolynomial::Coefficient *
FqPolynomial::
coefficient_data()
  const
{
  return this->coefficients_.data();
}


FqPolynomial::Coefficient
FqPolynomial::
evaluate(
    const FqElement & point
    )
  const
{
  if ( !this->coefficients_.front().is_compatible_with(point) )
    throw std::invalid_argument(
        "FqPolynomial: evaluation point belongs to a different field");

  FqElement result(this->coefficients_.back());
  for ( std::size_t ix = this->coefficients_.size() - 1; ix > 0; --ix ) {
    result *= point;
    result += this->coefficients_[ix - 1];
  }
  return result;
}

} // namespace hycu
