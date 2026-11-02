/*============================================================================

    This file is part of HyCu.

===============================================================================*/

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "flint_fq_polynomial.hh"


BOOST_AUTO_TEST_CASE( flint_fq_polynomial_evaluates_over_prime_field )
{
  const hycu::FqField field(7, 1);
  const hycu::FqPolynomial polynomial(
      field,
      {field.element(1), field.element(1), field.zero(), field.zero(),
       field.zero(), field.element(1)});

  // f(X) = 1 + X + X^5. Thus f(2) = 35 = 0 and f(3) = 247 = 2 in F_7.
  BOOST_CHECK(polynomial.evaluate(field.element(2)).is_zero());
  BOOST_CHECK(polynomial.evaluate(field.element(3)) == field.element(2));
  for ( ulong value = 0; value < 7; ++value ) {
    const hycu::FqElement point = field.element(value);
    const hycu::FqElement expected = field.one() + point + point.pow(5);
    BOOST_CHECK(polynomial.evaluate(point) == expected);
  }

  BOOST_CHECK_EQUAL(polynomial.degree(), 5U);
  BOOST_REQUIRE_EQUAL(polynomial.coefficient_count(), 6U);
  BOOST_CHECK(polynomial.coefficient_data() != 0);
  BOOST_CHECK(polynomial.coefficient(5) == field.one());
}


BOOST_AUTO_TEST_CASE( flint_fq_polynomial_evaluates_over_extension_field )
{
  // F_9 = F_3[A] / (A^2 + 1), so A^2 = 2.
  const hycu::FqField field(3, std::vector<ulong>{1, 0, 1});
  const hycu::FqElement a = field.x();
  const hycu::FqPolynomial polynomial(
      field, {a, field.one(), a + field.one()});

  // f(A) = A + A + (A + 1)A^2 = A + 2.
  BOOST_CHECK(polynomial.evaluate(a) == field.element({2, 1}));
}


BOOST_AUTO_TEST_CASE( flint_fq_polynomial_normalises_zero_coefficients )
{
  const hycu::FqField field(5, 1);
  const hycu::FqPolynomial zero(field);
  const hycu::FqPolynomial constant(
      field, {field.element(3), field.zero(), field.zero()});

  BOOST_CHECK(zero.is_zero());
  BOOST_CHECK(zero.evaluate(field.element(4)).is_zero());
  BOOST_CHECK_THROW(zero.degree(), std::domain_error);

  BOOST_CHECK_EQUAL(constant.degree(), 0U);
  BOOST_CHECK_EQUAL(constant.coefficient_count(), 1U);
  BOOST_CHECK(constant.evaluate(field.element(4)) == field.element(3));
  BOOST_CHECK_THROW(constant.coefficient(1), std::out_of_range);
}


BOOST_AUTO_TEST_CASE( flint_fq_polynomial_rejects_incompatible_fields )
{
  const hycu::FqField f7(7, 1);
  const hycu::FqField f11(11, 1);
  const hycu::FqPolynomial polynomial(f7, {f7.one(), f7.one()});
  hycu::FqPolynomial reassigned(f11, {f11.one()});

  BOOST_CHECK_THROW(
      hycu::FqPolynomial(f7, {f7.one(), f11.one()}),
      std::invalid_argument);
  BOOST_CHECK_THROW(polynomial.evaluate(f11.one()), std::invalid_argument);

  // Assignment must safely replace both the field and its coefficients.
  reassigned = polynomial;
  BOOST_CHECK(reassigned.field().is_compatible_with(f7));
  BOOST_CHECK(reassigned.evaluate(f7.element(2)) == f7.element(3));
}
