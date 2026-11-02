/*============================================================================

    This file is part of HyCu.

===============================================================================*/

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "flint_fq.hh"


BOOST_AUTO_TEST_CASE( flint_fq_polynomial_quotient_arithmetic )
{
  // F_9 = F_3[X] / (X^2 + 1), since X^2 + 1 is irreducible over F_3.
  const hycu::FqField field(3, std::vector<ulong>{1, 0, 1});
  const hycu::FqElement x = field.x();
  const hycu::FqElement one = field.one();

  BOOST_CHECK_EQUAL(field.prime(), 3U);
  BOOST_CHECK_EQUAL(field.degree(), 2U);
  BOOST_CHECK(field.modulus_coefficients()
              == (std::vector<ulong>{1, 0, 1}));

  // X^2 = -1 = 2 in the quotient field.
  const hycu::FqElement x_squared = x * x;
  BOOST_CHECK(x_squared == field.element(2));

  // (X + 1)(X + 2) = X^2 + 2 = 1 in F_9.
  BOOST_CHECK((x + one) * (x + field.element(2)) == one);
  BOOST_CHECK(x * x.inverse() == one);
  BOOST_CHECK(x.pow(8) == one);
}


BOOST_AUTO_TEST_CASE( flint_fq_automatic_modulus_defines_the_quotient )
{
  const hycu::FqField field(5, 2);
  const std::vector<ulong> modulus = field.modulus_coefficients();

  BOOST_REQUIRE_EQUAL(modulus.size(), 3U);
  BOOST_CHECK_EQUAL(modulus.back(), 1U);
  BOOST_CHECK(field.element(modulus).is_zero());
}


BOOST_AUTO_TEST_CASE( flint_fq_reduces_input_and_exposes_coefficients )
{
  const hycu::FqField field(3, std::vector<ulong>{1, 0, 1});

  // The constructor accepts an arbitrary polynomial and reduces it modulo C.
  const hycu::FqElement value = field.element({4, 5, 1});
  BOOST_CHECK(value == field.element({0, 2}));

  const std::vector<ulong> coefficients = value.coefficients();
  BOOST_REQUIRE_EQUAL(coefficients.size(), 2U);
  BOOST_CHECK_EQUAL(coefficients[0], 0U);
  BOOST_CHECK_EQUAL(coefficients[1], 2U);

  BOOST_REQUIRE_EQUAL(value.coefficient_count(), 2U);
  BOOST_REQUIRE(value.coefficient_data() != 0);
  BOOST_CHECK_EQUAL(value.coefficient_data()[0], 0U);
  BOOST_CHECK_EQUAL(value.coefficient_data()[1], 2U);
  BOOST_CHECK(value.native_handle() != 0);
  BOOST_CHECK(value.native_context() == field.native_handle());
}


BOOST_AUTO_TEST_CASE( flint_fq_rejects_invalid_fields_and_operations )
{
  BOOST_CHECK_THROW(hycu::FqField(4, 2), std::invalid_argument);
  BOOST_CHECK_THROW(
      hycu::FqField(3, std::vector<ulong>{1, 1, 1}),
      std::invalid_argument);

  const hycu::FqField f9(3, std::vector<ulong>{1, 0, 1});
  const hycu::FqField f25(5, std::vector<ulong>{2, 0, 1});

  BOOST_CHECK_THROW(f9.one() + f25.one(), std::invalid_argument);
  BOOST_CHECK_THROW(f9.one() / f9.zero(), std::domain_error);
  BOOST_CHECK_THROW(f9.zero().inverse(), std::domain_error);
}
