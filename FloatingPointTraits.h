/**
* Dedicated operations for floating points
*
* Dan Israel Malta
**/
#pragma once
#ifndef _FLOATING_POINT_TRAITS_H_
#define _FLOATING_POINT_TRAITS_H_
// clang-format off
#include <cmath>
#include <type_traits>

/**
* Precision when testing floats for equality (5 digits)
*/
#ifndef FLOAT_EQUALITY_PRECISION
#define FLOAT_EQUALITY_PRECISION 1.0e-5f
#endif

/**
* Precision when testing doubles for equality (14 digits)
*/
#ifndef DOUBLE_EQUALITY_PRECISION
#define DOUBLE_EQUALITY_PRECISION 1.0e-14
#endif

/**
* Precision when testing long doubles for equality (17 digits)
* (except for MSVC, since they are same as doubles (https://msdn.microsoft.com/en-us/library/9cx8xs15.aspx))
**/
#ifndef LONG_DOUBLE_EQUALITY_PRECISION
    #if !defined(_MSC_VER)
        #define LONG_DOUBLE_EQUALITY_PRECISION 1.0e-17
    #else
        #define LONG_DOUBLE_EQUALITY_PRECISION 1.0e-14
    #endif
#endif

/*************************/
/* Floating Point Traits */
/*************************/
namespace Numeric {

    /**
    * default numeric behavior (defined for integers, overridden later for floating point)
    **/
    template<typename T> struct FloatingPointTrait {
        FloatingPointTrait() = delete;

        // minimal difference between two values to be considered unequal.
        constexpr static T epsilon() { return T(1); };

        // '=='
        constexpr static bool Equals(const T xi_lhs, const T xi_rhs) {
            return (xi_lhs == xi_rhs);
        }
    };

    // float epsilon
    constexpr float FloatingPointTrait<float>::epsilon() {
        return FLOAT_EQUALITY_PRECISION;
    }

    // double epsilon
    constexpr double FloatingPointTrait<double>::epsilon() {
        return DOUBLE_EQUALITY_PRECISION;
    }

    // long double epsilon
    constexpr long double FloatingPointTrait<long double>::epsilon() {
        return LONG_DOUBLE_EQUALITY_PRECISION;
    }

    // floating point equality test using absolute error and relative error boundary test
    constexpr bool FloatingPointTrait<float>::Equals(const float xi_lhs, const float xi_rhs) {
        // binary equal
        if (xi_lhs == xi_rhs) {
            return true;
        }

        // locals
		const float absLhs{ std::abs(xi_lhs) },
                    absRhs{ std::abs(xi_rhs) },
                    diff{ std::abs(xi_lhs - xi_rhs) };

        // absolute error test (since one or both of inputs are zero)
        if ((xi_lhs == 0) || (xi_rhs == 0) || (diff < FloatingPointTrait<float>::epsilon())) {
            return (diff < FloatingPointTrait<float>::epsilon());
        }

        // relative error
        return (diff / (absLhs + absRhs) < FloatingPointTrait<float>::epsilon());
    }

    // double equality test using absolute error and relative error boundary test
    constexpr bool FloatingPointTrait<double>::Equals(const double xi_lhs, const double xi_rhs) {
        // binary equal
        if (xi_lhs == xi_rhs) {
            return true;
        }

        // locals
		const double absLhs{ std::abs(xi_lhs) },
			absRhs{ std::abs(xi_rhs) },
			diff{ std::abs(xi_lhs - xi_rhs) };

        // absolute error test (since one or both of inputs are zero)
        if ((xi_lhs == 0) || (xi_rhs == 0) || (diff < FloatingPointTrait<double>::epsilon())) {
            return (diff < FloatingPointTrait<double>::epsilon());
        }

        // relative error
        return (diff / (absLhs + absRhs) < FloatingPointTrait<double>::epsilon());
    }

    // long double equality test using absolute error and relative error boundary test
    constexpr bool FloatingPointTrait<long double>::Equals(const long double xi_lhs, const long double xi_rhs) {
        // binary equal
        if (xi_lhs == xi_rhs) {
            return true;
        }

        // locals
		const long double absLhs{ std::abs(xi_lhs) },
						  absRhs{ std::abs(xi_rhs) },
						  diff{ std::abs(xi_lhs - xi_rhs) };

        // absolute error test (since one or both of inputs are zero)
        if ((xi_lhs == 0) || (xi_rhs == 0) || (diff < FloatingPointTrait<long double>::epsilon())) {
            return (diff < FloatingPointTrait<long double>::epsilon());
        }

        // relative error
        return (diff / (absLhs + absRhs) < FloatingPointTrait<long double>::epsilon());
    }
};
// clang-format on
#endif // !_FLOATING_POINT_TRAITS_H_