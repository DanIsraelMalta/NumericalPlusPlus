/**
* Common numerical operations
*
* Dan Israel Malta
**/
#pragma once
#include "Constants.h"
#include "FloatingPointTraits.h"
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <numeric>
#include <functional>
#include <math.h>

/*********************/
/* Common Operations */
/*********************/
namespace Numeric {

	/**
    * \brief return true if value is positive (larger then zero)
    *
    * @param {T,    in}  value
    * @param {bool, out} true if value is positive
    **/
    template<typename T> constexpr inline bool IsPositive(const T xi_value) noexcept {
        return (xi_value >= FloatingPointTrait<T>::epsilon());
    }
    
    /**
    * \brief return true if value is negative (smaller then zero)
    *
    * @param {T,    in}  value
    * @param {bool, out} true if value is negative
    **/
    template<typename T> constexpr inline bool IsNegative(const T xi_value) noexcept {
        return (xi_value <= -FloatingPointTrait<T>::epsilon());
    }

	/**
    * \brief return true if value is zero (within floating point trait accuracy)
    *
    * @param {T,    in}  value
    * @param {bool, out} true if value is zero
    **/
    template<typename T> constexpr inline bool IsZero(const T xi_value) noexcept {
        return FloatingPointTrait<T>::Equals(xi_value, T{});
    }

    /**
    * \brief return the sign of a value
    *
    * @param {T, in}  value
    * @param {T, out} -1 if negative, otherwise: +1
    **/
    template<typename T> constexpr inline T Sign(const T xi_value) noexcept {
        return (xi_value < T{}) ? (T(-1)) : (T(1));
    }

    /**
    * \brief return the reminder of a division between two numbers
    *
    * @param {T, in}  a
    * @param {T, in}  b
    * @param {T, out} reminder of (a/b)
    **/
    template<typename T> constexpr inline T Remainder(const T xi_num, const T xi_den) noexcept {
        T xo_rem = static_cast<T>(std::fmod(static_cast<double>(xi_num), static_cast<double>(xi_den)));

        if (IsNegative(xo_rem)) {
            xo_rem += std::abs(xi_den);
        }

        return xo_rem;
    }

    /**
    * \brief limit a value into a given boundary
    *
    * @param {T, in}  value to limit
    * @param {T, in}  boundary
    * @param {T, in}  boundary
    * @param {T, out} limited value
    **/
    template<typename T> T constexpr inline Clamp(const T xio_value, const T xi_boundary1, const T xi_boundary2) noexcept {
        // boundary direction
        T low{ xi_boundary1 },
          high{ xi_boundary2 };
        if (high < low) {
            std::swap(high, low);
        }

        // clamp
        return Numeric::Min(high, Numeric::Max(low, xio_value));
    }

    /**
    * \brief limit a value into a symmetric (about zero) boundary
    *
    * @param {T, in}  value to limit
    * @param {T, in}  boundary
    * @param {T, out} limited value
    **/
    template<typename T> T constexpr inline ClampSym(const T xio_value, const T xi_boundary) noexcept {
        const T boundary{ std::abs(xi_boundary) };
        return Numeric::Min(boundary, Numeric::Max(-boundary, xio_value));
    }  

    /**
    * \brief limit a value into a given boundary, in a circular motion
    * 
    * @param {T, in}  value to limit
    * @param {T, in}  boundary lower value
    * @param {T, in}  boundary upper value
    * @param {T, out} limited value
    **/
    template<typename T> constexpr inline T ClampCircular(T xio_value, const T xi_min, const T xi_max) noexcept {
        return (xi_min + Remainder(xio_value - xi_min, xi_max - xi_min + T(1)));
    }

    /**
    * \brief normalize an angle to a 2*PI interval around a given center.
    *        example usage:
    *        NormalizeAngle(a, Constants::PI())  ->  [0,                 Constants::TAU()]
    *        NormalizeAngle(a, 0)                ->  [-Constants::PI(),  Constants::PI()]
    *        NormalizeAngle(a, -Constants::PI()) ->  [-Constants::TAU(), 0]
    *
    * @param {T, in}  angle to normalize
    * @param {T, in}  interval center
    * @param {T, out} normalized angle
    **/
    template<typename T> constexpr inline T NormalizeAngle(const T xi_angle, const T xi_center) noexcept {
        return (xi_angle - Constants<T>::TAU() * std::floor((xi_angle + Constants<T>::PI() - xi_center) / Constants<T>::TAU()));
    }

	/**
    * \brief return the minimal value among a list of values
    * 
    * @param {T, in}  values
    * @param {T, out} minimal value
    **/
	template<typename T> constexpr inline T Min(const T xi_a) noexcept {
		return (xi_a);
	}

    template<typename T> constexpr inline T Min(const T xi_a, const T xi_b) noexcept {
        return ((xi_a < xi_b) ? (xi_a) : (xi_b));
    }

    template<typename T, typename... TS> constexpr inline T Min(const T xi_a, const T xi_b, const TS... args) noexcept {
        return Min(Min(xi_a, xi_b), args...);
    }

    /**
    * \brief return the maximal value among a list of values
    *
    * @param {T, in}  values
    * @param {T, out} maximal value
    **/
	template<typename T> constexpr inline T Max(const T xi_a) noexcept {
		return (xi_a);
	}

    template<typename T> constexpr inline T Max(const T xi_a, const T xi_b) noexcept {
        return ((xi_a > xi_b) ? (xi_a) : (xi_b));
    }

    template<typename T, typename... TS> constexpr inline T Max(const T xi_a, const T xi_b, const TS... args) noexcept {
        return Max(Max(xi_a, xi_b), args...);
    }

	/**
    * \brief return the linear interpolation between two values
    *
    * @param {T, in}  value1
    * @param {T, in}  value2
    * @param {T, in}  interpolation value [0, 1]
    * @param {T, out} value1 * (1 - interpolator) + value2 * interpolator
    **/
    template<typename T> constexpr inline T Lerp(const T xi_lhs, const T xi_rhs, const T xi_interpolant) noexcept {
        return (xi_lhs * (static_cast<T>(1) - xi_interpolant) + xi_rhs * xi_interpolant);
    }

    /**
    * \brief stable numeric solution of a quadratic equation (a*x^2 + b*x + c = 0)
    * 
    * @param {T,    in}  a
    * @param {T,    in}  b
    * @param {T,    in}  c
    * @param {T,    out} smaller root (x1)
    * @param {T,    out} larger root (x2)
    * @param {bool, out} true if a solution exists - false otherwise
    **/
    template<typename T> constexpr bool SolveQuadratic(const T xi_a, const T xi_b, const T xi_c, T& xo_x1, T& xo_x2) noexcept {
        // trivial solution
        if (IsZero(xi_a) && IsZero(xi_b)) {
            xo_x1 = xo_x2 = T{};
            return true;
        }

        const T discriminant{ xi_b * xi_b - T(4) * xi_a * xi_c };

        if (IsNegative(discriminant)) {
            xo_x1 = T{};
            xo_x2 = T{};
            return false;
        }

        // solution
        const T t{ static_cast<T>(-0.5) * (xi_b + Sign(xi_b) * std::sqrt(discriminant)) };
        xo_x1 = t / xi_a;
        xo_x2 = xi_c / t;

        if (xo_x1 > xo_x2) {
            std::swap(xo_x1, xo_x2);
        }
        return true;
    }

    /**
    * \brief stable numeric solution of a cubic equation (x^3 + b*x^2 + c*x + d = 0)
    * 
    * @param {T,        in}  b
    * @param {T,        in}  c
    * @param {T,        in}  d
    * @param {uint32_t, out} number of real roots
    * @param {vector,   out} a 1x6 vector holding three paired solutions in the form (real solution #1, imag solution #1, ...)
    *                        if 1 real root exists: {real root 1, 0, Re(root 2), Im(root 2), Re(root 3), Im(root 3)}
    *                        if 3 real root exists: xo_roots[0] = real root 1, xo_roots[2] = real root 2, xo_roots[4] = real root 3
    **/
    template<typename T> constexpr uint32_t SolveCubic(const T xi_b, const T xi_c, const T xi_d, std::vector<T>& xo_roots) noexcept {
        std::fill(xo_roots.begin(), xo_roots.end(), T{});

        // transform to: x^3 + p*x + q = 0
        const T ov3{ static_cast<T>(1) / static_cast<T>(3) },
			    ov27{ static_cast<T>(1) / static_cast<T>(27) },
			    ovsqrt27{ static_cast<T>(1) / std::sqrt(static_cast<T>(27)) },
				bSqr{ xi_b * xi_b },
				p{ (static_cast<T>(3) * xi_c - bSqr) * ov3 },
				q{ (static_cast<T>(9) * xi_b * xi_c - static_cast<T>(27) * xi_d - static_cast<T>(2) * bSqr * xi_b) * ov27 };

        // x = w - (p / (3 * w))
        // (w^3)^2 - q*(w^3) - (p^3)/27 = 0
        T h{ q * q * static_cast<T>(0.25) + p * p * p  * ov27 };

        // one single real solution
        if (IsPositive(h)) {
            h = std::sqrt(h);

            const T qHalf{ q * static_cast<T>(0.5) },
				    bThird{ xi_b * ov3 },
					r{ qHalf + h },
					t{ qHalf - h },
					s{ std::cbrt(r) },
					u{ std::cbrt(t) },
					re{ -(s + u) * T(0.5) - bThird },
					im{  (s - u) * Constants<T>::SQRT3() * T(0.5) };

            // real root
            xo_roots[0] = (s + u) - bThird;

            // first complex root
            xo_roots[2] = re;
            xo_roots[3] = im;

            // second complex root
            xo_roots[4] = re;
            xo_roots[5] = -im;

            // one real solution
            return 1;
        }  // three real solutions
        else {            
            const T i{ p * std::sqrt(-p) * ovsqrt27 },     // p is negative (since h is positive)
					j{ std::cbrt(i) },
					k{ ov3 * std::acos((q / (static_cast<T>(2) * i))) },
					m{ std::cos(k) },
					n{ std::sin(k) * Constants<T>::SQRT3() },
					s{ -xi_b * ov3 };

            // roots
            xo_roots[0] = static_cast<T>(2) * j * m + s;
            xo_roots[2] = -j * (m + n) + s;
            xo_roots[4] = -j * (m - n) + s;

            // 3 real roots
            return 3;
        }
    }
};