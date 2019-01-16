/**
* An interval (an algebraic field defining a region - {infimum, supremum})
* Interval is called "empty" if its infimum' is greater then its supremum'.
* Interval is called "single" if its infimum' is exactly its 'supremum'.
*
* Remarks:
* > If a mathematical operation results in an empty interval, it will be of the form (std::numeric_limits<T>::max(), std::numeric_limits<T>::min()).
* > This module focuses on handling correctly empty intervals, unbounded intervals and attempts to always return 
*   the smallest enclosing interval. i.e.:
*   >> (1,  2) & (3, 4)        = (std::numeric_limits<T>::max(), std::numeric_limits<T>::min())
*   >> (1, 20) / (0, 1)        = (1, std::numeric_limits<T>::max())
*   >> (-x, 0) * (y, infinity) = (std::numeric_limits<T>::min(), 0)
* > A numerical operation on an empty interval shall not be performed, and the output shall be an empty interval
*
* Dan Israel Malta
**/

#pragma once

#include "Common.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <string>
#include <functional>

/************/
/* Interval */
/************/
namespace Numeric {

    template<typename T> class Interval {
        static_assert(std::is_arithmetic<T>::value, "Interval<T> - T must be of numerical type.");

        // properties
        private:
            T data[2];  // {interval lower bound, interval upper bound}

        // constructors
        public:

            // trivial constructor
            constexpr Interval() { data[0] = T{}; data[1] = T{}; };

            // component-wise constructors
            explicit constexpr Interval(const T xi_inf, const T xi_sup) { data[0] = xi_inf; data[1] = xi_sup; };
            explicit constexpr Interval(const T xi_value) : Interval(xi_value, xi_value) {}

            // construct from array
            explicit constexpr Interval(const T* xi_p) : Interval({ xi_p[0], xi_p[1] }) {}

            // construct using list initializer (throws "array iterator + offset out of range if more then 2 elements are entered")
            explicit constexpr Interval(std::initializer_list<T>&& xi_list) { std::move(std::begin(xi_list), std::end(xi_list), data); }

            // copy/move constructor
            Interval(const Interval&)     = default;
            Interval(Interval&&) noexcept = default;

            // copy/move assignment
            Interval& operator=(const Interval&)     = default;
            Interval& operator=(Interval&&) noexcept = default;

            // specialized copy assignments (based upon specialized constructors)
            constexpr Interval<T>& operator = (std::initializer_list<T>&& xi_list) {
                std::move(std::begin(xi_list), std::end(xi_list), data);
                return *this;
            }

        // modifiers
        public:

            // interval lower bound
            constexpr T& Inf()       noexcept { return data[0]; }
            constexpr T  Inf() const noexcept { return data[0]; }

            // interval upper bound
            constexpr T& Sup()       noexcept { return data[1]; }
            constexpr T  Sup() const noexcept { return data[1]; }

        // friends
        public:
            friend std::ostream& operator<<(std::ostream& xio_stream, const Interval& xi_int) {
                return xio_stream << "{" << xi_int.Inf() << ", " << xi_int.Sup() << "}";
            }
    };

    // if an interval is empty, invalidate it to the maximum
    template<typename T> constexpr void Validate(Interval<T>& i) noexcept {
        if (i.Inf() > i.Sup()) {
            i.Inf() = std::numeric_limits<T>::max();
            i.Sup() = std::numeric_limits<T>::min();
        }
    }

    /*********************/
    /* operator overload */
    /*********************/

    /**
    * '-' unary operator overload
    **/
    template<typename T> constexpr inline Interval<T> operator - (const Interval<T>& xi_interval) {        
        return Interval<T>(-xi_interval.Sup(), -xi_interval.Inf());
    }

    /**
    * '+' operator overload
    **/
    template<typename T> constexpr inline Interval<T> operator + (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) {
        return Interval<T>(xi_lhs.Inf() + xi_rhs.Inf(), xi_lhs.Sup() + xi_rhs.Sup());
    }

    template<typename T, typename U> constexpr inline Interval<T> operator + (const Interval<T>& xi_lhs, const U& xi_rhs) {
        return Interval<T>(xi_lhs.Inf() + static_cast<T>(xi_rhs), xi_lhs.Sup() + static_cast<T>(xi_rhs));
    }

    template<typename T, typename U> constexpr inline Interval<T> operator + (const U& xi_lhs, const Interval<T>& xi_rhs) {
        return Interval<T>(xi_rhs.Inf() + static_cast<T>(xi_lhs), xi_rhs.Sup() + static_cast<T>(xi_lhs));
    }

    /**
    * '-' operator overload
    **/
    template<typename T> constexpr inline Interval<T> operator - (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) {
        return Interval<T>(xi_lhs.Inf() - xi_rhs.Sup(), xi_lhs.Sup() - xi_rhs.Inf());
    }

    template<typename T, typename U> constexpr inline Interval<T> operator - (const Interval<T>& xi_lhs, const U& xi_rhs) {
        return Interval<T>(xi_lhs.Inf() - static_cast<T>(xi_rhs), xi_lhs.Sup() - static_cast<T>(xi_rhs));
    }

    /**
    * '*' operator overload
    **/
    template<typename T> constexpr Interval<T> operator * (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) {

        // xi_lhs entirely negative
        if (IsNegative(xi_lhs.Sup())) {
            return (IsNegative(xi_rhs.Sup()))                             ? Interval<T>(xi_lhs.Sup() * xi_rhs.Sup(), xi_lhs.Inf() * xi_rhs.Inf()) : // xi_rhs entirely negative
                   (IsNegative(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup())) ? Interval<T>(xi_lhs.Inf() * xi_rhs.Sup(), xi_lhs.Inf() * xi_rhs.Inf()) : // xi_rhs include zero
                                                                            Interval<T>(xi_lhs.Inf() * xi_rhs.Sup(), xi_lhs.Sup() * xi_rhs.Inf());  // xi_rhs.inf entirely positive
        } // xi_lhs include zero
        else if (IsNegative(xi_lhs.Inf()) && IsPositive(xi_rhs.Sup())) {
            return (IsNegative(xi_rhs.Sup()))                             ? Interval<T>(xi_lhs.Sup() * xi_rhs.Inf(), xi_lhs.Inf() * xi_rhs.Inf())               :   // xi_rhs entirely negative
                   (IsNegative(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup())) ? Interval<T>(Numeric::Min(xi_lhs.Inf() * xi_rhs.Sup(), xi_lhs.Sup() * xi_rhs.Inf()),
                                                                                        Numeric::Max(xi_lhs.Inf() * xi_rhs.Inf(), xi_lhs.Sup() * xi_rhs.Sup())) :   // xi_rhs include zero
                                                                            Interval<T>(xi_lhs.Inf() * xi_rhs.Sup(), xi_lhs.Sup() * xi_rhs.Sup());                  // xi_rhs.inf entirely positive
        }   // xi_lhs.inf entirely positive
        else {
            return (IsNegative(xi_rhs.Sup()))                             ? Interval<T>(xi_lhs.Sup() * xi_rhs.Inf(), xi_lhs.Inf() * xi_rhs.Sup()) : // xi_rhs entirely negative
                   (IsNegative(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup())) ? Interval<T>(xi_lhs.Sup() * xi_rhs.Inf(), xi_lhs.Sup() * xi_rhs.Sup()) : // xi_rhs include zero
                                                                            Interval<T>(xi_lhs.Inf() * xi_rhs.Inf(), xi_lhs.Sup() * xi_rhs.Sup());  // xi_rhs.inf entirely positive
        }
    }

    template<typename T, typename U> constexpr inline Interval<T> operator * (const Interval<T>& xi_lhs, const U& xi_rhs) {
        return Interval<T>(xi_lhs.Inf() * static_cast<T>(xi_rhs), xi_lhs.Sup() * static_cast<T>(xi_rhs));
    }

    template<typename T, typename U> constexpr inline Interval<T> operator * (const U& xi_lhs, const Interval<T>& xi_rhs) {
        return Interval<T>(xi_rhs.Inf() * static_cast<T>(xi_lhs), xi_rhs.Sup() * static_cast<T>(xi_lhs));
    }

    /**
    * '/' operator overload
    **/
    template<typename T> constexpr Interval<T> operator / (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) {
        // xi_lhs entirely negative
        if (IsNegative(xi_lhs.Sup())) {
            return (IsNegative(xi_rhs.Sup()))                             ? Interval<T>(xi_lhs.Sup() / xi_rhs.Inf(), xi_lhs.Inf() / xi_rhs.Sup()) :     // xi_rhs entirely negative
                   (IsZero(xi_rhs.Inf()) && IsZero(xi_rhs.Sup()))         ? Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min()) : // xi_rhs is the point zero
                   (IsNegative(xi_rhs.Inf()) && IsZero(xi_rhs.Sup()))     ? Interval<T>(xi_lhs.Sup() / xi_rhs.Inf(), std::numeric_limits<T>::max()) :   // xi_rhs 'supremum' is zero
                   (IsZero(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup()))     ? Interval<T>(std::numeric_limits<T>::min(), xi_lhs.Sup() / xi_rhs.Sup()) :   // xi_rhs 'infimum' is zero
                   (IsNegative(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup())) ? Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min()) : // xi_rhs includes zero
                                                                            Interval<T>(xi_lhs.Inf() / xi_rhs.Sup(), xi_lhs.Sup() / xi_rhs.Sup());      // xi_rhs.inf entirely positive
        } // xi_lhs include zero
        else if (IsNegative(xi_lhs.Inf()) && IsPositive(xi_rhs.Sup())) {
            return (IsNegative(xi_rhs.Sup()))                             ? Interval<T>(xi_lhs.Sup() / xi_rhs.Sup(), xi_lhs.Inf() / xi_rhs.Sup()) :     // xi_rhs entirely negative
                   (IsZero(xi_rhs.Inf()) && IsZero(xi_rhs.Sup()))         ? Interval<T>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max()) : // xi_rhs is the point zero
                   (IsNegative(xi_rhs.Inf()) && IsZero(xi_rhs.Sup()))     ? Interval<T>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max()) : // xi_rhs 'supremum' is zero
                   (IsZero(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup()))     ? Interval<T>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max()) : // xi_rhs 'infimum' is zero
                   (IsNegative(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup())) ? Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min()) : // xi_rhs includes zero
                                                                            Interval<T>(xi_lhs.Inf() / xi_rhs.Sup(), xi_lhs.Sup() / xi_rhs.Inf());      // xi_rhs.inf entirely positive
        }   // xi_lhs.inf entirely positive
        else {
            return (IsNegative(xi_rhs.Sup()))                             ? Interval<T>(xi_lhs.Sup() / xi_rhs.Sup(), xi_lhs.Inf() / xi_rhs.Inf()) :     // xi_rhs entirely negative
                   (IsZero(xi_rhs.Inf()) && IsZero(xi_rhs.Sup()))         ? Interval<T>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max()) : // xi_rhs is the point zero
                   (IsNegative(xi_rhs.Inf()) && IsZero(xi_rhs.Sup()))     ? Interval<T>(std::numeric_limits<T>::min(), xi_lhs.Inf() / xi_rhs.Inf()) :   // xi_rhs 'supremum' is zero
                   (IsZero(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup()))     ? Interval<T>(xi_lhs.Inf() / xi_rhs.Sup(), std::numeric_limits<T>::max()) :   // xi_rhs 'infimum' is zero
                   (IsNegative(xi_rhs.Inf()) && IsPositive(xi_rhs.Sup())) ? Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min()) : // xi_rhs includes zero
                                                                            Interval<T>(xi_lhs.Inf() / xi_rhs.Sup(), xi_lhs.Sup() / xi_rhs.Inf());      // xi_rhs.inf entirely positive
        }
    }

    template<typename T, typename U> constexpr inline Interval<T> operator / (const Interval<T>& xi_lhs, const U& xi_rhs) {
        T rhs = static_cast<T>(xi_rhs);
        return Interval<T>(xi_lhs.Inf() / rhs, xi_lhs.Sup() / rhs);
    }

    /**
    * '<' operator overload
    **/
    template<typename T> constexpr inline bool operator < (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept { return (xi_lhs.Sup() < xi_rhs.Inf()); }
    template<typename T> constexpr inline bool operator < (const Interval<T>& xi_lhs, const T&           xi_rhs) noexcept { return (xi_lhs.Sup() < xi_rhs); }

    /**
    * '<=' operator overload
    **/
    template<typename T> constexpr inline bool operator <= (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept { return (xi_lhs.Sup() <= xi_rhs.Inf()); }
    template<typename T> constexpr inline bool operator <= (const Interval<T>& xi_lhs, const T&           xi_rhs) noexcept { return (xi_lhs.Sup() <= xi_rhs); }

    /**
    * '>' operator overload
    **/
    template<typename T> constexpr inline bool operator > (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept { return (xi_lhs.Inf() > xi_rhs.Sup()); }
    template<typename T> constexpr inline bool operator > (const Interval<T>& xi_lhs, const T&           xi_rhs) noexcept { return (xi_lhs.Inf() > xi_rhs); }

    /**
    * '>=' operator overload
    **/
    template<typename T> constexpr inline bool operator >= (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept { return (xi_lhs.Inf() >= xi_rhs.Sup()); }
    template<typename T> constexpr inline bool operator >= (const Interval<T>& xi_lhs, const T&           xi_rhs) noexcept { return (xi_lhs.Inf() >= xi_rhs); }

    /**
    * '==' operator overload
    **/
    template<typename T> constexpr inline bool operator == (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept { return ((xi_lhs.Inf() == xi_rhs.Inf()) && (xi_lhs.Sup() >= xi_rhs.Sup())); }
    template<typename T> constexpr inline bool operator == (const Interval<T>& xi_lhs, const T&           xi_rhs) noexcept { return ((xi_lhs.Inf() == xi_rhs)     && (xi_lhs.Sup() >= xi_rhs)); }

    /**
    * '!=' operator overload
    **/
    template<typename T> constexpr inline bool operator != (const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept { return !(xi_lhs == xi_rhs); }
    template<typename T> constexpr inline bool operator != (const Interval<T>& xi_lhs, const T&           xi_rhs) noexcept { return !(xi_lhs == xi_rhs); }

    /**
    * '<<' ostream operator overload
    **/
    template<typename T> constexpr inline std::ostream& operator << (std::ostream& xio_ostream, const Interval<T>& xi_int) {
        xio_ostream << "{" << std::to_string(xi_int.Inf()) << ", " << std::to_string(xi_int.Sup()) << "}";
        return xio_ostream;
    }

    /**********************/
    /* Boolean operations */
    /**********************/

    /**
    * \brief return union (or) of severals interval
    * 
    * @param {interval, in}  interval
    *            ....
    * @param {interval, in}  interval
    * @param {interval, out} union
    **/
    template<typename T> constexpr inline Interval<T> Union(const Interval<T>& xi_lhs) noexcept {
        if (IsEmpty(xi_lhs)) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());
        return xi_lhs;
    }

    template<typename T> constexpr inline Interval<T> Union(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept {
        if (IsEmpty(xi_lhs) || IsEmpty(xi_rhs)) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());

        Interval<T> xo_union(Numeric::Min(xi_lhs.Inf(), xi_rhs.Inf()),
                             Numeric::Max(xi_lhs.Sup(), xi_rhs.Sup()));
        Validate(xo_union);

        return xo_union;
    }

    template<typename T, typename... TS> constexpr  Interval<T> Union(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs, const TS... intervals) noexcept {
        return Union(Union(xi_lhs, xi_rhs), intervals...);
    }

    /**
    * \brief return intersection (and) of several intervals
    *
    * @param {interval, in}  interval
    *            ....
    * @param {interval, in}  interval
    * @param {interval, out} intersection
    **/
    template<typename T> constexpr inline Interval<T> Interesection(const Interval<T>& xi_lhs) noexcept {
        if (IsEmpty(xi_lhs)) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());
        return xi_lhs;
    }

    template<typename T> constexpr inline Interval<T> Interesection(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept {
        if (IsEmpty(xi_lhs) || IsEmpty(xi_rhs)) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());
        
        Interval<T> xo_intersection(Numeric::Max(xi_lhs.Inf(), xi_rhs.Inf()),
                                    Numeric::Min(xi_lhs.Sup(), xi_rhs.Sup()));
        Validate(xo_intersection);

        return xo_intersection;
    }

    template<typename T, typename... TS> constexpr Interval<T> Interesection(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs,  const TS... intervals) noexcept {
        return Interesection(Interesection(xi_lhs, xi_rhs), intervals...);
    }

    /**
    * \brief return the minimal interval enclosure
    * 
    * @param {interval, in}  interval
    *               ...
    * @param {interval, in}  interval
    * @param {interval, out} min
    **/
    template<typename T> constexpr inline Interval<T> IntervalMin(const Interval<T>& xi_lhs) noexcept {
        if (IsEmpty(xi_lhs)) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());
        return xi_lhs;
    }

    template<typename T> constexpr inline Interval<T> IntervalMin(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept {
        if (IsEmpty(xi_lhs) || IsEmpty(xi_rhs)) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());
        
        Interval<T> xo_min(Numeric::Min(xi_lhs.Inf(), xi_rhs.Inf()),
                           Numeric::Min(xi_lhs.Sup(), xi_rhs.Sup()));
        Validate(xo_min);

        return xo_min;
    }

    template<typename T, typename... TS> constexpr Interval<T> IntervalMin(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs, const TS... intervals) noexcept {
        return IntervalMin(IntervalMin(xi_lhs, xi_rhs), intervals...);
    }

    /**
    * \brief return the maximal interval enclosure
    * 
    * @param {interval, in}  interval
    * @param {interval, out} maximum enclosure
    **/
    template<typename T> constexpr inline Interval<T> IntervalMax(const Interval<T>& xi_lhs) noexcept {
        if (IsEmpty(xi_lhs)) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());
        return xi_lhs;
    }

    template<typename T> constexpr inline Interval<T> IntervalMax(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept {
        if (IsEmpty(xi_lhs) || IsEmpty(xi_rhs)) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());

        Interval<T> xo_max(Numeric::Max(xi_lhs.Inf(), xi_rhs.Inf()),
                           Numeric::Max(xi_lhs.Sup(), xi_rhs.Sup()));
        Validate(xo_max);

        return xo_max;
    }

    template<typename T, typename... TS> constexpr Interval<T> IntervalMax(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs, const TS... intervals) noexcept {
        return IntervalMax(IntervalMax(xi_lhs, xi_rhs), intervals...);
    }

    /***********/
    /* Queries */
    /***********/

    /**
    * \brief test if interval is empty
    *
    * @param {interval, in}  interval
    * @param {bool,     out} true if interval is empty
    **/
    template<typename T> constexpr inline bool IsEmpty(const Interval<T>& xi_int) noexcept {
        return (xi_int.Inf() > xi_int.Sup());
    }

    /**
    * \brief test if interval is infinte
    *
    * @param {interval, in}  interval
    * @param {bool,     out} true if interval is infinte
    **/
    template<typename T> constexpr inline bool IsInfinte(const Interval<T>& xi_int) noexcept {
        return ((xi_int.Inf() == std::numeric_limits<T>::min()) && (xi_int.Sup() == std::numeric_limits<T>::max()));
    }

    /**
    * \brief test if interval is 1D
    *
    * @param {interval, in}  interval
    * @param {bool,     out} true if interval is 1D
    **/
    template<typename T> constexpr inline bool IsSingle(const Interval<T>& xi_int) noexcept {
        return FloatingPointTrait<T>::Equals(xi_int.Inf(), xi_int.Sup());
    }

    /**
    * \brief test if scalar is inside interval
    *
    * @param {T,        in}  scalar
    * @param {interval, in}  interval
    * @param {bool,     out} true if interval is empty
    **/
    template<typename T> constexpr inline bool IsMember(const T& xi_value, const Interval<T>& xi_int) noexcept {
        return ((xi_value <= xi_int.Sup()) && (xi_value >= xi_int.Inf()));
    }

    /*********************/
    /* Interval specific */
    /*********************/

    /**
    * \brief return the middle point of interval (or  std::numeric_limits<T>::quiet_NaN() if interval is empty)
    *
    * @param {interval, in}  interval
    * @param {T,        out} middle point
    **/
    template<typename T> constexpr inline T Middle(const Interval<T>& xi_int) noexcept {
        // full interval
        if (IsInfinte(xi_int)) return std::numeric_limits<T>::max();
        
        // empty interval?
        if (IsEmpty(xi_int)) return std::numeric_limits<T>::quiet_NaN();

        return (static_cast<T>(0.5) * xi_int.Inf() + static_cast<T>(0.5) * xi_int.Sup());
    }

    /**
    * \brief return the width of interval (or std::numeric_limits<T>::quiet_NaN() if interval is empty, 
    *                                      or std::numeric_limits<T>::max() if it is the entire region)
    *
    * @param {interval, in}  interval
    * @param {T,        out} interval width
    **/
    template<typename T> constexpr inline T Width(const Interval<T>& xi_int) noexcept {
        // full interval
        if (IsInfinte(xi_int)) return std::numeric_limits<T>::max();
        
        // empty interval?
        if (IsEmpty(xi_int)) return std::numeric_limits<T>::quiet_NaN();

        return (xi_int.Sup() - xi_int.Inf());
    }

    /**
    * \brief the mignitude (minimum of absolute of infinmum/supremum) of interval (or std::numeric_limits<T>::quiet_NaN() if interval is empty,
    *                                                                              or std::numeric_limits<T>::max() if it is the entire region)
    *
    * @param {interval, in}  interval
    * @param {T,        out} interval mignitude
    **/
    template<typename T> constexpr inline T Mignitude(const Interval<T>& xi_int) noexcept {
        // full interval
        if (IsInfinte(xi_int)) return std::numeric_limits<T>::max();

        // empty interval?
        if (IsEmpty(xi_int)) return std::numeric_limits<T>::quiet_NaN();
        
        return Numeric::Min(std::abs(xi_int.Inf()), std::abs(xi_int.Sup()));
    }

    /**
    * \brief the magnitude (maximum of absolute of infinmum/supremum) of interval (or std::numeric_limits<T>::quiet_NaN() if interval is empty,
    *                                                                              or std::numeric_limits<T>::max() if it is the entire region)
    *
    * @param {interval, in}  interval
    * @param {T,        out} interval magnitude
    **/
    template<typename T> constexpr inline T Magnitude(const Interval<T>& xi_int) noexcept {
        // full interval
        if (IsInfinte(xi_int)) return std::numeric_limits<T>::max();

        // empty interval?
        if (IsEmpty(xi_int)) return std::numeric_limits<T>::quiet_NaN();
        
        return Numeric::Max(std::abs(xi_int.Inf()), std::abs(xi_int.Sup()));
    }

    /***************************/
    /* mathematical operations */
    /***************************/

    template<typename T> constexpr Interval<T> Abs(const Interval<T>& xi_int) noexcept {
        // negative interval
        if (xi_int.Sup() <= T{}) return Interval<T>(-xi_int.Sup(), -xi_int.Inf());
        
        // interval surrounds a zero
        if ((xi_int.Inf() <= T{}) && IsPositive(xi_int.Sup())) return Interval<T>(T{}, Numeric::Max(xi_int.Inf(), xi_int.Sup()));

        // positive interval
        return xi_int;
    }

    template<typename T> constexpr inline Interval<T> Round(const Interval<T>& xi_int) noexcept {
        return Interval<T>(std::round(xi_int.Inf()), std::round(xi_int.Sup())); 
    }

    template<typename T> constexpr inline Interval<T> Ceil(const Interval<T>& xi_int) noexcept {
        return Interval<T>(std::ceil(xi_int.Inf()), std::ceil(xi_int.Sup())); 
    }

    template<typename T> constexpr inline Interval<T> Floor(const Interval<T>& xi_int) noexcept {
        return Interval<T>(std::floor(xi_int.Inf()), std::floor(xi_int.Sup())); 
    }

    template<typename T> constexpr inline Interval<T> Sqrt(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int) || IsInfinte(xi_int)) return xi_int;
        
        // negative interval?
        if (IsNegative(xi_int.Inf()) || IsNegative(xi_int.Sup())) {
            return Interval<T>(std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN());
        }

        // positive interval
        return std::move(Interval<T>(std::sqrt(xi_int.Inf()), std::sqrt(xi_int.Sup())));
    }

    template<typename T> constexpr inline Interval<T> Hypot(const Interval<T>& xi_lhs, const Interval<T>& xi_rhs) noexcept {
        // magnitudes
        const T lhsMig{ Mignitude(xi_lhs) },
                lhsMag{ Magnitude(xi_lhs) },
                rhsMig{ Mignitude(xi_rhs) },
                rhsMag{ Magnitude(xi_rhs) };

        // output
        return Interval<T>(std::sqrt(lhsMig * lhsMig + rhsMig * rhsMig), 
                           std::sqrt(lhsMag * lhsMag + rhsMag * rhsMag));
    }

    template<typename T> constexpr Interval<T> Log(const Interval<T>& xi_int) noexcept {
        // negative interval?
        if (IsZero(xi_int.Sup())) return Interval<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::min());

        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;
        
        // transform interval to positive region [T{}, Infinity]
        const Interval<T> positiveInterval(T{}, std::numeric_limits<T>::max()),
                          intLocal( Interesection(xi_int, positiveInterval) );

        // output
        return Interval<T>(std::log(intLocal.Inf()), std::log(intLocal.Sup()));
    }

    template<typename T> constexpr inline Interval<T> Exp(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;
        
        return Interval<T>(exp(xi_int.Inf()), exp(xi_int.Sup()));
    }

    template<typename T> constexpr Interval<T> Sin(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;
        
        // housekeeping
        const T widthLocal{ Width(xi_int) },
                sinInf{ std::sin(xi_int.Inf()) },
                sinSup{ std::sin(xi_int.Sup()) },
                low{ (widthLocal >= Constants<T>::TAU()) ? static_cast<T>(-1.0) : Numeric::Min(sinInf, sinSup) },
                up{  (widthLocal >= Constants<T>::TAU()) ? static_cast<T>(1.0)  : Numeric::Max(sinInf, sinSup) };
        
        // derivative @ tau
        const int32_t cosSignLow = [&]() {
            if (!IsZero(std::cos(xi_int.Inf()))) return static_cast<int32_t>(std::copysign(1, std::cos(xi_int.Inf())));

            // crossover point
            return static_cast<int32_t>(std::copysign(1, low));
        }();

        // derivative @ 2*pi
        const int32_t cosSignUp = [&]() {
            if (!IsZero(std::cos(xi_int.Sup()))) return static_cast<int32_t>(std::copysign(T(1), std::cos(xi_int.Sup())));
            
            // crossover point
            return (-1 * static_cast<int32_t>(std::copysign(T(1), up)));
        }();

        // interval divergence
        if (((cosSignLow == -1) && (cosSignUp == 1)) || 
            ((cosSignLow == cosSignUp) && (widthLocal >= Constants<T>::PI()))) {
            low = static_cast<T>(-1.0);
        }

        if (((cosSignLow == -1) && (cosSignUp == -1)) || 
            ((cosSignLow == cosSignUp) && (widthLocal >= Constants<T>::PI()))) {
            up = static_cast<T>(-1.0);
        }

        // output
        return Interval<T>(low, up);
    }

    template<typename T> constexpr Interval<T> Cos(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;

        // housekeeping
        const T widthLocal{ Width(xi_int) },
                cosInf{ std::sin(xi_int.Inf()) },
                cosSup{ std::sin(xi_int.Sup()) },
                low{ (widthLocal >= Constants<T>::TAU()) ? static_cast<T>(-1.0) : Numeric::Min(cosInf, cosSup) },
                up{  (widthLocal >= Constants<T>::TAU()) ? static_cast<T>(1.0)  : Numeric::Max(cosInf, cosSup) };

        // derivative @ 2*pi
        const int32_t cosSignLow = [&]() {
            if (!IsZero(std::cos(xi_int.Inf()))) return static_cast<int32_t>(std::copysign(1, std::sin(xi_int.Inf())));

            // crossover point
            return (-1 * static_cast<int32_t>(std::copysign(1, low)));
        }();

        // derivative @ 2*pi
        const int32_t cosSignUp = [&]() {
            if (!IsZero(std::cos(xi_int.Sup()))) return static_cast<int32_t>(std::copysign(1, std::sin(xi_int.Sup())));
            
            // crossover point
            return static_cast<int32_t>(std::copysign(1, up));
        }();

        // interval divergence
        if ((((cosSignLow == -1) && (cosSignUp == 1)) ||
            ((cosSignLow == cosSignUp) && (widthLocal >= Constants<T>::PI()))) && ((!IsSingle(xi_int)) && xi_int.Inf() != T{})) {
            low = static_cast<T>(-1.0);
        }

        if (((cosSignLow == -1) && (cosSignUp == -1)) || ((cosSignLow == cosSignUp) && (widthLocal >= Constants<T>::PI()))) {
            up = static_cast<T>(-1.0);
        }

        // output
        return Interval<T>(low, up);
    }

    template<typename T> constexpr Interval<T> Tan(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;

        // housekeeping
        const T widthLocal{ Width(xi_int) }, 
                low{ (widthLocal >= Constants<T>::PI()) ? std::tan(xi_int.Inf()) : std::numeric_limits<T>::min() }, 
                up{ (widthLocal >= Constants<T>::PI()) ? std::tan(xi_int.Sup()) : std::numeric_limits<T>::max() };

        // divergence
        if ((low > up) ||
            (Max(std::abs(low), std::abs(up)) < static_cast<T>(1.0)) || ((widthLocal > Constants<T>::TAU()) &&
            (static_cast<int32_t>(std::copysign(1, low)) == static_cast<int32_t>(std::copysign(1, up))))) {
            low = std::numeric_limits<T>::min();
            up = std::numeric_limits<T>::max();
        }

        // output
        return Interval<T>(low, up);
    }

    template<typename T> constexpr Interval<T> Atan(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;

        Interval<T> xo_atan( std::atan(xi_int.Inf()), std::atan(xi_int.Sup()));
        xo_atan.Validate();

        return xo_atan;
    }

    template<typename T> constexpr Interval<T> Sinh(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;
        
        Interval<T> xo_sinh(std::sinh(xi_int.Inf()), std::sinh(xi_int.Sup()));
        xo_sinh.Validate();

        return xo_sinh;
    }

    template<typename T> constexpr Interval<T> Cosh(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;
        
        Interval<T> xo_cosh(std::cosh(xi_int.Inf()), std::cosh(xi_int.Sup()));
        xo_cosh.Validate();

        return xo_cosh;
    }

    template<typename T> constexpr Interval<T> Tanh(const Interval<T>& xi_int) noexcept {
        // empty interval?
        if (IsEmpty(xi_int)) return xi_int;
        
        Interval<T> xo_tanh(std::cosh(xi_int.Inf()), std::cosh(xi_int.Sup()));
        xo_tanh.Validate();

        return xo_tanh;
    }
};