/**
* fixed size numerical vector (although general, designed for small vectors).
*
* Dan Israel Malta
**/
#pragma once

#include "FloatingPointTraits.h"
#include "Constants.h"
#include "static_for.h"
#include "Common.h"
#include <array>
#include <vector>
#include <functional>
#include <future>
#include <assert.h>
#include <string>

namespace Numeric {
    /*********************************/
    /* Vector Accessors/Constructors */
    /*********************************/

    // vector constructor type
    enum cRising   { rising };              // construct a vector whose elements are {0, 1, 2, ..., N-1}.
    enum cSequence { sequence };            // create a vector holding values from 'start' and continue with a given 'stride'.
                                            // (example: 'VectorN<int32_t, 4> a(Sequence, 3, -2)' will fill the vector in the following manner: 3, 1, -1, -3.
    enum cFibonacci { fibonacci  };         // build a Fibonacci series given two initial elements in the series.

    /**
    * \brief fixed size numerical vector
    *
    * @param {T, in} underlying type
    * @param {N, in} size
    **/
    template<typename T, std::size_t N> class VectorN {
        static_assert(std::is_arithmetic<T>::value, "VectorN<T,N> - T must be of numerical type.");

        // properties
    private:
        alignas(std::aligned_storage<sizeof(T), alignof(T)>::type) T m_data[N];

        // constructors
    public:

        // default constructor
        explicit constexpr VectorN() : m_data() {}

        // construct using a single value
        explicit constexpr VectorN(const T xi_value) {
            static_for<0, N>([&](std::size_t i) {
                m_data[i] = xi_value;
            });
        }

        // construct from a moveable array (usage: VectorN<float, 3> v{{1.0f, 2.0f, 3.0f}}; )
        explicit constexpr VectorN(T(&&xi_array)[N]) : m_data(std::make_move_iterator(std::begin(xi_array)), std::make_move_iterator(std::end(xi_array))) {}

        // construct from a parameter pack (usage: VectorN<int, 3> v(0, 1, 2); )
        template<class ... TS> explicit constexpr VectorN(T first, TS&&... values) : VectorN{ first, std::forward<T>(static_cast<T>(values))... } {}

        // construct a 'rising' vector
        explicit constexpr VectorN(cRising) {
            static_for<0, N>([&](std::size_t i) {
                m_data[i] = static_cast<T>(i);
            });
        }

        // construct a linear 'sequence' vector {xi_start, xi_start + 0 * xi_stride, xi_start + 1 * xi_stride, ..., xi_start + (N - 1) * xi_stride}
        explicit constexpr VectorN(cSequence, const T xi_start, const T xi_stride) {
            static_for<0, N>([&](std::size_t i) {
                m_data[i] = xi_start + static_cast<T>(i) * xi_stride;
            });
        }

        // construct a 'Fibonacci sequence' vector
        explicit constexpr VectorN(cFibonacci, const T xi_xo, const T xi_x1) {
            static_assert(N >= 2, "Fibonacci sequence requires at least two elements to be known.");

            m_data[0] = xi_xo;
            m_data[1] = xi_x1;

            static_for<2, N>([&](std::size_t i) {
                m_data[i] = m_data[i - 1] + m_data[i - 2];
            });
        }

        // construct using list initializer
        explicit constexpr VectorN(const std::initializer_list<T>&& xi_list) {
            // throws "array iterator + offset out of range" if more then N elements are entered
            std::move(std::begin(xi_list), std::end(xi_list), m_data);
        }

        // construct using std::array
        explicit constexpr VectorN(const std::array<T, N>& xi_arr) {
            std::move(std::begin(xi_arr), std::end(xi_arr), m_data);
        }

        explicit constexpr VectorN(std::array<T, N>&& xi_arr) {
            std::move(std::begin(xi_arr), std::end(xi_arr), m_data);
        }

        // construct using std::vector
        explicit constexpr VectorN(const std::vector<T>& xi_vec) {
            assert(N == xi_vec.size());
            std::move(std::begin(xi_vec), std::end(xi_vec), m_data);
        }

        explicit constexpr VectorN(std::vector<T>&& xi_vec) {
            assert(N == xi_vec.size());
            std::move(std::begin(xi_vec), std::end(xi_vec), m_data);
        }

        // copy semantics
        VectorN(const VectorN&) = default;
        VectorN& operator=(const VectorN&) = default;

        // move semantics
        VectorN(VectorN&&)            noexcept = default;
        VectorN& operator=(VectorN&&) noexcept = default;

        // assignment operations
    public:

        // assign an element
        constexpr VectorN& operator=(const T xi_value) noexcept {
            static_for<0, N>([&](std::size_t i) {
                m_data[i] = xi_value;
            });
            return *this;
        }

        // assign from a list
        constexpr VectorN& operator=(const std::initializer_list<T>& xi_list) {
            // throws "array iterator + offset out of range" if more then N elements are entered
            std::move(std::begin(xi_list), std::end(xi_list), m_data);
            return *this;
        }

        constexpr VectorN& operator=(std::initializer_list<T>&& xi_list) {
            // throws "array iterator + offset out of range" if more then N elements are entered
            std::move(std::begin(xi_list), std::end(xi_list), m_data);
            return *this;
        }

        // assign from array
        constexpr VectorN& operator=(const T(&&xi_array)[N]) {
            std::move(&xi_array[0], &xi_array[N], m_data);
            return *this;
        };

        // assign from std::array
        constexpr VectorN& operator=(const std::array<T, N>& xi_arr) {
            std::move(std::begin(xi_arr), std::end(xi_arr), m_data);
            return *this;
        }

        constexpr VectorN& operator=(std::array<T, N>&& xi_arr) {
            std::move(std::begin(xi_arr), std::end(xi_arr), m_data);
            xi_arr = nullptr;
            return *this;
        }

        // assign from std::vector
        constexpr VectorN& operator=(const std::vector<T>& xi_vec) {
            assert(N == xi_vec.size());
            std::move(std::begin(xi_vec), std::end(xi_vec), m_data);
            return *this;
        }

        constexpr VectorN& operator=(std::vector<T>&& xi_vec) {
            assert(N == xi_vec.size());
            std::move(std::begin(xi_vec), std::end(xi_vec), m_data);
            xi_vec = nullptr;
            return *this;
        }

        // iterators
    public:
        constexpr       T* begin()       { return &m_data[0]; }
        constexpr const T* begin() const { return &m_data[0]; }

        constexpr       T* end()       { return begin() + N; }
        constexpr const T* end() const { return begin() + N; }

        // set/get operations
    public:

        // '[]' element access
        constexpr T  operator[](const std::size_t i) const { return m_data[i]; }
        constexpr T& operator[](const std::size_t i)       { return m_data[i]; }

        // extract vector size
        constexpr std::size_t size() const noexcept { return N; }

        // get pointer to vector internal storage
        constexpr       T* data()       { return m_data; }
        constexpr const T* data() const { return m_data; }

        // --- special setters/getters for 2/3/4 elements vector ---

        template<typename = typename std::enable_if<N >= 1>::type> constexpr T& X()       noexcept { return m_data[0]; }
        template<typename = typename std::enable_if<N >= 1>::type> constexpr T  X() const noexcept { return m_data[0]; }

        template<typename = typename std::enable_if<N >= 2>::type> constexpr T& Y()       noexcept { return m_data[1]; }
        template<typename = typename std::enable_if<N >= 2>::type> constexpr T  Y() const noexcept { return m_data[1]; }

        template<typename = typename std::enable_if<N >= 3>::type> constexpr T& Z()       noexcept { return m_data[2]; }
        template<typename = typename std::enable_if<N >= 3>::type> constexpr T  Z() const noexcept { return m_data[2]; }

        template<typename = typename std::enable_if<N >= 4>::type> constexpr T& W()       noexcept { return m_data[3]; }
        template<typename = typename std::enable_if<N >= 4>::type> constexpr T  W() const noexcept { return m_data[3]; }

        template<typename = typename std::enable_if<N >= 1>::type> constexpr T& R()       noexcept { return m_data[0]; }
        template<typename = typename std::enable_if<N >= 1>::type> constexpr T  R() const noexcept { return m_data[0]; }

        template<typename = typename std::enable_if<N >= 2>::type> constexpr T& G()       noexcept { return m_data[1]; }
        template<typename = typename std::enable_if<N >= 2>::type> constexpr T  G() const noexcept { return m_data[1]; }

        template<typename = typename std::enable_if<N >= 3>::type> constexpr T& B()       noexcept { return m_data[2]; }
        template<typename = typename std::enable_if<N >= 3>::type> constexpr T  B() const noexcept { return m_data[2]; }

        template<typename = typename std::enable_if<N >= 4>::type> constexpr T& A()       noexcept { return m_data[3]; }
        template<typename = typename std::enable_if<N >= 4>::type> constexpr T  A() const noexcept { return m_data[3]; }

        // --- compile time 'swizzle like getters' ---

        // default is xyz() swizzling.
        // usage: VectorN<T, 3> zxy = vec.swizzle<2,0,1>();
        template<std::size_t i = 0, std::size_t j = 1, std::size_t k = 2, typename = typename std::enable_if<N >= 3>::type>
        constexpr VectorN<T, 3> swizzle() noexcept { return VectorN<T, 3>{ m_data[i], m_data[j], m_data[k] }; }

        template<std::size_t i = 0, std::size_t j = 1, std::size_t k = 2, typename = typename std::enable_if<N >= 3>::type>
        constexpr VectorN<T, 3> swizzle() const noexcept { return VectorN<T, 3>{ m_data[i], m_data[j], m_data[k] }; }

        // numerical assignment operator overloading
    public:

#define M_OPERATOR(OP)                                                 \
        constexpr VectorN& operator OP (const T xi_value) {            \
            static_for<0, N>([&](std::size_t i) {                      \
                m_data[i] OP xi_value;                                 \
            });                                                        \
            return *this;                                              \
        }                                                              \
        constexpr VectorN& operator OP (const VectorN& xi_vector) {    \
            assert(N == xi_vector.size());                             \
            static_for<0, N>([&](std::size_t i) {                      \
                m_data[i] OP xi_vector[i];                             \
            });                                                        \
            return *this;                                              \
        }                                                              \
        constexpr VectorN& operator OP (VectorN&& xi_vector) {         \
            assert(N == xi_vector.size());                             \
            static_for<0, N>([&](std::size_t i) {                      \
                m_data[i] OP std::move(xi_vector[i]);                  \
            });                                                        \
            return *this;                                              \
        }

        M_OPERATOR(-= );
        M_OPERATOR(+= );
        M_OPERATOR(*= );
        M_OPERATOR(/= );
        M_OPERATOR(&= );
        M_OPERATOR(|= );
        M_OPERATOR(^= );
        M_OPERATOR(>>= );
        M_OPERATOR(<<= );
#undef M_OPERATOR

        // stream operator overloading
    public:

        // output vector to stream
        friend std::ostream& operator<<(std::ostream& xio_stream, const VectorN& xi_vector) {
            xio_stream << "{";

            static_for<0, N - 1>([&](std::size_t i) {
                xio_stream << std::to_string(xi_vector[i]) << ", ";
            });
            xio_stream << xi_vector[N - 1] << "}";

            return xio_stream;
        }

        // read the space-separated components of a vector from a stream
        friend std::istream& operator>>(std::istream& is, const VectorN& xi_vector) {
            static_for<0, N>([&](std::size_t i) {
                is >> xi_vector[i];
            });

            return is;
        }

        // modifiers & queries
    public:

        // reverse a VectorN
        void Reverse() noexcept {
            static_for<0, N / 2>([&](std::size_t i) {
                std::swap(m_data[i], m_data[N - i - 1]);
            });
        }

        // rotate vector X steps to the left
        void RotateLeft(const std::size_t xi_steps) noexcept {
            std::rotate(&m_data[0], &m_data[xi_steps], &m_data[N]);
        }

        // sort VectorN according to a given predicate
        void Sort(std::function<bool(T, T)> xi_compare) noexcept {
            std::sort(&m_data[0], &m_data[N], xi_compare);
        }

        // test if VectorN is sorted according to a given predicate
        bool IsSorted(std::function<bool(T, T)> xi_compare) const noexcept {
            return std::is_sorted(&m_data[0], &m_data[N], xi_compare);
        }

        // returns the first element in VectorN which satisfy a certain predicate
        // (if an element was not found - returns the last element in vector)
        constexpr T FindIf(std::function<bool(const T)> xi_predicate) noexcept {
            std::size_t i{};
            bool predicate{ false };
            for (; (i < N) && !predicate; ++i) {
                predicate = xi_predicate(m_data[i]);
            }

            return (predicate) ? m_data[i - 1] : m_data[N - 1];
        }

        // randomly shuffle VectorN (using standard shuffling function)
        void Shuffle() {
            std::random_shuffle(&m_data[0], &m_data[N]);
        }

        // check if a unary predicate returns true for all elements in maVectorNtrix
        bool AllOf(std::function<bool(const T)> xi_predicate) const noexcept {
            return std::all_of(&m_data[0], &m_data[N], xi_predicate);
        }

        // check if a unary predicate returns true for at least one elements in VectorN
        bool AnyOf(std::function<bool(const T)> xi_predicate) const noexcept {
            return std::any_of(&m_data[0], &m_data[N], xi_predicate);
        }

        // check if a unary predicate returns true for no elements in VectorN
        bool NoneOf(std::function<bool(const T)> xi_predicate) const noexcept {
            return std::none_of(&m_data[0], &m_data[N], xi_predicate);
        }

        // return occurrence of a given value in a VectorN
        constexpr std::size_t Count(const T xi_value) const noexcept {
            return std::count(&m_data[0], &m_data[N], xi_value);
        }

        // element wise numerical operations
    public:

        /**
        * \brief limit all matrix elements into a given boundary
        *
        * @param {BaseMatrix, in|out} matrix to limit / limited matrix
        * @param {T,          in}     boundary
        * @param {T,          in}     boundary
        **/
        void Clamp(const T xi_boundary1, const T xi_boundary2) noexcept {
            static_for<0, N>([&](std::size_t i) {
                m_data[i] = Numeric::Clamp(m_data[i], xi_boundary1, xi_boundary2);
            });
        }

        /**
        * \brief limit all of matrix values into a symmetric (about zero) boundary
        *
        * @param {BaseMatrix, in|out} matrix to limit / limited matrix
        * @param {T,          in}     boundary
        **/
        void ClampSym(const T xi_boundary) noexcept {
            static_for<0, N>([&](std::size_t i) {
                m_data[i] = Numeric::ClampSym(m_data[i], xi_boundary);
            });
        }

        /**
        * \brief limit all of matrix values into a given boundary, in a circular motion
        *
        * @param {BaseMatrix, in|out} matrix to limit / limited matrix
        * @param {T,          in}     boundary lower value
        * @param {T,          in}     boundary upper value
        **/
        void ClampCircular(const T xi_min, const T xi_max) noexcept {
            static_for<0, N>([&](std::size_t i) {
                m_data[i] = Numeric::ClampSym(m_data[i], xi_min, xi_max);
            });
        }

        /**
        * \brief return the linear interpolation between two matrix
        *
        * @param {matrix, in}  matrix1
        * @param {matrix, in}  matrix2
        * @param {T,      in}  interpolation value [0, 1]
        * @param {matrix, out} linear interpolation from matrix1 to matrix2
        **/
        constexpr friend VectorN Lerp(const VectorN& xi_vec1, const VectorN& xi_vec2, const T& xi_interpolant) noexcept {
            assert(xi_vec1.size() == xi_vec2.size());
            assert(xi_interpolant >= T{});
            assert(xi_interpolant <= static_cast<T>(1));            

            // trivial cases
            if (FloatingPointTrait<T>::Equals(xi_interpolant, T{})) {
                return xi_vec1;
            }
            else if (FloatingPointTrait<T>::Equals(xi_interpolant, static_cast<T>(1))) {
                return xi_vec2;
            }

            VectorN xo_lerp;
            static_for<0, N>([&](std::size_t i) {
                xo_lerp[i] = Numeric::Lerp(xi_vec1[i], xi_vec2[i], xi_interpolant);
            });

            return xo_lerp;
        }

        // return the linear interpolation between two vectors
        constexpr friend VectorN Lerp(VectorN&& xi_vec1, VectorN&& xi_vec2, const T& xi_interpolant) noexcept {
            assert(xi_vec1.size() == xi_vec2.size());
            assert(xi_interpolant >= T{});
            assert(xi_interpolant <= static_cast<T>(1));

            // trivial cases
            if (FloatingPointTrait<T>::Equals(xi_interpolant, T{})) {
                return std::move(xi_vec1);
            }
            else if (FloatingPointTrait<T>::Equals(xi_interpolant, static_cast<T>(1))) {
                return std::move(xi_vec2);
            }

            static_for<0, N>([&](std::size_t i) {
                xi_vec1[i] = Numeric::Lerp(std::move(xi_vec1[i]), std::move(xi_vec2[i]), xi_interpolant);
            });

            return xi_vec1;
        }

        // return the minimal element-wise vector of several vectors
        constexpr friend VectorN Min(const VectorN& xi_vec) noexcept {
            return xi_vec;
        }

        constexpr friend VectorN Min(const VectorN& xi_vec1, const VectorN& xi_vec2) noexcept {
            assert(xi_vec1.size() == xi_vec2.size());
            VectorN xo_vec;

            static_for<0, N>([&](std::size_t i) {
                xo_vec[i] = Numeric::Min(xi_vec1[i], xi_vec2[i]);
            });

            return xo_vec;
        }

        constexpr friend VectorN Min(VectorN&& xi_vec1, VectorN&& xi_vec2) noexcept {
            assert(xi_vec1.size() == xi_vec2.size());

            static_for<0, N>([&](std::size_t i) {
                xi_vec1[i] = Numeric::Min(std::move(xi_vec1[i]), std::move(xi_vec2[i]));
            });

            return xi_vec1;
        }

        template<typename... TS> constexpr friend VectorN Min(const VectorN& xi_vec1, const VectorN& xi_vec2, const TS... args) noexcept {
            return Min(Min(xi_vec1, xi_vec2), args...);
        }

        template<typename... TS> constexpr friend VectorN Min(VectorN&& xi_vec1, VectorN&& xi_vec2, TS&&... args) noexcept {
            return Min(Min(xi_vec1, xi_vec2), args...);
        }

        // return the maximal element-wise vector of several vectors
        constexpr friend VectorN Max(const VectorN& xi_vec) noexcept {
            return xi_vec;
        }

        constexpr friend VectorN Max(const VectorN& xi_vec1, const VectorN& xi_vec2) noexcept {
            assert(xi_vec1.size() == xi_vec2.size());
            VectorN xo_vec;

            static_for<0, N>([&](std::size_t i) {
                xo_vec[i] = Numeric::Max(xi_vec1[i], xi_vec2[i]);
            });

            return xo_vec;
        }

        constexpr friend VectorN Max(VectorN&& xi_vec1, VectorN&& xi_vec2) noexcept {
            assert(xi_vec1.size() == xi_vec2.size());

            static_for<0, N>([&](std::size_t i) {
                xi_vec1[i] = Numeric::Max(std::move(xi_vec1[i]), std::move(xi_vec2[i]));
            });

            return xi_vec1;
        }

        template<typename... TS> constexpr friend VectorN Max(const VectorN& xi_vec1, const VectorN& xi_vec2, const TS... args) noexcept {
            return Max(Max(xi_vec1, xi_vec2), args...);
        }

        template<typename... TS> constexpr friend VectorN Max(VectorN&& xi_vec1, VectorN&& xi_vec2, TS&&... args) noexcept {
            return Max(Max(xi_vec1, xi_vec2), args...);
        }

        // return the sum of all elements in several vectors
        constexpr friend T Sum(const VectorN& xi_vec) noexcept {
            T xo_sum{};
            static_for<0, N>([&](std::size_t i) {
                xo_sum += xi_vec[i];
            });
            return xo_sum;
        }

        constexpr friend T Sum(VectorN&& xi_vec) noexcept {
            T xo_sum{};
            static_for<0, N>([&](std::size_t i) {
                xo_sum += std::move(xi_vec[i]);
            });
            return xo_sum;
        }

        template<typename... TS> constexpr friend T Sum(const VectorN& xi_vec1, const TS... args) noexcept {
            return Sum(xi_vec1) + Sum(args...);
        }

        template<typename... TS> constexpr friend T Sum(VectorN&& xi_vec1, const TS&&... args) noexcept {
            return Sum(xi_vec1) + Sum(args...);
        }

        // return the mean value of several vectors
        constexpr friend T Mean(const VectorN& xi_vec) noexcept {
            return Sum(xi_vec) / static_cast<T>(N);
        }

        constexpr friend T Mean(VectorN&& xi_vec) noexcept {
            return Sum(std::move(xi_vec)) / static_cast<T>(N);
        }

        template<typename... TS> constexpr friend T Mean(const VectorN& xi_vec1, const TS... args) noexcept {
            return Mean(xi_vec1) + Mean(args...);
        }

        template<typename... TS> constexpr friend T Mean(const VectorN&& xi_vec1, const TS&&... args) noexcept {
            return Mean(xi_vec1) + Mean(args...);
        }

        // product all the element in the VectorN
        constexpr friend T Prod(const VectorN& xi_vec) noexcept {
            T xo_prod{ 1 };
            static_for<0, N>([&](std::size_t i) {
                xo_prod *= xi_vec[i];
            });
            return xo_prod;
        }

        constexpr friend T Prod(VectorN&& xi_vec) noexcept {
            T xo_prod{ 1 };
            static_for<0, N>([&](std::size_t i) {
                xo_prod *= std::move(xi_vec[i]);
            });
            return xo_prod;
        }

        template<typename... TS> constexpr friend T Prod(const VectorN& xi_vec1, const TS... args) noexcept {
            return Prod(xi_vec1) * Prod(args...);
        }

        template<typename... TS> constexpr friend T Prod(VectorN&& xi_vec1, const TS&&... args) noexcept {
            return Prod(xi_vec1) * Prod(args...);
        }

        // dot product between two vectors
        constexpr friend T Dot(const VectorN& xi_a, const VectorN& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());

            T xo_dot{};
            static_for<0, N>([&](std::size_t i) {
                xo_dot += xi_a[i] * xi_b[i];
            });
            return xo_dot;
        }

        constexpr friend T Dot(VectorN&& xi_a, VectorN&& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());

            T xo_dot{};
            static_for<0, N>([&](std::size_t i) {
                xo_dot += std::move(xi_a[i]) * std::move(xi_b[i]);
            });
            return xo_dot;
        }

        // dot product between two vectors, on a given continuous range of indices
        constexpr friend T Dot(const VectorN& xi_a, const VectorN& xi_b, const std::size_t xi_index0, const std::size_t xi_index1) noexcept {
            assert(xi_index0 <= xi_index1);
            assert(xi_index1 <= xi_a.size());
            assert(xi_index1 <= xi_b.size());

            T xo_dot{};
            for (std::size_t i{ xi_index0 }; i < xi_index1; ++i) {
                xo_dot += xi_a[i] * xi_b[i];
            };
            return xo_dot;
        }

        constexpr friend T Dot(VectorN&& xi_a, VectorN&& xi_b, const std::size_t xi_index0, const std::size_t xi_index1) noexcept {
            assert(xi_index0 <= xi_index1);
            assert(xi_index1 <= xi_a.size());
            assert(xi_index1 <= xi_b.size());

            T xo_dot{};
            for (std::size_t i{ xi_index0 }; i < xi_index1; ++i) {
                xo_dot += std::move(xi_a[i]) * std::move(xi_b[i]);
            };
            return xo_dot;
        }

        // project first VectorN on second VectorN
        constexpr friend VectorN ProjectOn(const VectorN& xi_to, const VectorN& xi_on) noexcept {
            assert(xi_to.size() == xi_on.size());
            return (xi_on * (Dot(xi_to, xi_on) / Dot(xi_on, xi_on)));
        }

        friend VectorN ProjectOn(VectorN&& xi_to, VectorN&& xi_on) noexcept {
            assert(xi_to.size() == xi_on.size());
            return (xi_on * (Dot(xi_to, xi_on) / Dot(xi_on, xi_on)));
        }

        // euclidean distance between two VectorN
        constexpr friend T EuclideanDistanceBetween(const VectorN& xi_a, const VectorN& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            VectorN diff{ abs(xi_a - xi_b) };
            return diff.Norm();
        }

        friend T EuclideanDistanceBetween(VectorN&& xi_a, VectorN&& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            xi_a = abs(xi_a - xi_b);
            return xi_a.Norm();
        }

        // squared euclidean distance between two VectorN
        constexpr friend T SqueredEuclideanDistanceBetween(const VectorN& xi_a, const VectorN& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            VectorN diff{ abs(xi_a - xi_b) };
            return Dot(diff, diff);
        }

        friend T SqueredEuclideanDistanceBetween(VectorN&& xi_a, VectorN&& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            xi_a = abs(xi_a - xi_b);
            return Dot(xi_a, xi_a);
        }

        // return the Manhattan distance between two VectorN
        constexpr friend T ManhattanDistanceBetween(const VectorN& xi_a, const VectorN& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            VectorN diff{ abs(xi_a - xi_b) };
            return diff.Sum();
        }

        friend T ManhattanDistanceBetween(VectorN&& xi_a, VectorN&& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            xi_a = abs(xi_a - xi_b);
            return xi_a.Sum();
        }

        // return the Chebyshev distance between two VectorN
        constexpr friend T ChebyshevDistanceBetween(const VectorN& xi_a, const VectorN& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            VectorN diff{ abs(xi_a - xi_b) };
            return diff.MaxElement();
        }

        constexpr friend T ChebyshevDistanceBetween(VectorN&& xi_a, VectorN&& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            xi_a = abs(xi_a - xi_b);
            return xi_a.MaxElement();
        }

        // return the inverse Chebyshev distance between two VectorN
        constexpr friend T InverseChebyshevDistanceBetween(const VectorN& xi_a, const VectorN& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            VectorN diff{ abs(xi_a - xi_b) };
            return diff.MinElement();
        }

        constexpr friend T InverseChebyshevDistanceBetween(VectorN&& xi_a, VectorN&& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());
            xi_a = abs(xi_a - xi_b);
            return xi_a.MinElement();
        }

        // return the angle [radians] between two VectorN
        constexpr friend T Angle(const VectorN& xi_a, const VectorN& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());

            const T cosAngle{ Dot(xi_a, xi_b) / xi_a.Norm() / xi_b.Norm() };
            if (std::abs(cosAngle) <= static_cast<T>(1)) {
                return std::acos(cosAngle);
            }
            else {
                return std::numeric_limits<T>::quiet_NaN;
            }
        }

        constexpr friend T Angle(VectorN&& xi_a, VectorN&& xi_b) noexcept {
            assert(xi_a.size() == xi_b.size());

            const T cosAngle{ Dot(xi_a, xi_b) / xi_a.Norm() / xi_b.Norm() };
            if (std::abs(cosAngle) <= static_cast<T>(1)) {
                return std::acos(cosAngle);
            }
            else {
                return std::numeric_limits<T>::quiet_NaN;
            }
        }

        // general 'self' numerical operations
    public:
        // minimal element in the entire VectorN
        constexpr T MinElement() noexcept {
            return *std::min_element(&m_data[0], &m_data[N]);
        }       

        // maximal element in the entire VectorN
        constexpr T MaxElement() noexcept {
            return *std::max_element(&m_data[0], &m_data[N]);
        }

        // return VectorN norm (L2; its 'length')
        constexpr T Norm() noexcept {
            T xo_dot{};
            static_for<0, N>([&](std::size_t i) {
                xo_dot += m_data[i] * m_data[i];
            });

            return std::sqrt(xo_dot);
        }

        /**
        * \brief tests whether a vector is normalized.
        *
        * @param {double, in}  tolerance for normalization test (default is 2 * epsilon)
        * @param {bool,   out} true if squared length indicate a normalized algebric structure.
        **/
        constexpr bool IsNormalized(const T& xi_tol = static_cast<T>(2) * FloatingPointTrait<T>::epsilon()) const noexcept {
            return (std::abs(Norm() - static_cast<T>(1)) < static_cast<T>(2) * FloatingPointTrait<T>::epsilon());
        }

        // normalize VectorN
        void Normalize() noexcept {
            const T n{ Norm() };
            if (!FloatingPointTrait<T>::Equals(n, T{})) {
                const T nInv{ static_cast<T>(1) / n };
                static_for<0, N>([&](std::size_t i) {
                    m_data[i] *= nInv;
                });
            }
        }

        // Orthonormalize given VectorN using Gram-Schmidt process.
        void OrthoNormalize() noexcept {
            // normalize
            Normalize();

            static_for<0, N>([&](std::size_t i) {
                for (std::size_t j{}; j < i; ++j) {
                    const T dot{ m_data[i] * m_data[j] };
                    m_data[i] -= m_data[j] * dot;
                }
            });
        }
    };

    /**
    * numerical operator overload
    **/

    // unary minus
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> operator - (const VectorN<T, N>& xi_vec) {
        VectorN<T, N> xo_vec(xi_vec);
        xo_vec *= static_cast<T>(-1);
        return xo_vec;
    }

    template<typename T, std::size_t N> constexpr inline VectorN<T, N>& operator - (VectorN<T, N>&& xi_vec) {
        xi_vec *= static_cast<T>(-1);
        return xi_vec;
    }

#define M_OPERATOR(OP, AOP)                                                                                                                         \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> operator OP (const VectorN<T, N>& xi_lhs, const T& xi_rhs) {                 \
        VectorN<T, N> xo_vec(xi_lhs);                                                                                                               \
        xo_vec AOP xi_rhs;                                                                                                                          \
        return xo_vec;                                                                                                                              \
    }                                                                                                                                               \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N>& operator OP (VectorN<T, N>&& xi_lhs, const T& xi_rhs) {                     \
        xi_lhs AOP xi_rhs;                                                                                                                          \
        return xi_lhs;                                                                                                                              \
    }                                                                                                                                               \
    template<typename T, std::size_t N> constexpr  inline VectorN<T, N> operator OP (const T& xi_lhs, const VectorN<T, N>& xi_rhs) {                \
        VectorN<T, N> xo_vec(xi_rhs);                                                                                                               \
        xo_vec AOP xi_lhs;                                                                                                                          \
        return xo_vec;                                                                                                                              \
    }                                                                                                                                               \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N>& operator OP (const T& xi_lhs, VectorN<T, N>&& xi_rhs) {                     \
        xi_rhs AOP xi_lhs;                                                                                                                          \
        return xi_rhs;                                                                                                                              \
    }                                                                                                                                               \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> operator OP (const VectorN<T, N>& xi_lhs, const VectorN<T, N>& xi_rhs) {     \
        VectorN<T, N> xo_vec(xi_lhs);                                                                                                               \
        xo_vec AOP xi_rhs;                                                                                                                          \
        return xo_vec;                                                                                                                              \
    }                                                                                                                                               \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N>& operator OP (VectorN<T, N>&& xi_lhs, const VectorN<T, N>& xi_rhs) {         \
        xi_lhs AOP xi_rhs;                                                                                                                          \
        return xi_lhs;                                                                                                                              \
    }                                                                                                                                               \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N>& operator OP (const VectorN<T, N>& xi_lhs, VectorN<T, N>&& xi_rhs) {         \
        xi_rhs AOP xi_lhs;                                                                                                                          \
        return xi_rhs;                                                                                                                              \
    }                                                                                                                                               \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N>& operator OP (VectorN<T, N>&& xi_lhs, VectorN<T, N>&& xi_rhs) {              \
        xi_lhs AOP xi_rhs;                                                                                                                          \
        return xi_lhs;                                                                                                                              \
    }

    M_OPERATOR(+, +=);
    M_OPERATOR(-, -=);  // yes, it is possible to implement "scalar-vector', the result will be a vector
    M_OPERATOR(*, *=);
    M_OPERATOR(/, /=);
    M_OPERATOR(&, &=);
    M_OPERATOR(|, |=);
    M_OPERATOR(^, ^=);
    M_OPERATOR(>>, >>=);
    M_OPERATOR(<<, <<=);

#undef M_OPERATOR

    /**
    * relational operator overload
    **/
#define M_OPERATOR(OP)                                                                                                                   \
    template<typename T, std::size_t N> constexpr inline bool operator OP (const VectorN<T, N>& xi_lhs, const VectorN<T, N>& xi_rhs) {   \
        bool xo_rational{ false };                                                                                                       \
        for (std::size_t i{}; (i < N) && !xo_rational; ++i) {                                                                            \
            xo_rational = xi_lhs[i] OP xi_rhs[i];                                                                                        \
        }                                                                                                                                \
        return xo_rational;                                                                                                              \
    }                                                                                                                                    \
    template<typename T, std::size_t N> constexpr inline bool operator OP (const T xi_lhs, const VectorN<T, N>& xi_rhs) {                \
        bool xo_rational{ false };                                                                                                       \
        for (std::size_t i{}; (i < N) && !xo_rational; ++i) {                                                                            \
            xo_rational = xi_lhs OP xi_rhs[i];                                                                                           \
        }                                                                                                                                \
        return xo_rational;                                                                                                              \
    }                                                                                                                                    \
    template<typename T, std::size_t N> constexpr inline bool operator OP (const VectorN<T, N>& xi_lhs, const T& xi_rhs) {               \
        bool xo_rational{ false };                                                                                                       \
        for (std::size_t i{}; (i < N) && !xo_rational; ++i) {                                                                            \
            xo_rational = xi_lhs[i] OP xi_rhs;                                                                                           \
        }                                                                                                                                \
        return xo_rational;                                                                                                              \
    }

    M_OPERATOR(== );
    M_OPERATOR(!= );
    M_OPERATOR(>= );
    M_OPERATOR(> );
    M_OPERATOR(<= );
    M_OPERATOR(< );

#undef M_OPERATOR

    /**
    * numerical functions
    **/

    // unary functions
#define M_UNARY_FUNCTION(NAME, STL_NAME)                                                                          \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> NAME(const VectorN<T, N>& xi_value) {      \
        VectorN<T, N> xo_return(xi_value);                                                                        \
        static_for<0, N>([&](std::size_t i) {                                                                     \
            xo_return[i] = STL_NAME(xo_return[i]);                                                                \
        });                                                                                                       \
        return xo_return;                                                                                         \
    }                                                                                                             \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N>& NAME(VectorN<T, N>&& xi_value) {          \
        static_for<0, N>([&](std::size_t i) {                                                                     \
            xi_value[i] = STL_NAME(xi_value[i]);                                                                  \
        });                                                                                                       \
        return xi_value;                                                                                          \
    }

    M_UNARY_FUNCTION(abs, std::abs);
    M_UNARY_FUNCTION(floor, std::floor);
    M_UNARY_FUNCTION(ceil, std::ceil);
    M_UNARY_FUNCTION(round, std::round);
    M_UNARY_FUNCTION(rint, std::rint);
    M_UNARY_FUNCTION(trunc, std::trunc);
    M_UNARY_FUNCTION(sqrt, std::sqrt);
    M_UNARY_FUNCTION(cbrt, std::cbrt);
    M_UNARY_FUNCTION(exp, std::exp);
    M_UNARY_FUNCTION(exp2, std::exp2);
    M_UNARY_FUNCTION(log, std::log);
    M_UNARY_FUNCTION(log2, std::log2);
    M_UNARY_FUNCTION(log10, std::log10);
    M_UNARY_FUNCTION(log1p, std::log1p);
    M_UNARY_FUNCTION(sin, std::sin);
    M_UNARY_FUNCTION(cos, std::cos);
    M_UNARY_FUNCTION(tan, std::tan);
    M_UNARY_FUNCTION(asin, std::asin);
    M_UNARY_FUNCTION(acos, std::acos);
    M_UNARY_FUNCTION(atan, std::atan);
    M_UNARY_FUNCTION(sinh, std::sinh);
    M_UNARY_FUNCTION(cosh, std::cosh);
    M_UNARY_FUNCTION(tanh, std::tanh);
    M_UNARY_FUNCTION(asinh, std::asinh);
    M_UNARY_FUNCTION(acosh, std::acosh);
    M_UNARY_FUNCTION(atanh, std::atanh);
    M_UNARY_FUNCTION(erf, std::erf);
    M_UNARY_FUNCTION(erfc, std::erfc);
    M_UNARY_FUNCTION(tgamma, std::tgamma);
    M_UNARY_FUNCTION(lgamma, std::lgamma);

#undef M_UNARY_FUNCTION

    // binary functions
#define M_BINARY_FUNCTION(NAME, STL_NAME)                                                                                                      \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> NAME(const VectorN<T, N>& xi_lhs, const T xi_rhs) {                     \
        VectorN<T, N> xo_return();                                                                                                             \
        static_for<0, N>([&](std::size_t i) {                                                                                                  \
            xo_return[i] = STL_NAME(xi_lhs[i], xi_rhs);                                                                                        \
        });                                                                                                                                    \
        return xo_return;                                                                                                                      \
    }                                                                                                                                          \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N>& NAME(VectorN<T, N>&& xi_lhs, const T xi_rhs) {                         \
        static_for<0, N>([&](std::size_t i) {                                                                                                  \
            xi_lhs[i] = STL_NAME(xi_lhs[i], xi_rhs);                                                                                           \
        });                                                                                                                                    \
        return xi_lhs;                                                                                                                         \
    }                                                                                                                                          \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> NAME(const VectorN<T, N>& xi_lhs, const VectorN<T, N>& xi_rhs) {        \
        VectorN<T, N> xo_return();                                                                                                             \
        static_for<0, N>([&](std::size_t i) {                                                                                                  \
            xo_return[i] = STL_NAME(xi_lhs[i], xi_rhs[i]);                                                                                     \
        });                                                                                                                                    \
        return xo_return;                                                                                                                      \
    }                                                                                                                                          \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> NAME(VectorN<T, N>&& xi_lhs, const VectorN<T, N>& xi_rhs) {             \
        static_for<0, N>([&](std::size_t i) {                                                                                                  \
            xi_lhs[i] = STL_NAME(xi_lhs[i], xi_rhs[i]);                                                                                        \
        });                                                                                                                                    \
        return xi_lhs;                                                                                                                         \
    }                                                                                                                                          \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> NAME(const VectorN<T, N>& xi_lhs, VectorN<T, N>&& xi_rhs) {             \
        static_for<0, N>([&](std::size_t i) {                                                                                                  \
            xi_rhs[i] = STL_NAME(xi_lhs[i], xi_rhs[i]);                                                                                        \
        });                                                                                                                                    \
        return xi_rhs;                                                                                                                         \
    }                                                                                                                                          \
    template<typename T, std::size_t N> constexpr inline VectorN<T, N> NAME(VectorN<T, N>&& xi_lhs, VectorN<T, N>&& xi_rhs) {                  \
        static_for<0, N>([&](std::size_t i) {                                                                                                  \
            xi_lhs[i] = STL_NAME(xi_lhs[i], xi_rhs[i]);                                                                                        \
        });                                                                                                                                    \
        return xi_lhs;                                                                                                                         \
    }

    M_BINARY_FUNCTION(pow, std::pow);
    M_BINARY_FUNCTION(hypot, std::hypot);
    M_BINARY_FUNCTION(atan2, std::atan2);
    M_BINARY_FUNCTION(reminder, std::reminder);
    M_BINARY_FUNCTION(fmod, std::fmod);

#undef M_BINARY_FUNCTION

    /**
    * type queries
    **/
    template<typename>                  struct is_VectorN : public std::false_type {};
    template<typename T, std::size_t N> struct is_VectorN<VectorN<T, N>> : public std::true_type {};
    template<typename U> constexpr bool isVectorN(const U&) { return is_VectorN<U>::value; }
}
