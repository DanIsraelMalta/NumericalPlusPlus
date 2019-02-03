/**
* fixed size rectangular (i.e. ROW x COL) numerical matrix (although general, designed for small matrix).
*
* Dan Israel Malta
**/
#pragma once
#include "FloatingPointTraits.h"
#include "Constants.h"
#include "static_for.h"
#include "Common.h"
#include "VectorN.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>
#include <assert.h>
#include <algorithm>
#include <array>
#include <vector>
#include <type_traits>
#include <utility>
#include <limits>


namespace Numeric {

    /*****************/
    /* Matrix Layout */
    /*****************/

    // access matrix as row-major
    template<std::size_t ROW, std::size_t COL> struct RowMajor {
        static const bool m_RowMajor{ true };

        constexpr static std::size_t Index(const std::size_t xi_row, const std::size_t xi_col) {
            return (xi_row * COL + xi_col);
        }
    };

    // access matrix as column-major
    template<std::size_t ROW, std::size_t COL> struct ColumnMajor  {
        static const bool m_RowMajor{ false };   // column major...

        constexpr static std::size_t Index(const std::size_t xi_row, const std::size_t xi_col) {
            return (xi_col * COL + xi_row);
        }
    };

    // type trait to see if the property 'm_RowMajor' is included in a struct
    template<typename T, typename = void> struct has_RowMajor : std::false_type { };
    template<typename T> struct has_RowMajor<T, decltype(std::declval<T>().m_RowMajor, void())> : std::true_type { };

    // type trait to see if the method 'Index' is included in a struct
    template<typename T, typename = void> struct has_Index : std::false_type { };
    template<typename T> struct has_Index<T, decltype(std::declval<T>().Index, void())> : std::true_type { };

    /*********************************/
    /* Matrix Accessors/Constructors */
    /*********************************/

    // 1D matrix accessors
    struct Column { std::size_t m_Value; Column(const std::size_t& xi_value) : m_Value(xi_value) {}; }; // get a specific column from the matrix
    struct Row    { std::size_t m_Value; Row(const std::size_t& xi_value)    : m_Value(xi_value) {}; }; // get a specific row from the matrix
    enum aDiagonal { Diagonal };                                                                        // get matrix diagonal

    // 2D matrix accessor (syntactic sugar...)
    enum aLowerTriangular { LowerTriangular };  // get matrix lower triangular portion
    enum aUpperTriangular { UpperTriangular };  // get matrix upper triangular portion

    // matrix constructor type
    enum cIdentity     { Eye };             // construct an identity matrix
    enum cOuterProduct { Outerproduct };    // construct a matrix as outer product
    enum cRowWise      { RowWise };         // construct matrix row-wise
    enum cColumnWise   { ColumnWise };      // construct matrix column-wise
    enum cDiagonalWise { DiagonalWise };    // construct matrix diagonal-wise
    enum cVanDerMonde  { VanDerMonde };     // construct a Van-Der-Monde matrix
    enum cToeplitz     { Toeplitz };        // construct a Toeplitz matrix
    enum cAxisAngle    { AxisAngle };       // construct a rotation matrix (3x3) from an axis (normalized unit vector) and rotation angle (given by its sine & cosine components)

    /**
    * \brief fixed size numerical matrix
    *
    * @param {T,      in} underlying type
    * @param {ROW,    in} number of rows
    * @param {COL,    in} number of columns
    * @param {LAYOUT, in} matrix memory layout (row major or column major; by default it is row major)
    **/
    template<typename T, std::size_t ROW, std::size_t COL = ROW, class LAYOUT = RowMajor<ROW, COL>> class MatrixNM {
        static_assert(std::is_arithmetic<T>::value, "MatrixNM<T,ROW, COL> - T must be of numerical type.");
        static_assert(ROW != 0, "MatrixNM<T,ROW, COL> - ROW/COL parameters must be positive.");
        static_assert(COL != 0, "MatrixNM<T,ROW, COL> - ROW/COL parameters must be positive.");
        static_assert(std::is_class<LAYOUT>::value, "MatrixNM<T,ROW, COL, LAYOUT> - LAYOUT must be a struct of type 'RowMajor<ROW, COL>' or 'ColMajor<ROW, COL>'.");
        static_assert(has_RowMajor<LAYOUT>::value, "MatrixNM<T,ROW, COL, LAYOUT> - LAYOUT must include the property 'm_RowMajor'.");
        static_assert(has_Index<LAYOUT>::value, "MatrixNM<T,ROW, COL, LAYOUT> - LAYOUT must include the method 'Index'.");

        // properties
        private:
            const enum : bool { m_RowMajor = (LAYOUT::m_RowMajor) ? true : false }; // layout type
            const enum : std::size_t { SIZE = ROW * COL };                          // number of elements in matrix
            VectorN<T, SIZE> m_data;                                                // matrix data

        // constructors (from various collections/values)
        public:

            // default constructor
            explicit constexpr MatrixNM() : m_data() {}

            // construct using a single value
            explicit constexpr MatrixNM(const T xi_value) : m_data(xi_value) {}

            // construct from VectorN
            explicit constexpr MatrixNM(const VectorN<T, SIZE>& xi_vec) : m_data(xi_vec) {}

            // construct from a moveable array (usage: MatrixNM<float, 3, 2> v{{ 1.0f, 2.0f, 3.0f, 
            //                                                                   1.0f, 2.0f, 3.0f }}; )
            explicit constexpr MatrixNM(T(&&xi_array)[SIZE]) : m_data(std::make_move_iterator(std::begin(xi_array)), std::make_move_iterator(std::end(xi_array))) {
                // adjust matrix for ColumnMajor style
                if constexpr (!m_RowMajor) {
                    static_for<0, ROW>([&](std::size_t i) {
                        static_for<0, COL>([&](std::size_t j) {
                            std::swap(m_data[LAYOUT::Index(i, j)], m_data[LAYOUT::Index(j, i)]);
                        });
                    });
                }
            }

            // construct using list initializer
            // throws "array iterator + offset out of range" if more then SIZE elements are entered
            explicit constexpr MatrixNM(const std::initializer_list<T>&& xi_list) : m_data(std::move(xi_list)) {
                // adjust matrix for ColumnMajor style
                if constexpr (!m_RowMajor) {
                    static_for<0, ROW>([&](std::size_t i) {
                        static_for<0, COL>([&](std::size_t j) {
                            std::swap(m_data[LAYOUT::Index(i, j)], m_data[LAYOUT::Index(j, i)]);
                        });
                    });
                }
            }

            // construct using std::array
            explicit constexpr MatrixNM(const std::array<T, SIZE>& xi_arr) : m_data(xi_arr) {
                // adjust matrix for ColumnMajor style
                if constexpr (!m_RowMajor) {
                    static_for<0, ROW>([&](std::size_t i) {
                        static_for<0, COL>([&](std::size_t j) {
                            std::swap(m_data[LAYOUT::Index(i, j)], m_data[LAYOUT::Index(j, i)]);
                        });
                    });
                }
            }

            explicit constexpr MatrixNM(std::array<T, SIZE>&& xi_arr) : m_data(std::move(xi_arr)) {
                xi_arr = nullptr;

                // adjust matrix for ColumnMajor style
                if constexpr (!m_RowMajor) {
                    static_for<0, ROW>([&](std::size_t i) {
                        static_for<0, COL>([&](std::size_t j) {
                            std::swap(m_data[LAYOUT::Index(i, j)], m_data[LAYOUT::Index(j, i)]);
                        });
                    });
                }
            }

            // construct using std::vector
            explicit constexpr MatrixNM(const std::vector<T>& xi_vec) : m_data(xi_vec) {
                assert(SIZE == xi_vec.size());

                // adjust matrix for ColumnMajor style
                if constexpr (!m_RowMajor) {
                    static_for<0, ROW>([&](std::size_t i) {
                        static_for<0, COL>([&](std::size_t j) {
                            std::swap(m_data[LAYOUT::Index(i, j)], m_data[LAYOUT::Index(j, i)]);
                        });
                    });
                }
            }

            explicit constexpr MatrixNM(std::vector<T>&& xi_vec) : m_data(std::move(xi_vec)) {
                assert(SIZE == xi_vec.size());
                xi_vec = nullptr;

                // adjust matrix for ColumnMajor style
                if constexpr (!m_RowMajor) {
                    static_for<0, ROW>([&](std::size_t i) {
                        static_for<0, COL>([&](std::size_t j) {
                            std::swap(m_data[LAYOUT::Index(i, j)], m_data[LAYOUT::Index(j, i)]);
                        });
                    });
                }
            }

            // construct an I matrix
            explicit constexpr MatrixNM(cIdentity) : m_data(T{}) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    m_data[LAYOUT::Index(i, i)] = static_cast<T>(1);
                });
            }

            // construct the outer product from two vectors (usage: MatrixNM(Outerproduct, vec1, vec2))
            explicit constexpr MatrixNM(cOuterProduct, const VectorN<T, ROW>& xi_x, const VectorN<T, ROW>& xi_y) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, COL>([&](std::size_t j) {
                        m_data[LAYOUT::Index(i, j)] = xi_x[i] * xi_y[j];
                    });
                });
            }

            explicit constexpr MatrixNM(cOuterProduct, VectorN<T, ROW>&& xi_x, VectorN<T, ROW>&& xi_y) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, COL>([&](std::size_t j) {
                        m_data[LAYOUT::Index(i, j)] = std::move(xi_x[i]) * std::move(xi_y[j]);
                    });
                });
            }

            // construct a matrix such that all its rows are the same vector
            explicit constexpr MatrixNM(cRowWise, const VectorN<T, COL>& xi_vec) {
                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, COL>([&](std::size_t j) {
                        m_data[LAYOUT::Index(i, j)] = xi_vec[j];
                    });
                });
            }

            explicit constexpr MatrixNM(cRowWise, VectorN<T, COL>&& xi_vec) {
                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, COL>([&](std::size_t j) {
                        m_data[LAYOUT::Index(i, j)] = std::move(xi_vec[j]);
                    });
                });
            }

            // construct a matrix such that all its columns are the same vector
            explicit constexpr MatrixNM(cColumnWise, const VectorN<T, ROW>& xi_vec) {
                static_for<0, COL>([&](std::size_t i) {
                    static_for<0, ROW>([&](std::size_t j) {
                        m_data[LAYOUT::Index(j, i)] = xi_vec[j];
                    });
                });
            }

            explicit constexpr MatrixNM(cColumnWise, VectorN<T, ROW>&& xi_vec) {
                static_for<0, COL>([&](std::size_t i) {
                    static_for<0, ROW>([&](std::size_t j) {
                        m_data[LAYOUT::Index(j, i)] = std::move(xi_vec[j]);
                    });
                });
            }

            // construct a matrix such that its diagonal is same vector
            explicit constexpr MatrixNM(cDiagonalWise, const VectorN<T, ROW>& xi_vec) : m_data(T{}) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    m_data[LAYOUT::Index(i, i)] = xi_vec[i];
                });
            }

            explicit constexpr MatrixNM(cDiagonalWise, VectorN<T, ROW>&& xi_vec) : m_data(T{}) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    m_data[LAYOUT::Index(i, i)] = std::move(xi_vec[i]);
                });
            }

            // construct a Van-Der-Monde matrix
            explicit constexpr MatrixNM(cVanDerMonde, const VectorN<T, ROW>& xi_vec) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    for (std::size_t j{}; j < ROW; ++j) {
                        const T power{ static_cast<T>(ROW - j - 1) };
                        m_data[LAYOUT::Index(i, j)] = static_cast<T>(std::pow(xi_vec[i], power));
                    }
                });
            }

            explicit constexpr MatrixNM(cVanDerMonde, VectorN<T, ROW>&& xi_vec) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    for (std::size_t j{}; j < ROW; ++j) {
                        const T power{ static_cast<T>(ROW - j - 1) };
                        m_data[LAYOUT::Index(i, j)] = static_cast<T>(std::pow(std::move(xi_vec[i]), power));
                    }
                });
            }

            // construct a Toeplitz matrix
            explicit constexpr MatrixNM(cToeplitz, const VectorN<T, ROW>& xi_vec) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, ROW>([&](std::size_t j) {
                        m_data[LAYOUT::Index(i, j)] = (i >= j) ? xi_vec[i - j] : xi_vec[j - i];
                    });
                });
            }

            explicit constexpr MatrixNM(cToeplitz, VectorN<T, ROW>&& xi_vec) {
                static_assert(ROW == COL, "MatrixNM must be cubic.");

                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, ROW>([&](std::size_t j) {
                        m_data[LAYOUT::Index(i, j)] = (i >= j) ? std::move(xi_vec[i - j]) : std::move(xi_vec[j - i]);
                    });
                });
            }

            // construct a rotation matrix from an axis (normalized) and an angle (given by its sine & cosine components)
            explicit constexpr MatrixNM(cAxisAngle, const VectorN<T, 3>& xi_axis, const T xi_sine, const T xi_cosine) {
                static_assert(ROW == COL, "Rotation matrix must be cubic.");
                static_assert(ROW == 3, "Rotation matrix must be 3x3.");
                
                // locals
                const T oneMinusCosine{ T(1) - xi_cosine },
                        xx{ xi_axis.X() * xi_axis.X() },
                        xy{ xi_axis.X() * xi_axis.Y() },
                        xz{ xi_axis.X() * xi_axis.Z() },
                        yy{ xi_axis.Y() * xi_axis.Y() },
                        yz{ xi_axis.Y() * xi_axis.Z() },
                        zz{ xi_axis.Z() * xi_axis.Z() };

                if constexpr (m_RowMajor) {
                    m_data = { xi_cosine + xx * oneMinusCosine,             xy * oneMinusCosine + xi_axis.Z() * xi_sine, xz * oneMinusCosine - xi_axis.Y() * xi_sine,
                               xy * oneMinusCosine - xi_axis.Z() * xi_sine, xi_cosine + yy * oneMinusCosine,             yz * oneMinusCosine + xi_axis.X() * xi_sine,
                               xz * oneMinusCosine + xi_axis.Y() * xi_sine, yz * oneMinusCosine - xi_axis.X() * xi_sine, xi_cosine + zz * oneMinusCosine };
                }
                else {
                    m_data = { xi_cosine + xx * oneMinusCosine,             xy * oneMinusCosine - xi_axis.Z() * xi_sine, xz * oneMinusCosine + xi_axis.Y() * xi_sine,
                               xy * oneMinusCosine + xi_axis.Z() * xi_sine, xi_cosine + yy * oneMinusCosine,             yz * oneMinusCosine - xi_axis.X() * xi_sine,
                               xz * oneMinusCosine - xi_axis.Y() * xi_sine, yz * oneMinusCosine + xi_axis.X() * xi_sine, xi_cosine + zz * oneMinusCosine };
                }
            }
            
            // copy semantics
            MatrixNM(const MatrixNM&)            = default;
            MatrixNM& operator=(const MatrixNM&) = default;

            // move semantics
            MatrixNM(MatrixNM&&)            noexcept = default;
            MatrixNM& operator=(MatrixNM&&) noexcept = default;

        // assignment operator
        public:

            // assign an element
            constexpr MatrixNM& operator=(const T xi_value) noexcept {
                m_data = xi_value;
                return *this;
            }

            // assign a VectorN
            constexpr MatrixNM& operator=(const VectorN<T, SIZE>& xi_vec) noexcept {
                m_data = xi_vec;
                return *this;
            }

            constexpr MatrixNM& operator=(VectorN<T, SIZE>&& xi_vec) noexcept {
                m_data = std::move(xi_vec);
                return *this;
            }

            // assign from a list
            constexpr MatrixNM& operator=(const std::initializer_list<T>& xi_list) {
                // throws "array iterator + offset out of range" if more then N elements are entered
                m_data = xi_list;
                return *this;
            }

            constexpr MatrixNM& operator=(std::initializer_list<T>&& xi_list) {
                // throws "array iterator + offset out of range" if more then N elements are entered
                m_data = std::move(xi_list);
                return *this;
            }

            // assign from array
            constexpr MatrixNM& operator=(const T(&&xi_array)[SIZE]) {
                std::move(&xi_array[0], &xi_array[SIZE], m_data);
                return *this;
            };

            // assign from std::array
            constexpr MatrixNM& operator=(const std::array<T, SIZE>& xi_arr) {
                m_data = xi_arr;
                return *this;
            }

            constexpr MatrixNM& operator=(std::array<T, SIZE>&& xi_arr) {
                m_data = std::move(xi_arr);
                return *this;
            }

            // assign from std::vector
            constexpr MatrixNM& operator=(const std::vector<T>& xi_vec) {
                assert(SIZE == xi_vec.size());
                m_data = xi_vec;
                return *this;
            }

            constexpr MatrixNM& operator=(std::vector<T>&& xi_vec) {
                assert(SIZE == xi_vec.size());
                m_data = std::move(xi_vec);
                return *this;
            }

        // set/get operations
        public:

            // '[]' per-index element access (should not be used outside this header!)
            constexpr T  operator[](const std::size_t i) const { return m_data[i]; }
            constexpr T& operator[](const std::size_t i)         { return m_data[i]; }

            // '(row, col)' element access
            constexpr T  operator()(const std::size_t row, const std::size_t col) const { return m_data[LAYOUT::Index(row, col)]; }
            constexpr T& operator()(const std::size_t row, const std::size_t col)         { return m_data[LAYOUT::Index(row, col)]; }

            // get a specific row
            VectorN<T, COL> constexpr operator()(const Row& xi_row) const {
                VectorN<T, COL> xo_row;

                const std::size_t row{ xi_row.m_Value };
                static_for<0, COL>([&](std::size_t i) {
                    xo_row[i] = m_data[LAYOUT::Index(row, i)];
                });

                return xo_row;
            }

            // get a specific column
            VectorN<T, ROW> constexpr operator()(const Column& xi_col) const {
                VectorN<T, ROW> xo_col;

                const std::size_t col{ xi_col.m_Value };
                static_for<0, ROW>([&](std::size_t i) {
                    xo_col[i] = m_data[LAYOUT::Index(i, col)];
                });

                return xo_col;
            }

            // get diagonal
            VectorN<T, ROW> constexpr operator()(aDiagonal) const {
                static_assert(ROW == COL, "matrix must be cubic.");

                VectorN<T, ROW> xo_diag;

                static_for<0, ROW>([&](std::size_t i) {
                    xo_diag[i] = m_data[LAYOUT::Index(i, i)];
                });

                return xo_diag;
            }

            // get lower triangular
            MatrixNM<T, ROW, COL> constexpr operator()(aLowerTriangular) const {
                MatrixNM<T, ROW, COL> xo_lower(0.0f);

                static_for<0, ROW>([&](std::size_t i) {
                    for (std::size_t j{}; j <= i; ++j) {
                        xo_lower(i, j) = m_data[LAYOUT::Index(i, j)];
                    }
                });

                return xo_lower;
            }

            MatrixNM<T, ROW, COL> constexpr operator()(aLowerTriangular) {
                MatrixNM<T, ROW, COL> xo_lower(0.0f);

                static_for<0, ROW>([&](std::size_t i) {
                    for (std::size_t j{}; j <= i; ++j) {
                        xo_lower(i, j) = m_data[LAYOUT::Index(i, j)];
                    }
                });

                return xo_lower;
            }

            // get upper triangular
            MatrixNM<T, ROW, COL> constexpr operator()(aUpperTriangular) const {
                MatrixNM<T, ROW, COL> xo_upper(0.0f);

                static_for<0, ROW>([&](std::size_t i) {
                    for (std::size_t j{ i }; j < COL; ++j) {
                        xo_upper(i, j) = m_data[LAYOUT::Index(i, j)];
                    }
                });

                return xo_upper;
            }

            // get a block of matrix {row start, row end, column start, column end} (the block size is known at compile time)
            template<std::size_t ROW_MIN, std::size_t ROW_MAX, std::size_t COL_MIN, std::size_t COL_MAX>
            constexpr MatrixNM<T, ROW_MAX - ROW_MIN, COL_MAX - COL_MIN> GetRegion() const noexcept {
                static_assert(ROW_MIN < ROW_MAX, "GetRegion(): ROW_MIN < ROW_MAX.");
                static_assert(COL_MIN < COL_MAX, "GetRegion(): COL_MIN < COL_MAX.");
                static_assert(ROW_MAX <= ROW, "GetRegion(): ROW_MAX < ROW.");
                static_assert(COL_MAX <= COL, "GetRegion(): COL_MAX < COL.");

                MatrixNM<T, ROW_MAX - ROW_MIN, COL_MAX - COL_MIN> xo_block(T{});

                static_for<0, ROW_MAX - ROW_MIN>([&](std::size_t i) {
                    static_for<0, COL_MAX - COL_MIN>([&](std::size_t j) {
                        xo_block(i, j) = m_data[LAYOUT::Index(i + ROW_MIN, j + COL_MIN)];
                    });
                });

                return xo_block;
            }

            // --- special setters/getters for 2x2/3x3/4x4 elements matrix ---

            // 2x2
            template<typename = typename std::enable_if<(ROW == COL) && (ROW == 2)>::type> constexpr VectorN<T, 2> X() const noexcept { return VectorN<T, 2>{ m_data[LAYOUT::Index(0, 0)], 
                                                                                                                                                              m_data[LAYOUT::Index(0, 1)] }; }
            template<typename = typename std::enable_if<(ROW == COL) && (ROW == 2)>::type> constexpr VectorN<T, 2> Y() const noexcept { return VectorN<T, 2>{ m_data[LAYOUT::Index(1, 0)], 
                                                                                                                                                              m_data[LAYOUT::Index(1, 1)] }; }

            // 3x3 (extract Euler angles from rotation matrix)
            template<typename = typename std::enable_if<(ROW == COL) && (ROW == 3)>::type>
            constexpr VectorN<T, 3> Euler() const noexcept {
                return VectorN<T, 3>{ std::atan2(-this->operator()(1,2), this->operator()(2,2)) ,
                                      std::atan2(this->operator()(0,2),  std::hypot(this->operator()(1,2), this->operator()(2,2))),
                                      std::atan2(this->operator()(0,1),  this->operator()(0,0)) };
            }

            // 3x3, 4x4
            template<typename = typename std::enable_if<(ROW == COL) && ((ROW == 3) || (ROW == 4))>::type> constexpr VectorN<T, 3> X() const noexcept { return VectorN<T, 3>{ m_data[LAYOUT::Index(0, 0)],
                                                                                                                                                                              m_data[LAYOUT::Index(0, 1)],
                                                                                                                                                                              m_data[LAYOUT::Index(0, 2)] }; }

            template<typename = typename std::enable_if<(ROW == COL) && ((ROW == 3) || (ROW == 4))>::type> constexpr VectorN<T, 3> Y() const noexcept { return VectorN<T, 3>{ m_data[LAYOUT::Index(1, 0)],
                                                                                                                                                                              m_data[LAYOUT::Index(1, 1)],
                                                                                                                                                                              m_data[LAYOUT::Index(1, 2)] };
            }

            template<typename = typename std::enable_if<(ROW == COL) && ((ROW == 3) || (ROW == 4))>::type> constexpr VectorN<T, 3> Z() const noexcept { return VectorN<T, 3>{ m_data[LAYOUT::Index(2, 0)],
                                                                                                                                                                              m_data[LAYOUT::Index(2, 1)],
                                                                                                                                                                              m_data[LAYOUT::Index(2, 2)] };
            }

            // 4X4 (affine matrix translation vector (either last row or last column))
            template<typename = typename std::enable_if<(ROW == COL) && (ROW == 4)>::type> constexpr VectorN<T, 3> Translation() const noexcept { return VectorN<T, 3>{ m_data[LAYOUT::Index(3, 0)],
                                                                                                                                                                        m_data[LAYOUT::Index(3, 1)],
                                                                                                                                                                        m_data[LAYOUT::Index(3, 2)] };
            }

            // 4x4 (affine matrix scaling vector (its the diagonal of the 3x3 block))
            template<typename = typename std::enable_if<(ROW == COL) && (ROW == 4)>::type> constexpr VectorN<T, 3> Scale() const noexcept { return VectorN<T, 3>{ this->operator()(0, 0),
                                                                                                                                                                  this->operator()(1, 1),
                                                                                                                                                                  this->operator()(2, 2) };
            }
            
            // 4x4 (affine matrix transformation matrix)
            template<typename = typename std::enable_if<(ROW == COL) && (ROW == 4)>::type> constexpr MatrixNM<T, 3, 3> DCM() const noexcept { return this->GetRegion<0, 3, 0, 3>(); }
            
        // numerical assignment operator overloading
        public:

            // operations with scalars
#define M_OPERATOR(OP)                                           \
            constexpr MatrixNM& operator OP (const T xi_value) { \
                m_data OP xi_value;                              \
                return *this;                                    \
            }

            M_OPERATOR(+=);
            M_OPERATOR(-=);
            M_OPERATOR(*=);
            M_OPERATOR(/=);
            M_OPERATOR(&=);
            M_OPERATOR(|=);
            M_OPERATOR(^=);
            M_OPERATOR(>>=);
            M_OPERATOR(<<=);

#undef M_OPERATOR

            // operations with equally size matrix
#define M_OPERATOR(OP, AOP)                                                            \
            constexpr MatrixNM& operator OP (const MatrixNM<T, ROW, COL>& xi_mat) {    \
                static_for<0, SIZE>([&](std::size_t i) {                               \
                    m_data[i] AOP xi_mat[i];                                           \
                });                                                                    \
                return *this;                                                          \
            }                                                                          \
            constexpr MatrixNM& operator OP (MatrixNM<T, ROW, COL>&& xi_mat) {         \
                static_for<0, SIZE>([&](std::size_t i) {                               \
                    m_data[i] AOP std::move(xi_mat[i]);                                \
                });                                                                    \
                return *this;                                                          \
            }

            M_OPERATOR(+, +=);
            M_OPERATOR(-, -=);
            M_OPERATOR(&, &=);
            M_OPERATOR(| , |=);
            M_OPERATOR(^, ^=);
            M_OPERATOR(>> , >>=);
            M_OPERATOR(<< , <<=);

#undef M_OPERATOR

            // cubic matrix multiplication
            constexpr MatrixNM& operator *= (const MatrixNM<T, ROW, COL>& xi_mat) {
                static_assert(ROW == COL, "MatrixNM *= MatrixNM operations is only allowed when both matrix are cubic.");
                static_assert(m_RowMajor == xi_mat.IsRowMajor(), "MatrixNM *= MatrixNM operations is only allowed when both matrix layout is identical.");

                MatrixNM<T, ROW, COL>&& xo_mul;

                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, COL>([&](std::size_t j) {
                        T sum{};

                        static_for<0, ROW>([&](std::size_t ii) {
                            sum += m_data[LAYOUT::Index(i, ii)] * xi_mat(i, j);
                        });
                        xo_mul[LAYOUT::Index(i, j)] = sum;
                    });
                });

                *this = std::move(xo_mul);
                return *this;
            }

            constexpr MatrixNM& operator *= (MatrixNM<T, ROW, COL>&& xi_mat) {
                static_assert(ROW == COL, "MatrixNM *= MatrixNM operations is only allowed when both matrix are cubic.");
                static_assert(m_RowMajor == xi_mat.IsRowMajor(), "MatrixNM *= MatrixNM operations is only allowed when both matrix layout is identical.");

                MatrixNM<T, ROW, COL>&& xo_mul;

                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, COL>([&](std::size_t j) {
                        T sum{};

                        static_for<0, ROW>([&](std::size_t ii) {
                            sum += m_data[LAYOUT::Index(i, ii)] * std::move(xi_mat(i, j));
                        });
                        xo_mul[LAYOUT::Index(i, j)] = sum;
                    });
                });

                *this = std::move(xo_mul);
                return *this;
            }

        // stream operator overloading
        public:

            // output matrix elements to a stream
            friend std::ostream& operator<<(std::ostream& xio_stream, const MatrixNM& xi_mat) {
                xio_stream << "\n{\n";

                static_for<0, ROW>([&](std::size_t i) {
                    xio_stream << "{";
                    static_for<0, COL - 1>([&](std::size_t j) {
                        xio_stream << xi_mat(i, j) << ", ";
                    });
                    xio_stream << xi_mat(i, COL - 1) << "}\n";
                });

                return xio_stream << "\n}\n";
            }

        // element wise iterators
        public:

                  T* begin()       { return &m_data[0]; }
            const T* begin() const { return &m_data[0]; }

                  T* end()       { return begin() + SIZE; }
            const T* end() const { return begin() + SIZE; }         

        // modifiers
        public:

            // swap rows
            void SwapRows(const std::size_t a, const std::size_t b) {
                static_for<0, COL>([&](std::size_t i) {
                    std::swap(m_data[LAYOUT::Index(a, i)], m_data[LAYOUT::Index(b, i)]);
                });
            }

            // swap columns
            void SwapColumns(const std::size_t a, const std::size_t b) {
                static_for<0, ROW>([&](std::size_t i) {
                    std::swap(m_data[LAYOUT::Index(i, a)], m_data[LAYOUT::Index(i, b)]);
                });
            }

            // set row
            void SetRow(const std::size_t xi_inedx, const T xi_value) {
                static_for<0, COL>([&](std::size_t i) {
                    m_data[LAYOUT::Index(xi_inedx, i)] = xi_value;
                });
            }

            void SetRow(const std::size_t xi_inedx, const VectorN<T, COL> xi_vec) {
                static_for<0, COL>([&](std::size_t i) {
                    m_data[LAYOUT::Index(xi_inedx, i)] = xi_vec[i];
                });
            }

            // set column
            void SetColumn(const std::size_t xi_inedx, const T xi_value) {
                static_for<0, ROW>([&](std::size_t i) {
                    m_data[LAYOUT::Index(i, xi_inedx)] = xi_value;
                });
            }

            void SetColumn(const std::size_t xi_inedx, const VectorN<T, ROW> xi_vec) {
                static_for<0, ROW>([&](std::size_t i) {
                    m_data[LAYOUT::Index(i, xi_inedx)] = xi_vec[i];
                });
            }

            // return a transposed matrix
            constexpr MatrixNM<T, ROW, COL> Transpose() noexcept {
                MatrixNM<T, ROW, COL> xo_transposed;

                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, COL>([&](std::size_t j) {
                        xo_transposed(j, i) = m_data[LAYOUT::Index(i, j)];
                    });
                });

                return xo_transposed;
            }

            // return the (i,j) minor (matrix, not determinant value) of a given matrix
            template<std::size_t I, std::size_t J> constexpr MatrixNM<T, ROW - 1, COL - 1> Minor() noexcept {
                static_assert(I < ROW, "MatrixNM<T, COL, ROW>::Minor() index are out of bound");
                static_assert(J < COL, "MatrixNM<T, COL, ROW>::Minor() index are out of bound");

                MatrixNM<T, ROW - 1, COL - 1> xo_minor;

                // top left region
                static_for<0, I>([&](std::size_t i) {
                    static_for<0, J>([&](std::size_t j) {
                        xo_minor(i, j) = m_data[LAYOUT::Index(i, j)];
                    });
                });

                // top right region
                static_for<0, I>([&](std::size_t i) {
                    static_for<J + 1, COL>([&](std::size_t j) {
                        xo_minor(i, j - 1) = m_data[LAYOUT::Index(i, j)];
                    });
                });

                // bottom left region
                static_for<I + 1, ROW>([&](std::size_t i) {
                    static_for<0, J>([&](std::size_t j) {
                        xo_minor(i - 1, j) = m_data[LAYOUT::Index(i, j)];
                    });
                });

                // bottom right region
                static_for<I + 1, ROW>([&](std::size_t i) {
                    static_for<J + 1, COL>([&](std::size_t j) {
                        xo_minor(i - 1, j - 1) = m_data[LAYOUT::Index(i, j)];
                    });
                });

                return xo_minor;
            }

        // queries
        public:

            // minimal element in the entire matrix
            constexpr T MinElement() const noexcept {
                return *std::min_element(&m_data[0], &m_data[SIZE]);
            }

            // maximal element in the entire matrix
            constexpr T MaxElement() const noexcept {
                return *std::max_element(&m_data[0], &m_data[SIZE]);
            }

            // check if a unary predicate returns true for all elements in matrix
            bool AllOf(std::function<bool(const T)> xi_predicate) const noexcept {
                return std::all_of(&m_data[0], &m_data[SIZE], xi_predicate);
            }

            // check if a unary predicate returns true for at least one elements in matrix
            bool AnyOf(std::function<bool(const T)> xi_predicate) const noexcept {
                return std::any_of(&m_data[0], &m_data[SIZE], xi_predicate);
            }

            // check if a unary predicate returns true for no elements in matrix
            bool NoneOf(std::function<bool(const T)> xi_predicate) const noexcept {
                return std::none_of(&m_data[0], &m_data[SIZE], xi_predicate);
            }

            // return occurrence of a given value in a matrix
            constexpr std::size_t Count(const T& xi_value) const noexcept {
                return std::count(&m_data[0], &m_data[SIZE], xi_value);
            }

            // return MATRIX norm (L2; its 'length')
            constexpr T Norm() const noexcept {
                T xo_dot{};
                static_for<0, SIZE()>([&](std::size_t i) {
                    xo_dot += m_data[i] * m_data[i];
                });

                return std::sqrt(xo_dot);
            }

            /**
            * \brief tests whether MATRIX is normalized.
            *
            * @param {double, in}  tolerance for normalization test (default is 2 * epsilon)
            * @param {bool,   out} true if squared length indicate a normalized algebric structure.
            **/
            constexpr bool IsNormalized(const T& xi_tol = static_cast<T>(2) * FloatingPointTrait<T>::epsilon()) const noexcept {
                return (std::abs(Norm() - static_cast<T>(1)) < static_cast<T>(2) * FloatingPointTrait<T>::epsilon());
            }

            // return matrix size
            constexpr std::size_t Size() const noexcept { return SIZE; }

            // test if matrix is row major
            constexpr bool IsRowMajor() const noexcept { return m_RowMajor; }

            // test if matrix is square
            constexpr bool IsSqure() const noexcept { return (COL == ROW); }

            /**
            * \brief check if a CUBIC matrix is symmetric (around its diagonal; A = A^T)
            *
            * @param {bool, out} true if matrix is symmetric (around its diagonal)
            **/
            constexpr bool IsSymmetric() const noexcept {
                static_assert(ROW == COL, "MatrixNM must be cubic");

                bool xo_symmetric{ true };
                for (std::size_t j{}; (j < ROW) && xo_symmetric; ++j) {
                    for (std::size_t i{}; (i < ROW) && xo_symmetric; ++i) {
                        xo_symmetric = FloatingPointTrait<T>::Equals(m_data[LAYOUT::Index(i, j)], m_data[LAYOUT::Index(j, i)]);
                    }
                }
            
                return xo_symmetric;
            }
            
            /**
            * \brief check if a CUBIC matrix is skew-symmetric (around its diagonal; -A = A^T)
            *
            * @param {bool, out} true if matrix is skew-symmetric (around its diagonal)
            **/
            constexpr bool IsSkewSymmetric() const noexcept {
                static_assert(ROW == COL, "MatrixNM must be cubic");

                bool xo_symmetric{ true };
                for (std::size_t j{}; (j < ROW) && xo_symmetric; ++j) {
                    for (std::size_t i{}; (i < ROW) && xo_symmetric; ++i) {
                        xo_symmetric = FloatingPointTrait<T>::Equals(m_data[LAYOUT::Index(i, j)], -m_data[LAYOUT::Index(j, i)]);
                    }
                }
            
                return xo_symmetric;
            }
            
            /**
            * \brief check if a CUBIC matrix is upper triangular
            *
            * @param {bool, out} true if matrix is upper triangular
            **/
            constexpr bool IsUpperTriangular() const noexcept {
                static_assert(ROW == COL, "MatrixNM must be cubic");

                bool xo_triangular{ true };
                for (std::size_t i{}; (i < ROW) && xo_triangular; ++i) {
                    for (std::size_t j{}; (j < i) && xo_triangular; ++j) {
                        xo_triangular = IsZero(m_data[LAYOUT::Index(i, j)]);
                    }
                }
            
                return xo_triangular;
            }
            
            /**
            * \brief check if a CUBIC matrix is lower triangular
            *
            * @param {bool, out} true if matrix is lower triangular
            **/
            constexpr bool IsLowerTriangular() const noexcept {
                static_assert(ROW == COL, "MatrixNM must be cubic");

                bool xo_triangular{ true };
                for (std::size_t i{}; (i < ROW) && xo_triangular; ++i) {
                    for (std::size_t j{ i + 1 }; (j < ROW) && xo_triangular; ++j) {
                        xo_triangular = IsZero(m_data[LAYOUT::Index(i, j)]);
                    }
                }
            
                return xo_triangular;
            }
            
           
            /**
            * \brief check if a CUBIC matrix is diagonal
            *
            * @param {bool, out} true if matrix is diagonal
            **/
            constexpr bool IsDiagonal() const noexcept {
                static_assert(ROW == COL, "MatrixNM must be cubic");

                // diagonal
                bool xo_diagonal{ true };
                for (std::size_t i{}; (i < ROW) && xo_diagonal; ++i) {
                    xo_diagonal = !IsZero(m_data[LAYOUT::Index(i, i)]);
                }
            
                if (!xo_diagonal) {
                    return false;
                }

                // rest of matrix
                for (std::size_t i{}; (i < ROW) && xo_diagonal; ++i) {
                    for (std::size_t j{}; (j < COL) && xo_diagonal; ++j) {
                        if (i != j) {
                            xo_diagonal = IsZero(m_data[LAYOUT::Index(i, j)]);
                        }
                    }
                }
            
                return xo_diagonal;
            }
            
            /**
            * \brief check if a CUBIC matrix is permutation matrix (there is only one '1' in every column and row)
            *
            * @param {bool, out} true if matrix is permutation matrix
            **/
            constexpr bool IsPermutation() const noexcept {
                static_assert(ROW == COL, "MatrixNM must be cubic");

                const T rowTol{ static_cast<T>(1) + static_cast<T>(ROW) * FloatingPointTrait<T>::epsilon() },
                        colTol{ static_cast<T>(1) + static_cast<T>(COL) * FloatingPointTrait<T>::epsilon() };

                bool xo_permutation{ true };
                for (std::size_t i{}; (i < ROW) && xo_permutation; ++i) {                   
                    VectorN<T, COL>&& tempRow{ this->operator()(Row(i)) };
                    xo_permutation = !(Sum(tempRow) > rowTol);
            
                    VectorN<T, ROW>&& tempCol{ this->operator()(Column(i)) };
                    xo_permutation &= !(Sum(tempCol) > colTol);
                }
            
                return xo_permutation;
            }

            /**
            * \brief tests whether the matrix is per-column-normalized.
            *
            * @param {T,    in}  tolerance for normalization test (default is 2 * epsilon)
            * @param {bool, out} true if each of matrix columns is normalized, false otherwise
            **/
            constexpr bool IsNormalizedPerColumn(const T& xi_tol = T(2) * FloatingPointTrait<T>::epsilon()) const noexcept {
                bool xo_normalized{ true };

                for (std::size_t i{}; (i < COL) && xo_normalized; ++i) {
                    VectorN<T, ROW>&& column{ this->operator()(Column(i)) };
                    xo_normalized = std::move(column).IsNormalized(xi_tol);
                }

                return xo_normalized;
            }

            /**
            * \brief tests whether the matrix is per-row-normalized.
            *
            * @param {T,    in}  tolerance for normalization test (default is 2 * epsilon)
            * @param {bool, out} true if each of matrix rows is normalized, false otherwise
            **/
            constexpr bool IsNormalizedPerRow(const T& xi_tol = T(2) * FloatingPointTrait<T>::epsilon()) const noexcept {
                bool xo_normalized{ true };
                for (std::size_t i{}; (i < ROW) && xo_normalized; ++i) {
                    VectorN<T, COL>&& row{ this->operator()(Row(i)) };
                    xo_normalized = std::move(row).IsNormalized(xi_tol);
                }

                return xo_normalized;
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
                static_for<0, SIZE>([&](std::size_t i) {
                    m_data[i] = Numeric::Clamp(m_data[i], xi_boundary1,xi_boundary2);
                });
            }

            /**
            * \brief limit all of matrix values into a symmetric (about zero) boundary
            *
            * @param {BaseMatrix, in|out} matrix to limit / limited matrix
            * @param {T,          in}     boundary
            **/
            void ClampSym(const T xi_boundary) noexcept {
                static_for<0, SIZE>([&](std::size_t i) {
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
                static_for<0, SIZE>([&](std::size_t i) {
                    m_data[i] = Numeric::ClampSym(m_data[i], xi_min, xi_max);
                });
            }

        // matrix decomposition
        public:

            /**
            * \brief perform Lower-Upper decomposition for a CUBIC matrix (using "Doolittle" algorithm)
            *        lower triangulated matrix is also unit matrix.
            *
            * @param {matrixnm,  out} matrix whos lower and upper portions are the decomposition outcome
            * @param {vectorn,   out} decomposition pivot vector (row vector, i.e - VectorN<std::size_t, COL>)
            * @param {int32_t,   out} pivot sign
            **/
            void LU(MatrixNM<T, ROW, COL>& xo_lu, VectorN<std::size_t, COL>& xo_pivot, int32_t& xo_sign) noexcept {
                static_assert(ROW == COL, "LU decomposition can be performed only on a cubic matrix.");

                // housekeeping
                xo_sign = 1;
                static_for<0, COL>([&](std::size_t i) {
                    xo_pivot[i] = i;
                });

                MatrixNM<T, ROW, COL>&& self{ m_data };
                xo_lu = std::move(self);

                // Outer loop
                static_for<0, COL>([&](std::size_t j) {

                    // j-th column
                    VectorN<T, ROW> luCol{ xo_lu(Column(j)) };

                    // Apply previous transformations
                    static_for<0, COL>([&](std::size_t i) {
                        // i-th row
                        VectorN<T, COL> luRow{ xo_lu(Row(i)) };

                        // "left looking" dot product
                        T s{ Dot(luRow, luCol, 0, std::min(i, j)) };

                        luRow[j] = luCol[i] -= s;

                        xo_lu.SetRow(i, luRow);
                    });

                    // Find pivot and exchange if necessary
                    std::size_t piv{ j };
                    for (std::size_t i{ j + 1 }; i < COL; ++i) {
                        piv = (std::abs(luCol[i]) > std::abs(luCol[piv])) ? i : piv;
                    }

                    // swap around pivot
                    if (piv != j) {
                        static_for<0, COL>([&](std::size_t k) {
                            std::swap(xo_lu(piv, k), xo_lu(j, k));
                        });

                        std::swap(xo_pivot[piv], xo_pivot[j]);
                        xo_sign = -xo_sign;
                    }

                    // calculate multipliers
                    if ((j < COL) && !IsZero(xo_lu(j, j))) {
                        for (std::size_t i{ j + 1 }; i < COL; ++i) {
                            xo_lu(i, j) /= xo_lu(j, j);
                        }
                    }
                });
            }

            /**
            * \brief perform QR decomposition using "Givens rotations".
            *
            * @param {matrix, in}  matrix to be decomposed (whose number of rows is either equal or larger then the number of columns)
            * @param {matrix, out} Q matrix (orthogonal matrix with orthogonal columns, i.e. - Q*Q^T = I; COLxROW)
            * @param {matrix, out} R matrix (upper triangular matrix; COLxCOL)
            **/
            void QR(MatrixNM<T, COL, ROW>& xo_Q, MatrixNM<T, COL, COL>& xo_R) noexcept {

                /**
                * "givens rotation" (return in order of {cosine, sine, radius})
                **/
                auto GivensRotationInternal = [](const T& a, const T& b) {
                    std::array<T, 3> xo_rotation;

                    // a == 0
                    if (IsZero(a)) {
                        xo_rotation = {T{}, std::copysign(T(1), b), std::abs(b)};
                    } // b == 0
                    else if (IsZero(b)) {
                        xo_rotation = {std::copysign(T(1), a), T{}, std::abs(a)};
                    } // a > b
                    else if (std::abs(a) > std::abs(b)) {
                        T t{ b / a },
                          u{ std::copysign(T(1), a) * std::sqrt(T(1) + t * t) },
                          c{ T(1) / u};

                        xo_rotation = { c, t * c, a * u};
                    }
                    else {
                        T t{ a / b },
                          u{std::copysign(T(1), b) * std::sqrt(T(1) + t * t)},
                          s{T(1) / u};

                        xo_rotation = {t * s, s, b * u};
                    }

                    return std::move(xo_rotation);
                };

                MatrixNM<T, ROW, COL>&& self{ m_data };
                MatrixNM<T, ROW, COL> R(std::move(self));

                MatrixNM<T, ROW, ROW> Q(T{});
                static_for<0, ROW>([&](std::size_t i) {
                    Q(i, i) = static_cast<T>(1.0);
                });

                static_for<0, COL>([&](std::size_t j) {
                    for (std::size_t i{ ROW - 1 }; i >= j + 1; --i) {
                        std::array<T, 3> CSR = GivensRotationInternal(R(i - 1, j), R(i, j));

                        // R' = G * R
                        static_for<0, COL>([&](std::size_t x) {
                            T temp1{ R(i - 1, x) },
                              temp2{ R(i, x) };
                            R(i - 1, x) =  temp1 * CSR[0] + temp2 * CSR[1];
                            R(i,     x) = -temp1 * CSR[1] + temp2 * CSR[0];
                        });
                        R(i - 1, j) = CSR[2];
                        R(i,     j) = T{};

                        // Q' = Q * G^T
                        static_for<0, ROW>([&](std::size_t x) {
                            T temp1{ Q(x, i - 1) },
                              temp2{ Q(x, i) };
                            Q(x, i - 1) =  temp1 * CSR[0] + temp2 * CSR[1];
                            Q(x, i)     = -temp1 * CSR[1] + temp2 * CSR[0];
                        });
                    }
                });

                // adjust Q
                if constexpr (ROW != COL) {
                    Q = Q.Transpose();
                }
                xo_Q = Q.GetRegion<0, COL, 0, ROW>();

                // adjust R
                xo_R = R.GetRegion<0, COL, 0, COL>();
            }

            /**
            * \brief perform singular value decomposition on a given matrix, i.e. - given A, return A = U * W * V^T
            *        where U & V columns are orthonormal (U^T*U = U*U^T = V^T*V = V*V^T = I).
            *
            *        important remarks:
            *        1) to avoid arithmetic complexity which is not needed for most SVD related operations:
            *           a) instead of U, U*W is returned.
            *              i.e - for a 3x3 matrix: 
            *              {U(0, 0) * sqrt(W(0, 1)), U(0, 1) * sqrt(W(1, 1)), U(0, 2) * sqrtW(2, 1));
            *               U(1, 0) * sqrt(W(0, 1)), U(1, 1) * sqrt(W(1, 1)), U(1, 2) * sqrtW(2, 1));
            *               U(2, 0) * sqrt(W(0, 1)), U(2, 1) * sqrt(W(1, 1)), U(2, 2) * sqrtW(2, 1))}
            *           b) W holds the SQUARE of singular values (i.e. - "eigenvalues").
            *           C) V is not transposed.
            *        2) SVD values are numerically correct, but their SIGN MIGHT BE WRONG. For further reading on this subject, see:
            *           SANDIA REPORT, SAND2007-6422, "Resolving the Sign Ambiguity in the Singular Value Decomposition", R. Bra, E. Acar, T. Kolda.
            *
            * @param {matrix, in}  matrix to be decomposed
            * @param {matrix, out} U*W [ROWxCOL]
            * @param {vector, out} VectorN holding SQUARE of singular values [COLx1]
            * @param {matrix, out} V (not V transposed!) [COLxCOL]
            * @param {bool,   out} true if SVD rutine converged, false otherwise
            **/
            bool SVD(MatrixNM<T, ROW, COL>& xo_UW, VectorN<T, COL>& xo_W2, MatrixNM<T, COL, COL>& xo_V) noexcept {
                static_assert(ROW >= COL, "MatrixNM::SVD can only be perofrmed on a matrix whose number of rows is equal or larger then the number of columns.");

                const std::size_t sweepLimit{ (COL < 120) ? 30 : (COL / 4) };
                std::size_t estimatedColumnRank{ COL },
                            rotationCount{ ROW },
                            sweepCount{};

                const T eps{ FloatingPointTrait<T>::epsilon() },
                        e2{ static_cast<T>(10.0) * ROW * eps * eps },
                        tol{ static_cast<T>(0.1) * eps },               // algorithm convergence tolerance
                        tolSqr{ tol * tol };

                // fill A
                MatrixNM<T, ROW + COL, COL> A(T{});
                static_for<0, ROW>([&](std::size_t i) {
                    static_for<0, COL>([&](std::size_t j) {
                        A(i, j) = m_data[LAYOUT::Index(i, j)];
                    });
                });

                static_for<0, COL>([&](std::size_t i) {
                    A(ROW + i, i) = static_cast<T>(1.0);
                });

                while ((rotationCount != 0) && (sweepCount++ <= sweepLimit)) {
                    rotationCount = estimatedColumnRank * (estimatedColumnRank - 1) / 2;

                    // Jacobi sweep
                    for (std::size_t j{}; j < estimatedColumnRank - 1; ++j) {
                        for (std::size_t k{ j + 1 }; k < estimatedColumnRank; ++k) {
                            T p{}, q{}, r{};

                            static_for<0, ROW>([&](std::size_t i) {
                                const T x0{ A(i, j) },
                                        y0{ A(i, k) };
                                p += x0 * y0;
                                q += x0 * x0;
                                r += y0 * y0;
                            });
                            xo_W2[j] = q;
                            xo_W2[k] = r;

                            // are columns ordered?
                            if (q >= r) {
                                if ((q <= e2 * xo_W2[0]) ||           // are columns 'norm' so small that rotation makes no sense?
                                    (std::abs(p) <= tol * q)) {       // is the inner product is small with respect to the larger of the columns?
                                                                      // then rotation angle is very small and we can skip it
                                    --rotationCount;
                                } // givens rotation
                                else {
                                    const T qInv{ static_cast<T>(1.0) / q };
                                    p *= qInv;
                                    r = static_cast<T>(1.0) - r * qInv;
                                    const T vt{ std::sqrt(static_cast<T>(4.0) * p * p + r * r) },
                                            c0{ std::sqrt(static_cast<T>(0.5) * (static_cast<T>(1.0) + r / vt)) },
                                            s0{ p / (vt * c0) };

                                    static_for<0, ROW + COL>([&](std::size_t i) {
                                        const T d1{ A(i, j) },
                                                d2{ A(i, k) };
                                        A(i, j) =  d1 * c0 + d2 * s0;
                                        A(i, k) = -d1 * s0 + d2 * c0;
                                    });
                                }
                            } // givens rotation
                            else {
                                const T rInv{ static_cast<T>(1.0) / r };
                                p *= rInv;
                                q = q * rInv - static_cast<T>(1.0);

                                const T vt{ std::sqrt(static_cast<T>(4.0) * p * p + q * q) };

                                T s0{ std::sqrt(static_cast<T>(0.5) * (static_cast<T>(1.0) - q / vt)) };
                                if (IsNegative(p)) {
                                    s0 = -s0;
                                }

                                const T c0{ p / (vt * s0) };

                                static_for<0, ROW + COL>([&](std::size_t i) {
                                    const T d1{ A(i, j) },
                                            d2{ A(i, k) };
                                    A(i, j) =  d1 * c0 + d2 * s0;
                                    A(i, k) = -d1 * s0 + d2 * c0;
                                });
                            }
                        }
                    }

                    const T bound{ xo_W2[0] * tol + tolSqr };
                    while ((estimatedColumnRank > 2) && (xo_W2[estimatedColumnRank - 1] <= bound)) {
                        estimatedColumnRank--;
                    }
                }

                // reached maximal number of sweeps
                if (sweepCount > sweepLimit) {
                    return false;
                }

                // extract V
                xo_V = A.GetRegion<ROW, ROW + COL, 0, COL>();

                // extract U*W
                xo_UW = A.GetRegion<0, ROW, 0, COL>();

                return true;
            }

            /**
            * \brief given matrix A, constructs a lower triangular matrix L such that L*L' = A.
            *        'A' must be symmetric positive definite.
            *        This decomposition is called Cholesky decomposition.
            *        (it's roughly TWICE as efficient as the LU decomposition)
            *
            * @param {matrix, in}  'A' (must be symmetric)
            * @param {matrix, out} lower decomposition
            **/
            void Cholesky(MatrixNM<T, ROW, COL>& xo_lower) noexcept {
                static_assert(ROW == COL, "");

                xo_lower = T{};

                static_for<0, COL>([&](std::size_t j) {
                    T d{};

                    for (std::size_t k{}; k < j; ++k) {
                        T s{};

                        for (std::size_t i{}; i < k; ++i) {
                            s += xo_lower(k, i) * xo_lower(j, i);
                        }

                        xo_lower(j, k) = s = (this->operator()(j, k) - s) / xo_lower(k, k);
                        d += s * s;
                    }

                    d = this->operator()(j, j) - d;

                    xo_lower(j, j) = IsPositive(d) ? (std::sqrt(d)) : (T{});
                });
            }

        // matrix properties
        public:

            /**
            * \brief return the determinant of a given matrix
            *        (notice, performs LU decomposition for matrix above 4x4)
            *
            * @param {matrix,  in} matrix
            * @param {T,      out} matrix determinant
            **/
            constexpr T Determinant() const noexcept {
                static_assert(ROW == COL, "Determinant operations can only be done on cubic matrix.");

                T xo_det{};

                if constexpr (ROW == 2) {
                    xo_det = this->operator()(0, 0) * this->operator()(1, 1)  - this->operator()(0, 1) * this->operator()(1, 0);
                }
                else if constexpr (ROW == 3) {
                    xo_det = this->operator()(0, 0) * (this->operator()(1, 1) * this->operator()(2, 2) - this->operator()(2, 1) * this->operator()(1, 2)) -
                             this->operator()(0, 1) * (this->operator()(1, 0) * this->operator()(2, 2) - this->operator()(2, 0) * this->operator()(1, 2)) +
                             this->operator()(0, 2) * (this->operator()(1, 0) * this->operator()(2, 1) - this->operator()(2, 0) * this->operator()(1, 1));
                }
                else if constexpr (ROW == 4) {

                    // determinant using cofactors (along first row)
                    MatrixNM<T, ROW, COL>&& self{ m_data };
                    xo_det = std::move(self(0, 0)) * std::move(self).Minor<0, 0>().Determinant() - std::move(self(0, 1)) * std::move(self).Minor<0, 1>().Determinant() +
                             std::move(self(0, 2)) * std::move(self).Minor<0, 2>().Determinant() - std::move(self(0, 3)) * std::move(self).Minor<0, 3>().Determinant();
                }
                else {

                    MatrixNM<T, ROW, COL>&& self{ m_data };

                    // LU decomposition
                    MatrixNM<T, ROW, COL> lowerUpper;
                    VectorN<std::size_t, COL> piv;
                    int32_t _sign;
                    self.LU(lowerUpper, piv, _sign);

                    // determinant calculation
                    xo_det = static_cast<T>(_sign);

                    static_for<0, COL>([&](std::size_t i) {
                        xo_det *= lowerUpper(i, i);
                    });                 
                }

                return xo_det;
            }

            // return a given cofactor of matrix
            template<std::size_t I, std::size_t J> constexpr T Cofactor() const noexcept {
                const T coeff{ std::pow(static_cast<T>(-1.0), static_cast<T>(I + J)) };
                MatrixNM<T, ROW, COL>&& self{ m_data };
                return (coeff * std::move(self).Minor<I, J>().Determinant());
            }

        // linear system solvers
        public:

            /**
            * \brief solve linear system A*x=b (A must be CUBIC).
            *        (notice, performs LU decomposition internally)
            *
            * @param {matrix, in}  A (CUBIC)
            * @param {vector, in}  b (column vector)
            * @param {vector, out} x (column vector)
            **/
            void SolveSquare(const VectorN<T, ROW>& xi_b, VectorN<T, ROW>& xo_x) noexcept {
                static_assert(ROW == COL, "MatrixNM::SolveSquare can only operate on cubic matrix (i.e. - ROW == COL).");

                // LU decomposition
                MatrixNM<T, ROW, COL>&& self{ m_data };
                MatrixNM<T, ROW, COL> lowerUpper;
                VectorN<std::size_t, COL> piv;
                int32_t _sign;
                self.LU(lowerUpper, piv, _sign);
                
                // x is the permuted copy of xi_b as xi_pivot
                static_for<0, COL>([&](std::size_t i) {
                    xo_x[i] = xi_b[piv[i]];
                });

                // Solve L*Y = B(pivoted)
                static_for<0, COL>([&](std::size_t k) {
                    for (std::size_t i{ k + 1 }; i < COL; ++i) {
                        xo_x[i] -= xo_x[k] * lowerUpper(i, k);
                    }
                });

                // Solve U*X = Y
                for (int64_t k{ COL - 1 }; k >= 0; k--) {
                    xo_x[static_cast<std::size_t>(k)] /= lowerUpper(static_cast<std::size_t>(k), static_cast<std::size_t>(k));

                    for (std::size_t i{}; i < static_cast<std::size_t>(k); ++i) {
                        xo_x[i] -= xo_x[static_cast<std::size_t>(k)] * lowerUpper(i, static_cast<std::size_t>(k));
                    }
                }
            }

            /**
            * \brief solve linear system A*x=b (using internally both QR & LU decomposition)
            *
            *        Notice that QR decomposition is used and not:
            *        > pseudo-inverse - to avoid increasing the output matrix condition number (happens when multiplying the matrix by its transpose),
            *        > SVD - high running time complexity.
            *
            * @param {matrix, in}  A (ROWxCOL, ROW >= COL)
            * @param {vector, in}  B (column vector, ROWx1)
            * @param {vector, out} X (column vector, COL, 1)
            **/
            void SolveRectangualr(const VectorN<T, ROW>& xi_B, VectorN<T, COL>& xo_X) noexcept {
                static_assert(ROW >= COL, "MatrixNM::SolveRectangualr can only be perofrmed on a matrix whose number of rows is equal or larger then the number of columns.");

                // QR decomposition
                MatrixNM<T, COL, ROW> Q;
                MatrixNM<T, COL, COL> R;
                this->QR(Q, R);
                
                // C = Q * B
                VectorN<T, COL> C{ Q * xi_B };

                // R*x = C
                MatrixNM<T, COL, COL> Rlu;
                VectorN<std::size_t, COL> piv;
                int32_t sign;
                R.LU(Rlu, piv, sign);

                Rlu.SolveSquare(C, xo_X);
            }

            /**
            * \brief solve linear system A*x=b using Cholesky decomposition.
            *        'A' must be CUBIC and symmetric positive definite. 
            *        (only the CUBIC requirement is enforced)
            *
            * @param {matrix, in}  A (CUBIC)
            * @param {vector, in}  b (column vector)
            * @param {vector, out} x (column vector)
            **/
            void SolveCholesky(const VectorN<T, ROW>& xi_b, VectorN<T, ROW>& xo_x) noexcept {
                static_assert(ROW == COL, "MatrixNM::SolveCholesky can only operate on cubic matrix (i.e. - ROW == COL).");
                MatrixNM<T, ROW, COL> L;
                this->Cholesky(L);

                xo_x = xi_b;

                // Solve L*y = b;
                static_for<0, COL>([&](std::size_t k) {
                    for (std::size_t i{}; i < k; ++i) {
                        xo_x[k] -= xo_x[i] * L(k, i);
                    }

                    xo_x[k] /= L(k, k);
                });

                // Solve L'*X = Y;
                for (int64_t k{ COL - 1 }; k >= 0; --k) {
                    const std::size_t ks{ static_cast<std::size_t>(k) };

                    for (std::size_t i{ ks + 1 }; i < COL; ++i) {
                        xo_x[ks] -= xo_x[i] * L(i, ks);
                    }

                    xo_x[ks] /= L(ks, ks);
                }
            }
    };

    /**
    * numerical operator overload
    **/

    // --- unary minus ---
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> operator - (const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat) {
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat(xi_mat);
        xo_mat *= static_cast<T>(-1);
        return xo_mat;
    }

    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator - (MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat) {
        static_for<0, ROW * COL>([&](std::size_t i) {
            xi_mat[i] = -xi_mat[i];
        });
        return xi_mat;
    }

    // binary operations without left hand side scalar
#define M_BINARY_OP_NO_LHS_SCALAR(OP, AOP)                                                                                                                          \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                       \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> operator OP (const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat, const T xi_value) {                                    \
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat(xi_mat);                                                                                                               \
        xo_mat AOP xi_value;                                                                                                                                        \
        return xo_mat;                                                                                                                                              \
    }                                                                                                                                                               \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                       \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator OP (MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat, const T xi_value) {                                        \
        xi_mat AOP xi_value;                                                                                                                                        \
        return xi_mat;                                                                                                                                              \
    }                                                                                                                                                               \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                       \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> operator OP (const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat1, const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat2) {       \
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat(xi_mat1);                                                                                                              \
        xo_mat AOP xi_mat2;                                                                                                                                         \
        return xo_mat;                                                                                                                                              \
    }                                                                                                                                                               \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                       \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator OP (MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat1, MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat2) {                \
        xi_mat1 AOP xi_mat2;                                                                                                                                        \
        return xi_mat1;                                                                                                                                             \
    }

    M_BINARY_OP_NO_LHS_SCALAR(-, -=);
    M_BINARY_OP_NO_LHS_SCALAR(/, /=);
    M_BINARY_OP_NO_LHS_SCALAR(&, &=);
    M_BINARY_OP_NO_LHS_SCALAR(|, |=);
    M_BINARY_OP_NO_LHS_SCALAR(^, ^=);
    M_BINARY_OP_NO_LHS_SCALAR(>>, >>=);
    M_BINARY_OP_NO_LHS_SCALAR(<<, <<=);

#undef M_BINARY_OP_NO_LHS_SCALAR

    // binary operations with left hand side scalar
#define M_BINARY_OP_WITH_LHS_SCALAR(OP, AOP)                                                                                                                          \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                         \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> operator OP (const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat, const T xi_value) {                                      \
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat(xi_mat);                                                                                                                 \
        xo_mat AOP xi_value;                                                                                                                                          \
        return xo_mat;                                                                                                                                                \
    }                                                                                                                                                                 \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                         \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator OP (MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat, const T xi_value) {                                          \
        xi_mat AOP xi_value;                                                                                                                                          \
        return xi_mat;                                                                                                                                                \
    }                                                                                                                                                                 \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                         \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> operator OP (const T xi_value, const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat) {                                      \
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat(xi_mat);                                                                                                                 \
        xo_mat AOP xi_value;                                                                                                                                          \
        return xo_mat;                                                                                                                                                \
    }                                                                                                                                                                 \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                         \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator OP (const T xi_value, MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat) {                                          \
        xi_mat AOP xi_value;                                                                                                                                          \
        return xi_mat;                                                                                                                                                \
    }                                                                                                                                                                 \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                         \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> operator OP (const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat1, const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat2) {         \
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat(xi_mat1);                                                                                                                \
        xo_mat AOP xi_mat2;                                                                                                                                           \
        return xo_mat;                                                                                                                                                \
    }                                                                                                                                                                 \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                                         \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator OP (MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat1, MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat2) {                  \
        xi_mat1 AOP xi_mat2;                                                                                                                                          \
        return xi_mat1;                                                                                                                                               \
    }

    M_BINARY_OP_WITH_LHS_SCALAR(+, +=);
    M_BINARY_OP_WITH_LHS_SCALAR(*, *=);

#undef M_BINARY_OP_WITH_LHS_SCALAR

    // matrix * vector (i.e - right multiply a matrix by a vector)
    // (row X 1) = (row X col) * (1 X COL)
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline VectorN<T, ROW> operator * (const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat, const VectorN<T, COL>& xi_vec) {
        VectorN<T, ROW> xo_vec(T{});

        static_for<0, ROW>([&](std::size_t i) {
            T _sum{};

            static_for<0, COL>([&](std::size_t j) {
                _sum += xi_mat(i, j) * xi_vec[j];
            });

            xo_vec[i] = _sum;
        });

        return xo_vec;
    }

    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    inline VectorN<T, ROW>& operator * (MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat, const VectorN<T, COL>& xi_vec) {
        VectorN<T, ROW> xo_vec(T{});

        static_for<0, ROW>([&](std::size_t i) {
            T _sum{};

            static_for<0, COL>([&](std::size_t j) {
                _sum += xi_mat(i, j) * xi_vec[j];
            });

            xo_vec[i] = _sum;
        });

        return xo_vec;
    }

    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    inline VectorN<T, ROW>& operator * (const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat, VectorN<T, COL>&& xi_vec) {
        VectorN<T, ROW> xo_vec(T{});

        static_for<0, ROW>([&](std::size_t i) {
            T _sum{};

            static_for<0, COL>([&](std::size_t j) {
                _sum += xi_mat(i, j) * xi_vec[j];
            });

            xo_vec[i] = _sum;
        });

        return xo_vec;
    }

    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    inline VectorN<T, ROW>& operator * (MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat, VectorN<T, COL>&& xi_vec) {
        VectorN<T, ROW> xo_vec(T{});

        static_for<0, ROW>([&](std::size_t i) {
            T _sum{};

            static_for<0, COL>([&](std::size_t j) {
                _sum += xi_mat(i, j) * xi_vec[j];
            });

            xo_vec[i] = _sum;
        });

        return xo_vec;
    }

    // vector * matrix (i.e. - left multiply a matrix by a vector)
    // (1 X col) = (1 X row) * (row X col)
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline VectorN<T, COL> operator * (const VectorN<T, ROW>& xi_vec, const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat) {
        VectorN<T, COL> xo_vec(T{});

        static_for<0, ROW>([&](std::size_t i) {
            T _sum{};

            static_for<0, COL>([&](std::size_t j) {
                _sum += xi_mat(i, j) * xi_vec[j];
            });

            xo_vec[i] = _sum;
        });

        return xo_vec;
    }

    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline VectorN<T, COL>& operator * (VectorN<T, ROW>&& xi_vec, const MatrixNM<T, ROW, COL, LAYOUT>& xi_mat) {
        VectorN<T, COL> xo_vec(T{});

        static_for<0, ROW>([&](std::size_t i) {
            T _sum{};

            static_for<0, COL>([&](std::size_t j) {
                _sum += xi_mat(i, j) * xi_vec[j];
            });

            xo_vec[i] = _sum;
        });

        return xo_vec;
    }

    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline VectorN<T, COL>& operator * (const VectorN<T, ROW>& xi_vec, MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat) {
        VectorN<T, COL> xo_vec(T{});

        static_for<0, ROW>([&](std::size_t i) {
            T _sum{};

            static_for<0, COL>([&](std::size_t j) {
                _sum += xi_mat(i, j) * xi_vec[j];
            });

            xo_vec[i] = _sum;
        });

        return xo_vec;
    }

    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline VectorN<T, COL>& operator * (VectorN<T, ROW>&& xi_vec, MatrixNM<T, ROW, COL, LAYOUT>&& xi_mat) {
        VectorN<T, COL> xo_vec(T{});

        static_for<0, ROW>([&](std::size_t i) {
            T _sum{};

            static_for<0, COL>([&](std::size_t j) {
                _sum += xi_mat(i, j) * xi_vec[j];
            });

            xo_vec[i] = _sum;
        });

        return xo_vec;
    }

    // matrix * matrix
    // (row X col) = (row X dim) * (dim X col)
    template<typename T, std::size_t ROW, std::size_t COL, std::size_t DIM, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> operator * (const MatrixNM<T, ROW, DIM>& xi_lhs, const MatrixNM<T, DIM, COL>& xi_rhs) {
        static_assert(xi_lhs.IsRowMajor() == xi_rhs.IsRowMajor(), "MatrixNM * MatrixNM operations is only allowed when both matrix layout is identical.");
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat;

        static_for<0, ROW>([&](std::size_t i) {
            static_for<0, COL>([&](std::size_t j) {
                T _sum{};

                static_for<0, DIM>([&](std::size_t d) {
                    _sum += xi_lhs(i, d) * xi_rhs(d, j);
                });

                xo_mat(i, j) = _sum;
            });
        });

        return xo_mat;
    }

    template<typename T, std::size_t ROW, std::size_t COL, std::size_t DIM, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator * (MatrixNM<T, ROW, DIM>&& xi_lhs, const MatrixNM<T, DIM, COL>& xi_rhs) {
        static_assert(xi_lhs.IsRowMajor() == xi_rhs.IsRowMajor(), "MatrixNM * MatrixNM operations is only allowed when both matrix layout is identical.");
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat;

        static_for<0, ROW>([&](std::size_t i) {
            static_for<0, COL>([&](std::size_t j) {
                T _sum{};

                static_for<0, DIM>([&](std::size_t d) {
                    _sum += xi_lhs(i, d) * xi_rhs(d, j);
                });

                xo_mat(i, j) = _sum;
            });
        });

        return xo_mat;
    }

    template<typename T, std::size_t ROW, std::size_t COL, std::size_t DIM, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator * (const MatrixNM<T, ROW, DIM>& xi_lhs, MatrixNM<T, DIM, COL>&& xi_rhs) {
        static_assert(xi_lhs.IsRowMajor() == xi_rhs.IsRowMajor(), "MatrixNM * MatrixNM operations is only allowed when both matrix layout is identical.");
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat;

        static_for<0, ROW>([&](std::size_t i) {
            static_for<0, COL>([&](std::size_t j) {
                T _sum{};

                static_for<0, DIM>([&](std::size_t d) {
                    _sum += xi_lhs(i, d) * xi_rhs(d, j);
                });

                xo_mat(i, j) = _sum;
            });
        });

        return xo_mat;
    }

    template<typename T, std::size_t ROW, std::size_t COL, std::size_t DIM, class LAYOUT = RowMajor<ROW, COL>>
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& operator * (MatrixNM<T, ROW, DIM>&& xi_lhs, MatrixNM<T, DIM, COL>&& xi_rhs) {
        static_assert(xi_lhs.IsRowMajor() == xi_rhs.IsRowMajor(), "MatrixNM * MatrixNM operations is only allowed when both matrix layout is identical.");
        MatrixNM<T, ROW, COL, LAYOUT> xo_mat;

        static_for<0, ROW>([&](std::size_t i) {
            static_for<0, COL>([&](std::size_t j) {
                T _sum{};

                static_for<0, DIM>([&](std::size_t d) {
                    _sum += xi_lhs(i, d) * xi_rhs(d, j);
                });

                xo_mat(i, j) = _sum;
            });
        });

        return xo_mat;
    }

    /**
    * relational operator overload
    **/
#define M_OPERATOR(OP)                                                                                                                   \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                            \
    constexpr inline bool operator OP (const MatrixNM<T, ROW, COL, LAYOUT>& xi_lhs, const MatrixNM<T, ROW, COL, LAYOUT>& xi_rhs) {       \
        bool xo_rational{ false };                                                                                                       \
        for (std::size_t i{}; (i < ROW * COL) && !xo_rational; ++i) {                                                                    \
            xo_rational = xi_lhs[i] OP xi_rhs[i];                                                                                        \
        }                                                                                                                                \
        return xo_rational;                                                                                                              \
    }                                                                                                                                    \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                            \
    constexpr inline bool operator OP (const T xi_lhs, const MatrixNM<T, ROW, COL, LAYOUT>& xi_rhs) {                                    \
        bool xo_rational{ false };                                                                                                       \
        for (std::size_t i{}; (i < ROW * COL) && !xo_rational; ++i) {                                                                    \
            xo_rational = xi_lhs OP xi_rhs[i];                                                                                           \
        }                                                                                                                                \
        return xo_rational;                                                                                                              \
    }                                                                                                                                    \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                            \
constexpr inline bool operator OP (const MatrixNM<T, ROW, COL, LAYOUT>& xi_lhs, const T& xi_rhs) {                                       \
        bool xo_rational{ false };                                                                                                       \
        for (std::size_t i{}; (i < ROW * COL) && !xo_rational; ++i) {                                                                    \
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
#define M_UNARY_FUNCTION(NAME, STL_NAME)                                                                   \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>              \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> NAME(const MatrixNM<T, ROW, COL, LAYOUT>& xi_value) {   \
        MatrixNM<T, ROW, COL, LAYOUT> xo_return(xi_value);                                                 \
        static_for<0, ROW * COL>([&](std::size_t i) {                                                      \
            xo_return[i] = STL_NAME(xo_return[i]);                                                         \
        });                                                                                                \
        return xo_return;                                                                                  \
    }                                                                                                      \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>              \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& NAME(MatrixNM<T, ROW, COL, LAYOUT>&& xi_value) {       \
        static_for<0, ROW * COL>([&](std::size_t i) {                                                      \
            xi_value[i] = STL_NAME(xi_value[i]);                                                           \
        });                                                                                                \
        return xi_value;                                                                                   \
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
#define M_BINARY_FUNCTION(NAME, STL_NAME)                                                                                                                  \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                              \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> NAME(const MatrixNM<T, ROW, COL, LAYOUT>& xi_lhs, const T xi_rhs) {                                     \
        MatrixNM<T, ROW, COL, LAYOUT> xo_return();                                                                                                         \
        static_for<0, ROW * COL>([&](std::size_t i) {                                                                                                      \
            xo_return[i] = STL_NAME(xi_lhs[i], xi_rhs);                                                                                                    \
        });                                                                                                                                                \
        return xo_return;                                                                                                                                  \
    }                                                                                                                                                      \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                              \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT>& NAME(MatrixNM<T, ROW, COL, LAYOUT>&& xi_lhs, const T xi_rhs) {                                         \
        static_for<0, ROW * COL>([&](std::size_t i) {                                                                                                      \
            xi_lhs[i] = STL_NAME(xi_lhs[i], xi_rhs);                                                                                                       \
        });                                                                                                                                                \
        return xi_lhs;                                                                                                                                     \
    }                                                                                                                                                      \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                              \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> NAME(const MatrixNM<T, ROW, COL, LAYOUT>& xi_lhs, const MatrixNM<T, ROW, COL, LAYOUT>& xi_rhs) {        \
        MatrixNM<T, ROW, COL, LAYOUT> xo_return();                                                                                                         \
        static_for<0, ROW * COL>([&](std::size_t i) {                                                                                                      \
            xo_return[i] = STL_NAME(xi_lhs[i], xi_rhs[i]);                                                                                                 \
        });                                                                                                                                                \
        return xo_return;                                                                                                                                  \
    }                                                                                                                                                      \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                              \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> NAME(MatrixNM<T, ROW, COL, LAYOUT>&& xi_lhs, const MatrixNM<T, ROW, COL, LAYOUT>& xi_rhs) {             \
        static_for<0, ROW * COL>([&](std::size_t i) {                                                                                                      \
            xi_lhs[i] = STL_NAME(xi_lhs[i], xi_rhs[i]);                                                                                                    \
        });                                                                                                                                                \
        return xi_lhs;                                                                                                                                     \
    }                                                                                                                                                      \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                              \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> NAME(const MatrixNM<T, ROW, COL, LAYOUT>& xi_lhs, MatrixNM<T, ROW, COL, LAYOUT>&& xi_rhs) {             \
        static_for<0, ROW * COL>([&](std::size_t i) {                                                                                                      \
            xi_rhs[i] = STL_NAME(xi_lhs[i], xi_rhs[i]);                                                                                                    \
        });                                                                                                                                                \
        return xi_rhs;                                                                                                                                     \
    }                                                                                                                                                      \
    template<typename T, std::size_t ROW, std::size_t COL, class LAYOUT = RowMajor<ROW, COL>>                                                              \
    constexpr inline MatrixNM<T, ROW, COL, LAYOUT> NAME(MatrixNM<T, ROW, COL, LAYOUT>&& xi_lhs, MatrixNM<T, ROW, COL, LAYOUT>&& xi_rhs) {                  \
        static_for<0, ROW * COL>([&](std::size_t i) {                                                                                                      \
            xi_lhs[i] = STL_NAME(xi_lhs[i], xi_rhs[i]);                                                                                                    \
        });                                                                                                                                                \
        return xi_lhs;                                                                                                                                     \
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
    template<typename>                                 struct is_MatrixNM                    : public std::false_type {};
    template<typename T, std::size_t N, std::size_t M> struct is_MatrixNM<MatrixNM<T, N, M>> : public std::true_type  {};
    template<typename U> constexpr bool isMatrixNM(const U&) { return is_MatrixNM<U>::value; }
};
