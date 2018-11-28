/**
* A collection of 2D/3D related operations
* (some functions located here are specialized MatrixNM functions for particular 2D/3D situations)
*
* Dan Israel Malta
**/
#pragma once
#include "FloatingPointTraits.h"
#include "Constants.h"
#include "static_for.h"
#include "Common.h"
#include "VectorN.h"
#include "MatrixNM.h"
#include "Quaternion.h"
#include <optional>

namespace Numeric {

	// --------------------------
	// --- general operations ---
	// --------------------------

	/**
	* return the cross product between two vectors
	**/
	template<typename T> constexpr inline VectorN<T, 3> Cross(const VectorN<T, 3>& xi_a, const VectorN<T, 3>& xi_b) noexcept {
		return VectorN<T, 3>{ xi_a[1] * xi_b[2] - xi_b[1] * xi_a[2],
							  xi_b[0] * xi_a[2] - xi_a[0] * xi_b[2],
							  xi_a[0] * xi_b[1] - xi_b[0] * xi_a[1] };
	}

	// ------------------------------------------
	// --- special operations with 2x2 matrix ---
	// ------------------------------------------

	/**
    * \brief return the angle which rotates a 2x2 matrix into diagonal form,
    *        i.e - elements (1, 0) = (0,1) = 0, and elements (0, 0) > (1, 1) are the eigenvalues
    *
    * @param {matrix, in}  matrix
    * @param {T,      out} angle [rad]
    **/
    template<typename T, class LAYOUT> constexpr inline T Diagonalizer(const MatrixNM<T, 2, 2, LAYOUT>& xi_mat) noexcept {
        const T diff{ xi_mat(1, 1) - xi_mat(0, 0 ) };
        return std::atan2(diff + std::sqrt(diff * diff + static_cast<T>(4) * xi_mat(0, 1) * xi_mat(1, 0)), 
						  static_cast<T>(2) * xi_mat(0, 1));
    }

	/**
    * \brief return matrix eigenvalues
    *
    * @param {matrix, in}  matrix
    * @param {T,      out} eigenvalue
    * @param {T,      out} eigenvalue
    **/
    template<typename T, class LAYOUT> void EigenValues(const MatrixNM<T, 2, 2, LAYOUT>& xi_mat, T& xo_eigen1, T& xo_eigen2) noexcept {
        const T diff{ xi_mat(0, 0 ) - xi_mat(1, 1) },
			    center{ xi_mat(0, 0 ) + xi_mat(1, 1) },
			    delta{ std::sqrt(diff * diff + T(4) * xi_mat(1, 0) * xi_mat(0, 1)) };

        xo_eigen1 = static_cast<T>(0.5) * (center + delta);
        xo_eigen2 = static_cast<T>(0.5) * (center - delta);
    }

	/**
    * \brief return matrix eigenvectors
    *
    * @param {matrix, in}  matrix
    * @param {vector, out} eigenvector #1
    * @param {vector, out} eigenvector #2
    **/
    template<typename T, class LAYOUT> void EigenVector(const MatrixNM<T, 2, 2, LAYOUT>& xi_mat, VectorN<T, 2>& xo_eigen1, VectorN<T, 2>& xo_eigen2) noexcept {
        T eigenvalue1, eigenvalue2;

        EigenValues(xi_mat, eigenvalue1, eigenvalue2);

        xo_eigen1 = { -xi_mat(1, 0), xi_mat(0, 0 ) - eigenvalue1 };
        xo_eigen2 = { -xi_mat(1, 0), xi_mat(0, 0 ) - eigenvalue2 };

        xo_eigen1.Normalize();
        xo_eigen2.Normalize();
    }

	/**
    * \brief test if 2x2 matrix is orthonormal
    *
    * @param {matrix, in}  matrix
    * @param {bool,   out} true if a matrix is orthonormal
    **/
     template<typename T, class LAYOUT> constexpr bool IsOrthonormal(const MatrixNM<T, 2, 2, LAYOUT>& xi_mat) noexcept {
		const VectorN<T, 2> x{ xi_mat(Column(0)) };

        if (!FloatingPointTrait<T>::Equals(Dot(x, x), T(1))) {
            return false;
        }

		const VectorN<T, 2> y{ xi_mat(Column(1)) };
        if (!FloatingPointTrait<T>::Equals(Dot(y, y), T(1))) {
            return false;
        }

        if (!IsZero(Dot(x, y))) {
            return false;
        }

        return true;
    }

	/**
    * \brief perform SVD decomposition of a SYMMETRIC matrix, i.e.:
    *           A      =      U    *     W     *   U^T
    *        [ A  B ]  =  [ c  -s ] [ r1   0  ] [  c  s ]
    *        [ B  C ]     [ s   c ] [  0   r2 ] [ -s  c ]
    *
    *        notice that as with any SVD decomposition which does not incorporate sign ambiguity logic,
    *        this method will output the correct values, but not necessarily the correct signs.
    *
    * @param {matrix, in}  SYMMETRIC matrix to be decomposed
    * @param {matrix, out} U
    * @param {vector, out} W (as a vector)
    **/
    template<typename T, class LAYOUT> void SVDsymmetric2x2(const MatrixNM<T, 2, 2, LAYOUT>& xi_mat, MatrixNM<T, 2, 2, LAYOUT>& xo_U, VectorN<T, 2>& xo_W) noexcept {
        assert(FloatingPointTrait<T>::Equals(xi_mat(0, 1), xi_mat(1, 0)));

        const T A{ xi_mat(0, 0 ) },
				B{ xi_mat(0, 1) },
				C{ xi_mat(1, 1) },
				traceSum{ A + C },
				traceDiff{ A - C },
				rt{ std::sqrt(traceDiff * traceDiff + T(4) * B * B) };
        T r1, r2, c, s, t;

        // eigenvalues
        if (IsPositive(traceSum)) {
            r1 = static_cast<T>(0.5) * (traceSum + rt);
            t  = static_cast<T>(1) / (r1);
            r2 = (A * t) * C - (B * t) * B;
        }
        else if (IsNegative(traceSum)) {
            r2 = static_cast<T>(0.5) * (traceSum - rt);
            t  = static_cast<T>(1) / (r2);
            r1 = (A * t) * C - (B * t) * B;
        }
        else {
            r1 =  static_cast<T>(0.5) * rt;
            r2 = -static_cast<T>(0.5) * rt;
        }
        xo_W[0] = r1;
        xo_W[1] = r2;

        // eigenvectors
        c = [traceDiff, rt]() {
            if (IsPositive(traceDiff)) {
                return (traceDiff + rt);
            }
            else {
                return (traceDiff - rt);
            }
        }();

        if (std::abs(c) > static_cast<T>(2) * std::abs(B)) {
            t = -static_cast<T>(2) * B / c;
            s = static_cast<T>(1) / std::sqrt(static_cast<T>(1) + t * t);
            c = t * s;
        }
        else if (FloatingPointTrait<double>::Equals(std::abs(B), 0.0)) {
            c = static_cast<T>(1);
            s = T{};
        }
        else {
            t = -static_cast<T>(0.5) * c / B;
            c = static_cast<T>(1) / std::sqrt(static_cast<T>(1) + t * t);
            s = t * c;
        }

        if (IsPositive(traceDiff)) {
            t =  c;
            c = -s;
            s =  t;
        }

        xo_U = { c, -s, 
				 s,  c };
    }

	/**
    * \brief perform polar decomposition (PD) on a given matrix, i.e. - A = R * S,
    *        where R is a rotation matrix in Givens form, S is symmetric describing deformations.
    *        (a granted negative sign on the small magnitude singular value)
    *
    * @param {matrix, in}  A
    * @param {matrix, out} R
    * @param {matrix, out} S
    **/
    template<typename T, class LAYOUT> void PolarDecomposition(const MatrixNM<T, 2, 2, LAYOUT>& xi_A, MatrixNM<T, 2, 2, LAYOUT>& xo_R, MatrixNM<T, 2, 2, LAYOUT>& xo_S) noexcept {
        const T x0{ xi_A(0, 0 ) + xi_A(1, 1) },
				x1{ xi_A(1, 0) - xi_A(0, 1) };
		T den{ std::sqrt(x0 * x0 + x1 * x1) },
		  c{ T(1) },
          s{};

        // R
        if (!FloatingPointTrait<double>::Equals(den, T{})) {
            den =  static_cast<T>(1) / den;
            c   =  x0 * den;
            s   = -x1 * den;
        }
        xo_R = { c, -s,
                 s,  c };

        // S
        xo_S = xi_A;
		static_for<0, 2>([&](std::size_t j) {
			const T tau1{ xo_S(0, j) },
					tau2{ xo_S(1, j) };
            xo_S(0, j) = c * tau1 - s * tau2;
            xo_S(1, j) = s * tau1 + c * tau2;
		});
    }

	// ----------------------------------------
	// --- special operations on 3x3 matrix ---
	// ----------------------------------------

	/**
    * \brief return the inverse of a given matrix (if its singular, return a null matrix)
    *
    * @param {matrix, in}  matrix
    * @param {matrix, out} matrix inverse
    **/
    template<typename T, class LAYOUT> constexpr MatrixNM<T, 3, 3, LAYOUT> Inverse(const MatrixNM<T, 3, 3, LAYOUT>& xi_mat) noexcept {
        // determinant
        const T det{ xi_mat.Determinant() };
		MatrixNM<T, 3, 3, LAYOUT> xo_inv{ T{} };

        if (!IsZero(det)) {
            const T detInv{ static_cast<T>(1) / det };

            // value at each position (values are the determinant of inner minors)
            xo_inv = {  (xi_mat(1, 1) * xi_mat(2, 2) - xi_mat(2, 1) * xi_mat(1, 2)),
                       -(xi_mat(1, 0) * xi_mat(2, 2) - xi_mat(1, 2) * xi_mat(2, 0)),
                        (xi_mat(1, 0) * xi_mat(2, 1) - xi_mat(2, 0) * xi_mat(1, 1)),
                       -(xi_mat(0, 1) * xi_mat(2, 2) - xi_mat(0, 2) * xi_mat(2, 1)),
                        (xi_mat(0, 0) * xi_mat(2, 2) - xi_mat(0, 2) * xi_mat(2, 0)),
                       -(xi_mat(0, 0) * xi_mat(2, 1) - xi_mat(2, 0) * xi_mat(0, 1)),
                        (xi_mat(0, 1) * xi_mat(1, 2) - xi_mat(0, 2) * xi_mat(1, 1)),
                       -(xi_mat(0, 0) * xi_mat(1, 2) - xi_mat(1, 0) * xi_mat(0, 2)),
                        (xi_mat(0, 0) * xi_mat(1, 1) - xi_mat(1, 0) * xi_mat(0, 1)) };

            // inverted matrix
            xo_inv *= detInv;
            if constexpr (LAYOUT::m_RowMajor) {
                xo_inv = xo_inv.Transpose();
            }   
        }

		return xo_inv;
    }

	/**
    * \brief test if 3x3 matrix is orthonormal
    *
    * @param {matrix, in}  matrix
    * @param {bool,   out} true if a 3x3 matrix is orthonormal
    **/
    template<typename T, class LAYOUT> constexpr bool IsOrthonormal(const MatrixNM<T, 3, 3, LAYOUT>& xi_mat) noexcept {
        const VectorN<T, 3> x = xi_mat(Column(0));

        if (!FloatingPointTrait<T>::Equals(Dot(x, x), static_cast<T>(1))) {
            return false;
        }

        const VectorN<T, 3> y = xi_mat(Column(1));
        if (!FloatingPointTrait<T>::Equals(Dot(y, y), static_cast<T>(1))) {
            return false;
        }

        const VectorN<T, 3> z = xi_mat(Column(2));
        if (!FloatingPointTrait<T>::Equals(Dot(z, z), static_cast<T>(1))) {
            return false;
        }

        if (!IsZero(Dot(x, y))) {
            return false;
        }

        if (!IsZero(Dot(x, z))) {
            return false;
        }

        if (!IsZero(Dot(y, z))) {
            return false;
        }

        return true;
    }

	/**
    * \brief calculate eigenvalues of a real 3x3 matrix (using its characteristic polynomial matrix)
    *
    * @param {matrix, in}  matrix
    * @param {vector, out} eigenvalues (empty vector if eigenvalues are complex)
    **/
    template<typename T, class LAYOUT> void Eigenvalues(const MatrixNM<T, 3, 3, LAYOUT>& xi_mat, VectorN<T, 3>& xo_eig) noexcept {
		const T trA{ xi_mat(0, 0) + xi_mat(1, 1) + xi_mat(2, 2) },
                detA{ Determinant(xi_mat) },
			    cofSum{ xi_mat(1, 1) * xi_mat(2, 2) - xi_mat(2, 1) * xi_mat(1, 2) +    // 1st minor
					    xi_mat(0, 0) * xi_mat(2, 2) - xi_mat(2, 0) * xi_mat(0, 2) +    // 2nd minor
					    xi_mat(0, 0) * xi_mat(1, 1) - xi_mat(1, 0) * xi_mat(0, 1) };   // 3rd minor

        std::vector<T> roots(6);
        int32_t numOfRoots = Numeric::SolveCubic<T>(-trA, cofSum, -detA, roots);
        
		xo_eig = { T{}, T{}, T{} };
		if (numOfRoots > 0) {
			xo_eig = { roots[0], roots[2], roots[4] };
		}
    }

	/**
    * \brief given a real 3x3 SYMMETRIC matrix, return its eigenvalues and eigenvectors
	*        Remark: work good in case eigenvalues are well separated.
    * 
    * @param {matrix, in}  matrix
    * @param {vector, out} eigenvalues (empty vector if eigenvalues are complex)
    * @param {matrix, out} matrix whos columns/rows are the normalized eigenvectors (empty matrix if eigenvalues are complex)
    **/
	template<typename T, class LAYOUT> void SymmetricEigen(const MatrixNM<T, 3, 3, LAYOUT>& xi_mat, VectorN<T, 3>& xo_values, MatrixNM<T, 3, 3, LAYOUT>& xo_vectors) noexcept {
        Eigenvalues(xi_mat, xo_values);
		
		xo_vectors = { T{}, T{}, T{},
					   T{}, T{}, T{},
					   T{}, T{}, T{} };

        // Compute eigenvectors only if the eigenvalues are real
        if (IsPositive(xo_values.Sum())) {
			static_for<0, 3>([&](std::size_t i) {
				const VectorN<T, 3> r1{ xi_mat(0, 0) - xo_values[i], xi_mat(0, 1),                xi_mat(0, 2)			        },
									r2{ xi_mat(0, 1),				   xi_mat(1, 1) - xo_values[i], xi_mat(1, 2)				},
									r3{ xi_mat(0, 2),				   xi_mat(1, 2),			    xi_mat(2, 2) - xo_values[i] },
									e1{ Cross(r1, r2) };
				VectorN<T, 3> e2{ Cross(r2, r3) },
							  e3{ Cross(r3, r1) };
               
				// Make e2 and e3 point in the same direction as e1
				if (IsNegative(Dot(e1, e2))) {
					e2 = -e2;
				}
				if (IsNegative(Dot(e1, e3))) {
					e3 = -e3;
				}

				VectorN<T, 3> eigvec{ e1 + e2 + e3 };
				eigvec.Normalize();

				xo_vectors(i, 0) = eigvec[0];
				xo_vectors(i, 1) = eigvec[1];
				xo_vectors(i, 2) = eigvec[2];
			});
        }
    }

	// ----------------------------------------
	// --- special operations on 4x4 matrix ---
	// ----------------------------------------

	/**
    * \brief return the inverse of a given matrix (if its singular, return a null matrix)
    *
    * @param {matrix, in}  matrix
    * @param {matrix, out} matrix inverse
    **/
    template<typename T, class LAYOUT> constexpr MatrixNM<T, 4, 4, LAYOUT> Inverse(const MatrixNM<T, 4, 4, LAYOUT>& xi_mat) noexcept {
        // determinant
        const T det{ xi_mat.Determinant() };
        
		MatrixNM<T, 4, 4, LAYOUT> xo_inv(T{});
        if (!IsZero(det)) {
            const T detInv{ static_cast<T>(1) / det };

            // value at each position (values are the determinant of inner minors)
            xo_inv = { xi_mat(1, 1) * xi_mat(2, 2) * xi_mat(3, 3) - xi_mat(1, 1) * xi_mat(2, 3) * xi_mat(3, 2) -
                       xi_mat(2, 1) * xi_mat(1, 2) * xi_mat(3, 3) + xi_mat(2, 1) * xi_mat(1, 3) * xi_mat(3, 2) +
                       xi_mat(3, 1) * xi_mat(1, 2) * xi_mat(2, 3) - xi_mat(3, 1) * xi_mat(1, 3) * xi_mat(2, 2),

                      -xi_mat(0, 1) * xi_mat(2, 2) * xi_mat(3, 3) + xi_mat(0, 1) * xi_mat(2, 3) * xi_mat(3, 2) +
                       xi_mat(2, 1) * xi_mat(0, 2) * xi_mat(3, 3) - xi_mat(2, 1) * xi_mat(0, 3) * xi_mat(3, 2) -
                       xi_mat(3, 1) * xi_mat(0, 2) * xi_mat(2, 3) + xi_mat(3, 1) * xi_mat(0, 3) * xi_mat(2, 2),

                       xi_mat(0, 1) * xi_mat(1, 2) * xi_mat(3, 3) - xi_mat(0, 1) * xi_mat(1, 3) * xi_mat(3, 2) -
                       xi_mat(1, 1) * xi_mat(0, 2) * xi_mat(3, 3) + xi_mat(1, 1) * xi_mat(0, 3) * xi_mat(3, 2) +
                       xi_mat(3, 1) * xi_mat(0, 2) * xi_mat(1, 3) - xi_mat(3, 1) * xi_mat(0, 3) * xi_mat(1, 2),

                      -xi_mat(0, 1) * xi_mat(1, 2) * xi_mat(2, 3) + xi_mat(0, 1) * xi_mat(1, 3) * xi_mat(2, 2) +
                       xi_mat(1, 1) * xi_mat(0, 2) * xi_mat(2, 3) - xi_mat(1, 1) * xi_mat(0, 3) * xi_mat(2, 2) -
                       xi_mat(2, 1) * xi_mat(0, 2) * xi_mat(1, 3) + xi_mat(2, 1) * xi_mat(0, 3) * xi_mat(1, 2),

                      -xi_mat(1, 0) * xi_mat(2, 2) * xi_mat(3, 3) + xi_mat(1, 0) * xi_mat(2, 3) * xi_mat(3, 2) +
                       xi_mat(2, 0) * xi_mat(1, 2) * xi_mat(3, 3) - xi_mat(2, 0) * xi_mat(1, 3) * xi_mat(3, 2) -
                       xi_mat(3, 0) * xi_mat(1, 2) * xi_mat(2, 3) + xi_mat(3, 0) * xi_mat(1, 3) * xi_mat(2, 2),

                       xi_mat(0, 0) * xi_mat(2, 2) * xi_mat(3, 3) - xi_mat(0, 0) * xi_mat(2, 3) * xi_mat(3, 2) -
                       xi_mat(2, 0) * xi_mat(0, 2) * xi_mat(3, 3) + xi_mat(2, 0) * xi_mat(0, 3) * xi_mat(3, 2) +
                       xi_mat(3, 0) * xi_mat(0, 2) * xi_mat(2, 3) - xi_mat(3, 0) * xi_mat(0, 3) * xi_mat(2, 2),

                      -xi_mat(0, 0) * xi_mat(1, 2) * xi_mat(3, 3) + xi_mat(0, 0) * xi_mat(1, 3) * xi_mat(3, 2) +
                       xi_mat(1, 0) * xi_mat(0, 2) * xi_mat(3, 3) - xi_mat(1, 0) * xi_mat(0, 3) * xi_mat(3, 2) -
                       xi_mat(3, 0) * xi_mat(0, 2) * xi_mat(1, 3) + xi_mat(3, 0) * xi_mat(0, 3) * xi_mat(1, 2),

                       xi_mat(0, 0) * xi_mat(1, 2) * xi_mat(2, 3) - xi_mat(0, 0) * xi_mat(1, 3) * xi_mat(2, 2) -
                       xi_mat(1, 0) * xi_mat(0, 2) * xi_mat(2, 3) + xi_mat(1, 0) * xi_mat(0, 3) * xi_mat(2, 2) +
                       xi_mat(2, 0) * xi_mat(0, 2) * xi_mat(1, 3) - xi_mat(2, 0) * xi_mat(0, 3) * xi_mat(1, 2),
                       
                       xi_mat(1, 0) * xi_mat(2, 1) * xi_mat(3, 3) - xi_mat(1, 0) * xi_mat(2, 3) * xi_mat(3, 1) -
                       xi_mat(2, 0) * xi_mat(1, 1) * xi_mat(3, 3) + xi_mat(2, 0) * xi_mat(1, 3) * xi_mat(3, 1) +
                       xi_mat(3, 0) * xi_mat(1, 1) * xi_mat(2, 3) - xi_mat(3, 0) * xi_mat(1, 3) * xi_mat(2, 1), 
                       
                      -xi_mat(0, 0) * xi_mat(2, 1) * xi_mat(3, 3) + xi_mat(0, 0) * xi_mat(2, 3) * xi_mat(3, 1) +
                       xi_mat(2, 0) * xi_mat(0, 1) * xi_mat(3, 3) - xi_mat(2, 0) * xi_mat(0, 3) * xi_mat(3, 1) -
                       xi_mat(3, 0) * xi_mat(0, 1) * xi_mat(2, 3) + xi_mat(3, 0) * xi_mat(0, 3) * xi_mat(2, 1),
                       
                       xi_mat(0, 0) * xi_mat(1, 1) * xi_mat(3, 3) - xi_mat(0, 0) * xi_mat(1, 3) * xi_mat(3, 1) -
                       xi_mat(1, 0) * xi_mat(0, 1) * xi_mat(3, 3) + xi_mat(1, 0) * xi_mat(0, 3) * xi_mat(3, 1) +
                       xi_mat(3, 0) * xi_mat(0, 1) * xi_mat(1, 3) - xi_mat(3, 0) * xi_mat(0, 3) * xi_mat(1, 1),

                      -xi_mat(0, 0) * xi_mat(1, 1) * xi_mat(2, 3) + xi_mat(0, 0) * xi_mat(1, 3) * xi_mat(2, 1) +
                       xi_mat(1, 0) * xi_mat(0, 1) * xi_mat(2, 3) - xi_mat(1, 0) * xi_mat(0, 3) * xi_mat(2, 1) -
                       xi_mat(2, 0) * xi_mat(0, 1) * xi_mat(1, 3) + xi_mat(2, 0) * xi_mat(0, 3) * xi_mat(1, 1),
                      
                      -xi_mat(1, 0) * xi_mat(2, 1) * xi_mat(3, 2) + xi_mat(1, 0) * xi_mat(2, 2) * xi_mat(3, 1) +
                       xi_mat(2, 0) * xi_mat(1, 1) * xi_mat(3, 2) - xi_mat(2, 0) * xi_mat(1, 2) * xi_mat(3, 1) -
                       xi_mat(3, 0) * xi_mat(1, 1) * xi_mat(2, 2) + xi_mat(3, 0) * xi_mat(1, 2) * xi_mat(2, 1),
                      
                       xi_mat(0, 0) * xi_mat(2, 1) * xi_mat(3, 2) - xi_mat(0, 0) * xi_mat(2, 2) * xi_mat(3, 1) -
                       xi_mat(2, 0) * xi_mat(0, 1) * xi_mat(3, 2) + xi_mat(2, 0) * xi_mat(0, 2) * xi_mat(3, 1) +
                       xi_mat(3, 0) * xi_mat(0, 1) * xi_mat(2, 2) - xi_mat(3, 0) * xi_mat(0, 2) * xi_mat(2, 1),
                      
                      -xi_mat(0, 0) * xi_mat(1, 1) * xi_mat(3, 2) + xi_mat(0, 0) * xi_mat(1, 2) * xi_mat(3, 1) +
                       xi_mat(1, 0) * xi_mat(0, 1) * xi_mat(3, 2) - xi_mat(1, 0) * xi_mat(0, 2) * xi_mat(3, 1) -
                       xi_mat(3, 0) * xi_mat(0, 1) * xi_mat(1, 2) + xi_mat(3, 0) * xi_mat(0, 2) * xi_mat(1, 1),
                      
                       xi_mat(0, 0) * xi_mat(1, 1) * xi_mat(2, 2) - xi_mat(0, 0) * xi_mat(1, 2) * xi_mat(2, 1) -
                       xi_mat(1, 0) * xi_mat(0, 1) * xi_mat(2, 2) + xi_mat(1, 0) * xi_mat(0, 2) * xi_mat(2, 1) +
                       xi_mat(2, 0) * xi_mat(0, 1) * xi_mat(1, 2) - xi_mat(2, 0) * xi_mat(0, 2) * xi_mat(1, 1) };

            // inverted matrix
            xo_inv *= detInv;
            if constexpr (LAYOUT::m_RowMajor) {
                xo_inv = xo_inv.Transpose();
            }
        }

        return xo_inv;
    }

	/**
    * \brief efficient inverse of Affine matrix representing rigid transformation
    *        (assuming matrix is not singular)
    *
    * @param {matrix, out} Affine inverse
    * @param {matrix, in}  Affine
    **/
    template<typename T, class LAYOUT> constexpr MatrixNM<T, 4, 4, LAYOUT> AffineInverse(const MatrixNM<T, 4, 4, LAYOUT>& xi_mat) noexcept {
        MatrixNM<T, 4, 4, LAYOUT> xo_inverse;

        // transpose DCM
		static_for<0, 3>([&](std::size_t icol) {
			static_for<0, 3>([&](std::size_t irow) {
                xo_inverse(icol, irow) = xi_mat(irow, icol);
			});
		});

        // multiply - translation by DCM transpose
        VectorN<T, 3> x, y, z;
        
        if constexpr (LAYOUT::m_RowMajor) {
            x = { xi_mat(0, 0), xi_mat(0, 1), xi_mat(0, 2) };
            y = { xi_mat(1, 0), xi_mat(1, 1), xi_mat(1, 2) };
            z = { xi_mat(2, 0), xi_mat(2, 1), xi_mat(2, 2) };
        }
        else {
            x = { xi_mat(0, 0), xi_mat(1, 0), xi_mat(2, 0) };
            y = { xi_mat(0, 1), xi_mat(1, 1), xi_mat(2, 1) };
            z = { xi_mat(0, 2), xi_mat(1, 2), xi_mat(2, 2) };
        }

        xo_inverse(3, 0) = -Dot(z, x);
        xo_inverse(3, 1) = -Dot(z, y);
        xo_inverse(3, 2) = -Dot(z, z);
        xo_inverse(3, 3) = static_cast<T>(1);

        return xo_inverse;
    }

	template<typename T, class LAYOUT> constexpr MatrixNM<T, 4, 4, LAYOUT>& AffineInverse(MatrixNM<T, 4, 4, LAYOUT>&& xi_mat) noexcept {
		// multiply - translation by DCM transpose
		VectorN<T, 3> x, y, z;

		if constexpr (xi_mat.m_RowMajor) {
			x = { xi_mat(0, 0), xi_mat(0, 1), xi_mat(0, 2) };
			y = { xi_mat(1, 0), xi_mat(1, 1), xi_mat(1, 2) };
			z = { xi_mat(2, 0), xi_mat(2, 1), xi_mat(2, 2) };
		}
		else {
			x = { xi_mat(0, 0), xi_mat(1, 0), xi_mat(2, 0) };
			y = { xi_mat(0, 1), xi_mat(1, 1), xi_mat(2, 1) };
			z = { xi_mat(0, 2), xi_mat(1, 2), xi_mat(2, 2) };
		}

		// transpose DCM
		static_for<0, 3>([&](std::size_t icol) {
			static_for<0, 3>([&](std::size_t irow) {
				xi_mat(icol, irow) = xi_mat(irow, icol);
			});
		});

		xi_mat(3, 0) = -Dot(z, x);
		xi_mat(3, 1) = -Dot(z, y);
		xi_mat(3, 2) = -Dot(z, z);
		xi_mat(3, 3) = static_cast<T>(1);

		return std::move(xi_mat);
	}

	/**
    * \brief given field of view, aspect ratio and clamp distance, return the perspective matrix
    *
    * @param {matrix,   out} perspective matrix
    * @param {T,        in}  field of view [rad]
    * @param {T,        in}  aspect ratio
    * @param {T,        in}  near field clamp (Z axis)
    * @param {T,        in}  far field clamp (Z axis)
    * @param {Handness, in}  coordinate system type
    **/
    template<typename T, class LAYOUT> void CreatePerspectiveMatrix(MatrixNM<T, 4, 4, LAYOUT>& xio_mat, const T xi_fov_y, const T xi_aspect_ratio, 
																	const T xi_near, const T xi_far, const Handness& xi_type) noexcept {
        const T y{ static_cast<T>(1) / std::tan(xi_fov_y * static_cast<T>(0.5)) },
				x{ y / xi_aspect_ratio },
				zdist{ xi_near - xi_far },
				zfar_per_zdist{ xi_far / zdist };

        if constexpr(LAYOUT::m_RowMajor) {
            xio_mat = { x,   T{}, T{},										T{},
                        T{}, y,   T{},										T{},
                        T{}, T{}, zfar_per_zdist * static_cast<T>(xi_type), -static_cast<T>(xi_type),
                        T{}, T{}, T(2) * xi_near * zfar_per_zdist,			T{} };
        }
        else {
            xio_mat = { x,   T{}, T{},										T{}, 
                        T{}, y,   T{},										T{}, 
                        T{}, T{}, zfar_per_zdist * static_cast<T>(xi_type), static_cast<T>(2) * xi_near * zfar_per_zdist,
                        T{}, T{}, -static_cast<T>(xi_type),                 T{} };
            }
    }

    /**
    * \brief given dimensions and distance, return the orthographic projection matrix
    *
    * @param {matrix,   out} orthographic matrix
    * @param {T,        in}  left boundary
    * @param {T,        in}  right boundary
    * @param {T,        in}  bottom boundary
    * @param {T,        in}  top boundary
    * @param {T,        in}  near field clamp (Z axis)
    * @param {T,        in}  far field clamp (Z axis)
    * @param {Handness, in}  coordinate system type
    **/
    template<typename T, class LAYOUT> void CreateOrthographicMatrix(MatrixNM<T, 4, 4, LAYOUT>& xio_mat,
																	 const T xi_left, const T xi_right, const T xi_bottom, const T xi_top,
																	 const T xi_near, const T xi_far, const Handness& xi_type) noexcept {
        const T sideDiffInv{ static_cast<T>(1)  / (xi_right - xi_left) },
				upDiffInv{ static_cast<T>(1)    / (xi_top - xi_bottom) },
				depthDiffInv{ static_cast<T>(1) / (xi_far - xi_near) };
        
        if constexpr (LAYOUT::m_RowMajor) {
            xio_mat = { static_cast<T>(2) * sideDiffInv,     T{},                                T{},														  T{},
                        T{},                                 static_cast<T>(2) * upDiffInv,     T{},														  T{},
                        T{},                                 T{},                               -static_cast<T>(xi_type) * static_cast<T>(2) * depthDiffInv,  T{},
                        -(xi_right + xi_left) * sideDiffInv, -(xi_top + xi_bottom) * upDiffInv, -(xi_far + xi_near) * depthDiffInv,							  static_cast<T>(1) };
        }
        else {
            xio_mat = { static_cast<T>(2) * sideDiffInv, T{},              T{},															-(xi_right + xi_left) * sideDiffInv,
                        T{},                static_cast<T>(2) * upDiffInv, T{},															-(xi_top + xi_bottom) * upDiffInv,              
                        T{},                T{},						   -static_cast<T>(xi_type) * static_cast<T>(2) * depthDiffInv, -(xi_far + xi_near) * depthDiffInv,             
                        T{},                T{},						   T{},															static_cast<T>(1) };
        }
    }

    /**
    * \brief given object position, camera position and "up" direction (as normalized vector), return the appropriate "view" matrix
    *
    * @param {matrix,   out} look at matrix
    * @param {matrix,   in}  object position (the "to look at" object; "target")
    * @param {matrix,   in}  camera ("eye") position
    * @param {matrix,   in}  "up" direction
    * @param {Handness, in}  coordinate system type
    **/
    template<typename T, class LAYOUT> void CreateLookAtMatrix(MatrixNM<T, 4, 4, LAYOUT>& xio_mat, const VectorN<T, 3>& xi_object, 
															   const VectorN<T, 3>& xi_eye, const VectorN<T, 3>& xi_up, const Handness& xi_type) noexcept {
        const T type{ static_cast<T>(xi_type) };
        const VectorN<T, 3> ZZ{ -type * (xi_object - xi_eye).Normalize() },
                            XX{ -type * Cross(xi_up, ZZ).Normalize() },
                            YY{ Cross(ZZ, XX) };

        if constexpr (LAYOUT::m_RowMajor) {
            xio_mat = { XX.X(),                 XX.Y(),           XX.Z(),                 T{},
                        YY.X(),                 YY.Y(),           YY.Z(),                 T{},
                        ZZ.X(),                 ZZ.Y(),           ZZ.Z(),                 T{},
                        type * Dot(XX, xi_eye), -Dot(YY, xi_eye), type * Dot(ZZ, xi_eye), static_cast<T>(1) };
        }
        else {
            xio_mat = { XX.X(), YY.X(), ZZ.X(), type * Dot(XX, xi_eye),
                        XX.Y(), YY.Y(), ZZ.Y(), -Dot(YY, xi_eye),
                        XX.Z(), YY.Z(), ZZ.Z(), type * Dot(ZZ, xi_eye),
                        T{},    T{},    T{},    static_cast<T>(1) };
        }
    }

	// ------------------------------------
	// --- various geometric operations ---
	// ------------------------------------

	/**
	* \brief calculate ray-plane intersection point
	* 
	* @param {vector,           in}  plane normal (normalized unit vector)
	* @param {vector,           in}  plane point
	* @param {vector,           in}  ray origin
	* @param {vector,           in}  ray direction  (normalized unit vector)
	* @param {optional<vector>, out} intersection point
	**/
	template<typename T> std::optional<VectorN<T, 3>> intersectPlane(const VectorN<T, 3>& xi_normal, const VectorN<T, 3>& xi_plane_point, 
																	 const VectorN<T, 3>& xi_ray_origin, const VectorN<T, 3>& xi_ray_direction) noexcept {
		T denom{ Dot(xi_normal, xi_ray_direction) };

		if (std::fabs(denom) > Numeric::FloatingPointTrait<T>::epsilon()) {
			const VectorN<T, 3> p0l0{ xi_plane_point - xi_ray_origin };
			const T t{ Dot(p0l0, xi_normal) / denom };
			if (IsPositive(t)) {
				return ( xi_ray_origin + t * xi_ray_direction );
			}
		}

		return {};
	}

	/**
	* \brief given two 'plane'-rays (as origin & direction), find their intersection point
	* 
	* @param vector, in}  plane-ray #1 origin
	* @param vector, in}  plane-ray #1 direction (normalized unit vector)
	* @param vector, in}  plane-ray #2 origin
	* @param vector, in}  plane-ray #3 direction (normalized unit vector)
	* @param vector, out} intersection point
	**/
	template<typename T> void intersectPRays(const VectorN<T, 3>& xi_ray1_origin, const VectorN<T, 3>& xi_ray1_direction,
											 const VectorN<T, 3>& xi_ray2_origin, const VectorN<T, 3>& xi_ray2_direction,
											 VectorN<T, 3>& xo_intersection) noexcept {
		// calculate interpolant's along rays to intersection
		VectorN<T, 3>&& w0{ xi_ray1_origin - xi_ray2_origin };
		//const VectorN<T, 3> w0{ xi_ray1_origin - xi_ray2_origin };
		const float a{ Dot(xi_ray1_direction, xi_ray1_direction) },
					b{ Dot(xi_ray1_direction, xi_ray2_direction) },
					c{ Dot(xi_ray2_direction, xi_ray2_direction) },
					d{ Dot(xi_ray1_direction, w0) },
					e{ Dot(xi_ray2_direction, w0) },
					den{ a * c - b * b },
					s{ (b * e - c * d) / den },	// interpolant along ray #1
					t{ (a * e - b * d) / den };	// interpolant along ray #2

		// intersection point
		xo_intersection = xi_ray1_origin + s * xi_ray1_direction;
	}

	/**
    * \brief return the angle between three points (law of cosine)
    *
    * @param {vector, in}  point 1 (edge point)
    * @param {vector, in}  point 2 (middle point)
    * @param {vector, in}  point 3 (edge point)
    * @param {T,      out} angle [rad]
    **/
	template<typename T> T AngleBetweenThreePoints(const VectorN<T, 3>& xi_a, const VectorN<T, 3>& xi_b, const VectorN<T, 3>& xi_c) noexcept {
		VectorN<T, 3> BminusA{ xi_b - xi_a },
					  CminusB{ xi_c - xi_b },
					  CminusA{ xi_c - xi_a };

		const T ABmagSqr{ Dot(BminusA, BminusA) },
				BCmagSqr{ Dot(CminusB, CminusB) },
				ACmagSqr{ Dot(CminusA, CminusA) },
				ABmag{ std::sqrt(ABmagSqr) },
				BCmag{ std::sqrt(BCmagSqr) },
				cosAngle{ (ABmagSqr + BCmagSqr - ACmagSqr) / (static_cast<T>(2.0) * ABmag * BCmag) };

		// have we gotten a false trigonometric value?
		return ((std::abs(cosAngle) <= static_cast<T>(1.0)) ? (std::acos(cosAngle) * Numeric::Constants<T>::RADIAN2DEGREE()) : T{});
    }

	/**
	* \brief given a point, project it onto a plane (given by point and normal)
	* 
	* @param {vector, in}  plane normal (normalized unit vector)
	* @param {vector, in}  plane point
	* @param {vector, in}  point to be projected
	* @param {vector, out} projected point
	**/
	template<typename T> inline void projectOnPlane(const VectorN<T, 3>& xi_normal, const VectorN<T, 3>& xi_plane_point, const VectorN<T, 3>& xi_point_to_project, VectorN<T, 3>& xo_projected_point) noexcept {
		const VectorN<T, 3> V{ xi_point_to_project - xi_plane_point },		// line from project point to plane
						    Vn{ Dot(V, xi_normal) * xi_normal };			// line projected on normal
		xo_projected_point = xi_point_to_project - Vn;
	}
};