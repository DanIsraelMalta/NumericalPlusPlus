/**
* A quaternion object.
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
#include "3D.h"
#include <array>
#include <vector>
#include <functional>
#include <future>
#include <assert.h>
#include <string>

namespace Numeric {

	/**
	* \brief quaternion {x, y, z, w}
	*
	* @param {T, in} underlying type
	**/
	template<typename T> class Quaternion {
		static_assert(std::is_arithmetic<T>::value, "Quaternion<T> - T must be of numerical type.");

		// properties
		private:
			alignas(std::aligned_storage<sizeof(T), alignof(T)>::type) T m_data[4];

		// constructors
		public:

			// default constructor
			explicit constexpr Quaternion() : m_data{ T{}, T{}, T{}, static_cast<T>(1.0) } {}

			// construct using a single value
			explicit constexpr Quaternion(const T xi_value) {
				static_for<0, 4>([&](std::size_t i) {
					m_data[i] = xi_value;
				});
			}

			// construct using a VectorN
			explicit constexpr Quaternion(const VectorN<T, 4>& xi_vector) {
				static_for<0, 4>([&](std::size_t i) {
					m_data[i] = xi_vector[i];
				});
			}

			explicit constexpr Quaternion(VectorN<T, 4>&& xi_vector) {
				static_for<0, 4>([&](std::size_t i) {
					m_data[i] = std::move(xi_vector[i]);
				});
			}

			// construct from a moveable array (usage: Quaternion<float> v{{1.0f, 2.0f, 3.0f, 0.5f}}; )
			explicit constexpr Quaternion(T(&&xi_array)[4]) : m_data(std::make_move_iterator(std::begin(xi_array)), std::make_move_iterator(std::end(xi_array))) {}

			// construct from a parameter pack (usage: VectorN<int, 3> v(0, 1, 2); )
			template<class ... N> explicit constexpr Quaternion(T first, N&&... values) : Quaternion{ first, std::forward<T>(static_cast<T>(values))... } {}

			// construct using list initializer
			explicit constexpr Quaternion(const std::initializer_list<T>&& xi_list) {
				// throws "array iterator + offset out of range" if more then 4 elements are entered
				std::move(std::begin(xi_list), std::end(xi_list), m_data);
			}

			// construct using std::array
			explicit constexpr Quaternion(const std::array<T, 4>& xi_arr) {
				std::move(std::begin(xi_arr), std::end(xi_arr), m_data);
			}

			explicit constexpr Quaternion(std::array<T, 4>&& xi_arr) {
				std::move(std::begin(xi_arr), std::end(xi_arr), m_data);
			}

			// construct using std::vector
			explicit constexpr Quaternion(const std::vector<T>& xi_vec) {
				assert(4 == xi_vec.size());
				std::move(std::begin(xi_vec), std::end(xi_vec), m_data);
			}

			explicit constexpr Quaternion(std::vector<T>&& xi_vec) {
				assert(4 == xi_vec.size());
				std::move(std::begin(xi_vec), std::end(xi_vec), m_data);
			}

			// copy semantics
			Quaternion(const Quaternion&)			 = default;
			Quaternion& operator=(const Quaternion&) = default;

			// move semantics
			Quaternion(Quaternion&&)            noexcept = default;
			Quaternion& operator=(Quaternion&&) noexcept = default;

			// construct from an axis (normalized) & angle [rad]
			explicit constexpr Quaternion(const VectorN<T, 3>& xi_axis, const T xi_angle) {
				const T halfAngle{ static_cast<T>(0.5) * xi_angle },
						WW{ std::cos(halfAngle) },
						sinHalf{ std::sin(halfAngle) };

				const VectorN<T, 3> axis{ sinHalf * xi_axis };

				// in case some components are close to zero
				VectorN<T, 4> quat{ axis.X(), axis.Y(), axis.Z(), WW };
				quat.Normalize();

				m_data[0] = quat[0];
				m_data[1] = quat[1];
				m_data[2] = quat[2];
				m_data[3] = quat[3];
			}

			explicit constexpr Quaternion(VectorN<T, 3>&& xi_axis, const T xi_angle) {
				const T halfAngle{ static_cast<T>(0.5) * xi_angle },
					    WW{ std::cos(halfAngle) },
					    sinHalf{ std::sin(halfAngle) };

				xi_axis *= sinHalf;

				// in case some components are close to zero
				VectorN<T, 4> quat{ xi_axis.X(), xi_axis.Y(), xi_axis.Z(), WW };
				quat.Normalize();

				m_data[0] = quat[0];
				m_data[1] = quat[1];
				m_data[2] = quat[2];
				m_data[3] = quat[3];
			}

			// construct from two (normalized) vectors (quaternion will describe rotation from first to second)
			explicit constexpr Quaternion(const VectorN<T, 3>& xi_from, const VectorN<T, 3>& xi_to) {
				const T cosAngle{ Dot(xi_from, xi_to) },
						angle{ std::acos(cosAngle) };

				VectorN<T, 3> axis{ Cross(xi_from, xi_to) };
				axis.Normalize();

				const Quaternion<T> fromAxisAngle(axis, angle);
				m_data[0] = fromAxisAngle[0];
				m_data[1] = fromAxisAngle[1];
				m_data[2] = fromAxisAngle[2];
				m_data[3] = fromAxisAngle[3];
			}

			explicit constexpr Quaternion(VectorN<T, 3>&& xi_from, VectorN<T, 3>&& xi_to) {
				const T cosAngle{ Dot(xi_from, xi_to) },
						angle{ std::acos(cosAngle) };

				xi_from = Cross(xi_from, std::move(xi_to));
				xi_from.Normalize();

				const Quaternion<T> fromAxisAngle(xi_from, angle);
				m_data[0] = fromAxisAngle[0];
				m_data[1] = fromAxisAngle[1];
				m_data[2] = fromAxisAngle[2];
				m_data[3] = fromAxisAngle[3];
			}

			// construct from rotation matrix (3x3)
			explicit constexpr Quaternion(const MatrixNM<T, 3, 3>& xi_mat) {
				// locals
				const T Tr{ static_cast<T>(1) + xi_mat(0, 0) + xi_mat(1, 1) + xi_mat(2, 2) };
				T S, Sinv, x, y, z, w;

				// to avoid large distortions (DCM trace should be positive)
				if (Numeric::IsPositive(Tr)) {
					S = static_cast<T>(2) * std::sqrt(Tr);
					Sinv = static_cast<T>(1) / S;

					w = T(0.25) * S;
					x = (xi_mat(2, 1) - xi_mat(1, 2)) * Sinv;
					y = (xi_mat(0, 2) - xi_mat(2, 0)) * Sinv;
					z = (xi_mat(1, 0) - xi_mat(0, 1)) * Sinv;
				}   // DCM trace is zero, find largest diagonal element, and build quaternion using it
				else {   // Column 0
					if ((xi_mat(0, 0) > xi_mat(1, 1)) &&
						(xi_mat(0, 0) > xi_mat(2, 2))) {
						S = static_cast<T>(2) * std::sqrt(static_cast<T>(1) + xi_mat(0, 0) - xi_mat(1, 1) - xi_mat(2, 2));
						Sinv = T(1) / S;

						x = static_cast<T>(0.25) * S;
						y = (xi_mat(1, 0) + xi_mat(0, 1)) * Sinv;
						z = (xi_mat(0, 2) + xi_mat(2, 0)) * Sinv;
						w = (xi_mat(2, 1) - xi_mat(1, 2)) * Sinv;
					}   // Column 1
					else if (xi_mat(1, 1) > xi_mat(2, 2))
					{
						S = static_cast<T>(2) * std::sqrt(static_cast<T>(1) + xi_mat(1, 1) - xi_mat(0, 0) - xi_mat(2, 2));
						Sinv = static_cast<T>(1) / S;

						x = (xi_mat(1, 0) + xi_mat(0, 1)) * Sinv;
						y = static_cast<T>(0.25) * S;
						z = (xi_mat(2, 1) + xi_mat(1, 2)) * Sinv;
						w = (xi_mat(0, 2) - xi_mat(2, 0)) * Sinv;
					}
					else
					{   // Column 2
						S = static_cast<T>(2) * std::sqrt(static_cast<T>(1) + xi_mat(2, 2) - xi_mat(0, 0) - xi_mat(1, 1));
						Sinv = static_cast<T>(1) / S;

						x = (xi_mat(0, 2) + xi_mat(2, 0)) * Sinv;
						y = (xi_mat(2, 1) + xi_mat(1, 2)) * Sinv;
						z = static_cast<T>(0.25) * S;
						w = (xi_mat(1, 0) - xi_mat(0, 1)) * Sinv;
					}
				}

				m_data[0] = -x;
				m_data[1] = -y;
				m_data[2] = -z;
				m_data[3] = w;
			}

		// assignment operations
		public:

			// assign an element
			constexpr Quaternion& operator=(const T xi_value) noexcept {
				static_for<0, 4>([&](std::size_t i) {
					m_data[i] = xi_value;
				});
				return *this;
			}

			// assign a VectorN
			constexpr Quaternion& operator=(const VectorN<T, 4>& xi_vector) noexcept {
				static_for<0, 4>([&](std::size_t i) {
					m_data[i] = xi_vector[i];
				});
				return *this;
			}

			// assign from a list
			constexpr Quaternion& operator=(const std::initializer_list<T>&& xi_list) {
				// throws "array iterator + offset out of range" if more then 4 elements are entered
				std::move(std::begin(xi_list), std::end(xi_list), m_data);
				return *this;
			}

			// assign from array
			constexpr Quaternion& operator=(const T(&&xi_array)[4]) {
				std::move(&xi_array[0], &xi_array[4], m_data);
				return *this;
			};

			// assign from std::array
			constexpr Quaternion& operator=(const std::array<T, 4>& xi_arr) {
				std::move(std::begin(xi_arr), std::end(xi_arr), m_data);
				return *this;
			}

			constexpr Quaternion& operator=(std::array<T, 4>&& xi_arr) {
				std::move(std::begin(xi_arr), std::end(xi_arr), m_data);
				xi_arr = nullptr;
				return *this;
			}

			// assign from std::vector
			constexpr Quaternion& operator=(const std::vector<T>& xi_vec) {
				assert(4 == xi_vec.size());
				std::move(std::begin(xi_vec), std::end(xi_vec), m_data);
				return *this;
			}

			constexpr Quaternion& operator=(std::vector<T>&& xi_vec) {
				assert(4 == xi_vec.size());
				std::move(std::begin(xi_vec), std::end(xi_vec), m_data);
				xi_vec = nullptr;
				return *this;
			}

			// assign from a matrix
			constexpr Quaternion& operator=(const MatrixNM<T, 3, 3>& xi_mat) {
				// locals
				const T Tr{ static_cast<T>(1) + xi_mat(0, 0) + xi_mat(1, 1) + xi_mat(2, 2) };
				T S, Sinv, x, y, z, w;

				// to avoid large distortions (DCM trace should be positive)
				if (Numeric::IsPositive(Tr)) {
					S = static_cast<T>(2) * std::sqrt(Tr);
					Sinv = static_cast<T>(1) / S;

					w = T(0.25) * S;
					x = (xi_mat(2, 1) - xi_mat(1, 2)) * Sinv;
					y = (xi_mat(0, 2) - xi_mat(2, 0)) * Sinv;
					z = (xi_mat(1, 0) - xi_mat(0, 1)) * Sinv;
				}   // DCM trace is zero, find largest diagonal element, and build quaternion using it
				else {   // Column 0
					if ((xi_mat(0, 0) > xi_mat(1, 1)) &&
						(xi_mat(0, 0) > xi_mat(2, 2))) {
						S = static_cast<T>(2) * std::sqrt(static_cast<T>(1) + xi_mat(0, 0) - xi_mat(1, 1) - xi_mat(2, 2));
						Sinv = T(1) / S;

						x = static_cast<T>(0.25) * S;
						y = (xi_mat(1, 0) + xi_mat(0, 1)) * Sinv;
						z = (xi_mat(0, 2) + xi_mat(2, 0)) * Sinv;
						w = (xi_mat(2, 1) - xi_mat(1, 2)) * Sinv;
					}   // Column 1
					else if (xi_mat(1, 1) > xi_mat(2, 2))
					{
						S = static_cast<T>(2) * std::sqrt(static_cast<T>(1) + xi_mat(1, 1) - xi_mat(0, 0) - xi_mat(2, 2));
						Sinv = static_cast<T>(1) / S;

						x = (xi_mat(1, 0) + xi_mat(0, 1)) * Sinv;
						y = static_cast<T>(0.25) * S;
						z = (xi_mat(2, 1) + xi_mat(1, 2)) * Sinv;
						w = (xi_mat(0, 2) - xi_mat(2, 0)) * Sinv;
					}
					else
					{   // Column 2
						S = static_cast<T>(2) * std::sqrt(static_cast<T>(1) + xi_mat(2, 2) - xi_mat(0, 0) - xi_mat(1, 1));
						Sinv = static_cast<T>(1) / S;

						x = (xi_mat(0, 2) + xi_mat(2, 0)) * Sinv;
						y = (xi_mat(2, 1) + xi_mat(1, 2)) * Sinv;
						z = static_cast<T>(0.25) * S;
						w = (xi_mat(1, 0) - xi_mat(0, 1)) * Sinv;
					}
				}

				m_data = { -x, -y, -z, w };
				return *this;
			}

		// iterators
		public:
			constexpr T* begin()	   { return &m_data[0]; }
			constexpr T* begin() const { return &m_data[0]; }

			constexpr T* end()	     { return begin() + 4; }
			constexpr T* end() const { return begin() + 4; }

		// set/get operations
		public:

			// '[]' element access
			constexpr T  operator[](const std::size_t i) const { return m_data[i]; }
			constexpr T& operator[](const std::size_t i)	     { return m_data[i]; }

			// extract quaternion x/y/z/w
			constexpr T  X() const { return m_data[0]; }
			constexpr T& X()       { return m_data[0]; }

			constexpr T  Y() const { return m_data[1]; }
			constexpr T& Y()	   { return m_data[1]; }

			constexpr T  Z() const { return m_data[2]; }
			constexpr T& Z()	   { return m_data[2]; }

			constexpr T  W() const { return m_data[3]; }
			constexpr T& W()	   { return m_data[3]; }

			// extract quaternion angle (assuming the quaternion is normalized)
			constexpr T Angle() const noexcept { return static_cast<T>(2) * std::acos(this->W()); }

			// extract quaternion (normalized) axis (assuming the quaternion is normalized)
			constexpr VectorN<T, 3> Axis() const noexcept {
				const T coeff{ static_cast<T>(1) / std::sqrt(static_cast<T>(1) - this->W() * this->W()) };
				return VectorN<T, 3>{ coeff * this->X(), coeff * this->Y(), coeff * this->Z() };
			}

			// extract a transformation matrix from quaternion
			template<class LAYOUT = RowMajor<3, 3>> constexpr MatrixNM<T, 3, 3, LAYOUT> DCM() const noexcept {
				const T q0q1{ m_data[0] * m_data[1] },
					q0q2{ m_data[0] * m_data[2] },
					q1q2{ m_data[1] * m_data[2] },
					q3q0{ m_data[3] * m_data[0] },
					q3q1{ m_data[3] * m_data[1] },
					q3q2{ m_data[3] * m_data[2] },
					q0Sqr{ m_data[0] * m_data[0] },
					q1Sqr{ m_data[1] * m_data[1] },
					q2Sqr{ m_data[2] * m_data[2] },
					q3Sqr{ m_data[3] * m_data[3] },
					one{ static_cast<T>(1) },
					two{ static_cast<T>(2) };

				if constexpr (LAYOUT::m_RowMajor) {
					return MatrixNM<T, 3, 3, LAYOUT>{ one - two * (q2Sqr + q1Sqr), two * (q0q1 + q3q2),         two * (q0q2 - q3q1),
									  two * (q0q1 - q3q2),         one - two * (q0Sqr + q2Sqr), two * (q1q2 + q3q0),
									  two * (q0q2 + q3q1),         two * (q1q2 - q3q0),         one - two * (q0Sqr + q1Sqr) };
				}
				else {
					return MatrixNM<T, 3, 3, LAYOUT>{ one - two * (q2Sqr + q1Sqr), two * (q0q1 - q3q2),         two * (q0q2 + q3q1),
									  two * (q0q1 + q3q2),         one - two * (q0Sqr + q2Sqr), two * (q1q2 - q3q0),
									  two * (q0q2 - q3q1),         two * (q1q2 + q3q0),         one - two * (q0Sqr + q1Sqr) };
				}
			}

			// get the conjugate quaternion
			constexpr Quaternion Conjugate() const noexcept {
				return Quaternion{ -m_data[0], -m_data[1], -m_data[2], m_data[3] };
			}

			// get the inverse of quaternion (if invertible - return the natural point)
			constexpr Quaternion Inverse() const noexcept {
				const T lenSqr{ m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2] + m_data[3] * m_data[3] };

				if (!IsZero(lenSqr)) {
					const T lenSqrInv{ static_cast<T>(1) / lenSqr };
					return Quaternion{ -m_data[0] * lenSqrInv, -m_data[1] * lenSqrInv, -m_data[2] * lenSqrInv, m_data[3] * lenSqrInv };
				}
				else {
					return Quaternion(T{});
				}
			}

		// special operations
		public:

			// return quaternion norm (L2; its 'length')
			constexpr T Norm() noexcept {
				T xo_dot{};
				static_for<0, 4>([&](std::size_t i) {
					xo_dot += m_data[i] * m_data[i];
				});

				return std::sqrt(xo_dot);
			}

			// normalize quaternion
			void Normalize() noexcept {
				const T n{ Norm() };
				if (!FloatingPointTrait<T>::Equals(n, T{})) {
					const T nInv{ static_cast<T>(1) / n };
					static_for<0, 4>([&](std::size_t i) {
						m_data[i] *= nInv;
					});
				}
			}

			// rotate a point using a (NORMALIZED) quaternion
			friend void RotatePointUsingQuaternion(VectorN<T, 3>& xio_point, const Quaternion<T>& xi_quat) noexcept {
				const VectorN<T, 3> quatAxis{ xi_quat.X(), xi_quat.Y(), xi_quat.Z() },
									temp{ static_cast<T>(2) * Cross(quatAxis, xio_point) };

				xio_point += temp * xi_quat.W() + Cross(quatAxis, temp);
			}

			friend void RotatePointUsingQuaternion(VectorN<T, 3>& xio_point, Quaternion<T>&& xi_quat) noexcept {
				const VectorN<T, 3> quatAxis{ xi_quat.X(), xi_quat.Y(), xi_quat.Z() },
									temp{ static_cast<T>(2) * Cross(quatAxis, xio_point) };

				xio_point += temp * xi_quat.W() + Cross(quatAxis, temp);
			}

		// numerical + assignment operations
		public:

			// quaternion multiplication
			constexpr Quaternion& operator*=(const Quaternion& xi_quat) {
				Quaternion current{ m_data[0], m_data[1], m_data[2], m_data[3] };

				m_data[0] = current.W() * xi_quat.X() + current.X() * xi_quat.W() + current.Y() * xi_quat.Z() - current.Z() * xi_quat.Y();
				m_data[1] = current.W() * xi_quat.Y() - current.X() * xi_quat.Z() + current.Y() * xi_quat.W() + current.Z() * xi_quat.X();
				m_data[2] = current.W() * xi_quat.Z() + current.X() * xi_quat.Y() - current.Y() * xi_quat.X() + current.Z() * xi_quat.W();
				m_data[3] = current.W() * xi_quat.W() - current.X() * xi_quat.X() - current.Y() * xi_quat.Y() - current.Z() * xi_quat.Z();

				return *this;
			}

		// stream operator
		public:

			// output vector to stream
			friend std::ostream& operator<<(std::ostream& xio_stream, const Quaternion& xi_quat) {
				xio_stream << "{" << xi_quat[0] << ", " << xi_quat[1] << ", " << xi_quat[2] << ", " << xi_quat[3] << "}";

				return xio_stream;
			}

			// read the space-separated components of a vector from a stream
			friend std::istream& operator>>(std::istream& is, const Quaternion& xi_quat) {
				is >> xi_quat[0];
				is >> xi_quat[1];
				is >> xi_quat[2];
				is >> xi_quat[3];

				return is;
			}
	};

	/**
	* operator overloading 
	**/

	template<typename T> constexpr inline Quaternion<T> operator * (const Quaternion<T>& xi_lhs, const Quaternion<T>& xi_rhs) {
		Quaternion<T> xo_quat(xi_lhs);
		xo_quat *= xi_rhs;
		return xo_quat;
	}

	template<typename T> constexpr inline Quaternion<T>& operator * (Quaternion<T>&& xi_lhs, const Quaternion<T>& xi_rhs) {
		xi_lhs *= xi_rhs;
		return xi_lhs;
	}

	template<typename T> constexpr inline Quaternion<T>& operator * (const Quaternion<T>& xi_lhs, Quaternion<T>&& xi_rhs) {
		xi_rhs *= xi_lhs;
		return xi_rhs;
	}

	template<typename T> constexpr inline Quaternion<T>& operator * (Quaternion<T>&& xi_lhs, Quaternion<T>&& xi_rhs) {
		xi_lhs *= xi_rhs;
		return xi_lhs;
	}
};
