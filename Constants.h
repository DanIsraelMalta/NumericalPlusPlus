/**
* This file holds all the
* 1) mathematical constants whose precision is defined per floating point type.
* 2) various 'Numeric' namespace constants/enumerations.
*
* Dan Israel Malta
**/
#pragma once
#ifndef _COSNTANTS_H_
#define _COSNTANTS_H_

#include <cmath>
#include <cinttypes>

#undef PI
#undef HALF_PI
#undef TAU
#undef SQRT2
#undef SQRT3

/*************/
/* Constants */
/*************/
namespace Numeric {

	/***************/
	/* Enumeration */
	/***************/

	// right/left hand side coordinate system type
	enum class Handness : int8_t {
	    left = -1,      // left hand side system
	    right = 1       // right hand side system
	};

	/***********************/
	/* numerical constants */
	/***********************/
	template<class T> struct Constants {
		Constants() = delete;

		// PI
		static constexpr T PI();

		// PI / 2
		static constexpr T HALF_PI();

		// 2 * PI
		static constexpr T TAU();

		// sqrt(2)
		static constexpr T SQRT2();

		// sqrt(3)
		static constexpr T SQRT3();

		// degrees to radians
		static constexpr T DEGREE2RADIAN();

		// radians to degrees
		static constexpr T RADIAN2DEGREE();
	};

	// single precision constants
	template<> struct Constants<float> {
		Constants() = delete;

		static constexpr float PI() { return 3.141592654f; }
		static constexpr float HALF_PI() { return 1.570796327f; }
		static constexpr float TAU() { return 6.283185307f; }
		static constexpr float SQRT2() { return 1.414213562f; }
		static constexpr float SQRT3() { return 1.732050808f; }
		static constexpr float DEGREE2RADIAN() { return 0.017453292f; }
		static constexpr float RADIAN2DEGREE() { return 57.29577951f; }
	};

	// double precision constants
	template<> struct Constants<double> {
		Constants() = delete;

		static constexpr double PI() { return 3.141592653589793; }
		static constexpr double HALF_PI() { return 1.570796326794897; }
		static constexpr double TAU() { return 6.283185307179586; }
		static constexpr double SQRT2() { return 1.414213562373095; }
		static constexpr double SQRT3() { return 1.732050807568877; }
		static constexpr double DEGREE2RADIAN() { return 0.017453292519943; }
		static constexpr double RADIAN2DEGREE() { return 57.29577951308232; }
	};

	// long double precision constants
	template<> struct Constants<long double> {
		Constants() = delete;

		static constexpr long double PI() { return 3.1415926535897932384626433832795; }
		static constexpr long double HALF_PI() { return 1.5707963267948966192313216916398; }
		static constexpr long double TAU() { return 6.283185307179586476925286766559; }
		static constexpr long double SQRT2() { return 1.4142135623730950488016887242097; }
		static constexpr long double SQRT3() { return 1.7320508075688772935274463415059; }
		static constexpr long double DEGREE2RADIAN() { return 0.01745329251994329576923690768489; }
		static constexpr long double RADIAN2DEGREE() { return 57.295779513082320876798154814105; }
	};
};

#endif //_COSNTANTS_H_
