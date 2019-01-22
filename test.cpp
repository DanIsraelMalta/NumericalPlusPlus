#include <iostream>
#include <array>

#include "static_for.h"
#include "VectorN.h"
#include "MatrixNM.h"
#include "Interval.h"
#include "Quaternion.h"
#include "Interval.h"
#include "Simplex.h"
#include "3D.h"
#include <future>
#include "IntegerSequence.h"


int main() {
	// --- VectorN tests ---
	{
		using namespace Numeric;

		VectorN<float, 3> a{ {0.0f, 1.0f, 2.0f} },
			b{ 0.0f, 1.0f, 2.0f },
			c(0.0f),
			cp(0, 1, 2),
			cr(rising);
		c += VectorN<float, 3>{ 0.0f, 1.0f, 2.0f };
		c -= 1.0f;

		std::cout << "c = " << c << "\n";

		const VectorN<float, 3> aa(2.0f);

		c += 2.0f;
		c *= aa;
		std::cout << "c = " << c << "\n";

		c += b - 2.0f * cp;
		std::cout << "c = " << c << "\n";

		c.Reverse();
		std::cout << "c = " << c << "\n";

		c.RotateLeft(1);
		std::cout << "c = " << c << "\n";

		c.Sort([](const auto &a, const auto &b) -> bool { return a < b; });
		std::cout << "c = " << c << "\n";

		float v = c.FindIf([](const auto& a) -> bool { return ((a > 1.0f) && (a <= 2.0f)); });	// v=2.0f

		const float m = c.MinElement();

		const float dd = EuclideanDistanceBetween(a, c);

		bool isc = isVectorN(c);	// true

		a.X() = -2.5f;
		std::cout << "a = " << a << "\n";

		// multiplication using 'move' semantics
		a = VectorN<float, 3>{1.0f, 22.0f, 3.0f} * VectorN<float, 3>{1.0f, 22.0f, 3.0f};

		// sin function is used with 'move' semantics
		a = sin(b * c);		

		auto lin = Lerp(VectorN<float, 3>(0.0f), VectorN<float, 3>{1.0f, 22.0f, 3.0f}, 0.5f);

		auto dotCopy = Dot(a, b);
		auto dotMove = Dot(VectorN<float, 3>(0.0f), VectorN<float, 3>{1.0f, 22.0f, 3.0f});

		auto distCopy = EuclideanDistanceBetween(a, b);
		auto distMove = EuclideanDistanceBetween(VectorN<float, 3>(0.0f), VectorN<float, 3>{1.0f, 22.0f, 3.0f});

		VectorN<int, 4> c0{{0, 1, 2, 3}},
						c2{{-1, 3, 3, 3}};
		auto cmin = Min(c0, c2);
		auto cmax = Max(c0, c2, cmin);
		int s3 = Sum(c0, cmin, cmax);
		int m3 = Mean(c0, cmin, cmax);
		int p3 = Prod(c0, cmin, cmax);
		std::cout << "sum = " << s3 << "\n";
		std::cout << "mean = " << m3 << "\n";

		VectorN<float, 4> seq(sequence, 3.0f, -2.0f);
		std::cout << "seq = " << seq << "\n";

		VectorN<float, 10> fib(fibonacci, 0.0f, 1.0f);
		std::cout << "fib = " << fib << "\n";

		// swizzle
		VectorN<float, 4> sw{ {0.0f, 1.0f, 2.0f, 3.0f} };
		VectorN<float, 3> yyz{ sw.swizzle<1, 1, 2>() },
						  zxy{ sw.swizzle<2, 0, 1>() };

		// *** asynchronous calculations ***
		std::future<VectorN<float, 3>> cf = std::async(std::launch::async, [&c]() { return sqrt(c); });
		std::cout << "cf = " << cf.get() << "\n";
	}

	// --- MatrixNM test ---
	{
		using namespace Numeric;

		MatrixNM<float, 2, 3> a{ { 1.0f, 2.0f, 3.0f,
								   4.0f, 5.0f, 6.0f } };
		std::cout << "a = " << a << "\n";

		MatrixNM<float, 3, 3> ar{ 1.0f, 2.0f, 3.0f,
								  4.0f, 5.0f, 6.0f,
								  7.0f, 8.0f, 9.0f };
		std::cout << "ar = " << ar << "\n";
		std::cout << "ar(1,2) = " << ar(1,2) << "\n";

		MatrixNM<float, 3, 3, ColumnMajor<3, 3>> ac{ 1.0f, 2.0f, 3.0f,
													 4.0f, 5.0f, 6.0f,
													 7.0f, 8.0f, 9.0f };
		std::cout << "ac = " << ac << "\n";
		std::cout << "ac(1,2) = " << ac(1,2) << "\n";

		ac(1, 2) = 1.2f;
		std::cout << "ac(1,2) = " << ac(1,2) << "\n";

		auto aca = -1.0f * abs(ac);
		std::cout << "aca = " << aca << "\n";

		MatrixNM<float, 3, 3> acLow(ac(LowerTriangular));
		std::cout << "ac (lower triangular) = " << acLow << "\n";

		MatrixNM<float, 3, 3> acUp(ac(UpperTriangular));
		std::cout << "ac (upper triangular) = " << acUp << "\n";

		VectorN<float, 3> row0(ar(Row(0))),
						  col0(ar(Column(0))),
						  diag(ar(Diagonal));
		std::cout << "ar (row 0) = " << row0 << "\n";
		std::cout << "ar (column 0) = " << col0 << "\n";
		std::cout << "ar (diag) = " << diag << "\n";

		MatrixNM<float, 4, 4> e(Eye),
			outerproduct(Outerproduct, VectorN<float, 4>{1.0f, 22.0f, 3.0f, 4.0f}, VectorN<float, 4>{1.0f, 22.0f, 3.0f, 4.0f}),
			rowWise(RowWise, VectorN<float, 4>{1.0f, 22.0f, 3.0f, 4.0f}),
			ColumnWise(ColumnWise, VectorN<float, 4>{1.0f, 22.0f, 3.0f, 4.0f}),
			diagonalWise(DiagonalWise, VectorN<float, 4>{1.0f, 22.0f, 3.0f, 4.0f}),
			vanDerMonde(VanDerMonde, VectorN<float, 4>{1.0f, 2.0f, 3.0f, 4.0f}),
			toeplitz(Toeplitz, VectorN<float, 4>{1.0f, 2.0f, 3.0f, 4.0f});
		MatrixNM<float, 3, 3> axisAngle(AxisAngle, VectorN<float, 3>{1.0f, 0.0f, 0.0f}, 1.0f, 0.0f),
							  dcm = e.DCM();
		MatrixNM<float, 5, 5> vanDerMonde5(VanDerMonde, VectorN<float, 5>{ 1.0f, 1.5f, 2.0f, 2.5f, 3.0f });

		std::cout << "outerproduct = " << outerproduct << "\n";
		std::cout << "rowWise = " << rowWise << "\n";
		std::cout << "ColumnWise = " << ColumnWise << "\n";
		std::cout << "diagonalWise = " << diagonalWise << "\n";
		std::cout << "vanDerMonde = " << vanDerMonde << "\n";
		std::cout << "toeplitz = " << toeplitz << "\n";
		std::cout << "axisAngle = " << axisAngle << "\n";
		std::cout << "dcm = " << dcm << "\n";
		std::cout << "vanDerMonde5 = " << vanDerMonde5 << "\n";

		ac = { 5.0f, 2.0f, 3.0f,
			   4.0f, 5.0f, 6.0f,
			   7.0f, 8.0f, 9.0f };
		std::cout << "ac = " << ac << "\n";

		MatrixNM<float, 2, 2> e2( e.GetRegion<0, 2, 0, 2>());
		std::cout << "e2 = " << e2 << "\n";

		VectorN<float, 2> e2c{ e2.X() };
		std::cout << "e2c = " << e2c << "\n";

		auto b = axisAngle + 2.0f * abs(dcm);
		std::cout << "b = " << b << "\n";

		auto minor22 = toeplitz.Minor<1, 1>();
		std::cout << "minor22 = " << minor22 << "\n";

		auto ar12 = ar.Minor<1, 2>();
		std::cout << "ac12 = " << ar12 << "\n";

		// determinant

		auto det1 = ar.Determinant();
		std::cout << "|ar| = " << det1 << " (0)\n";

		MatrixNM<double, 5> luDecomp{ 12, 13,  14, 15,  16,
									   16, 15,  14, 13,  12,
									   38, 75, -15,  5,  11,
									   92, 32,  27,  5, -31,
									   19,  8,  75,  5,  55 };
		double det3 = luDecomp.Determinant();
		std::cout << "det3 = " << det3 << " (7182784)\n";

		MatrixNM<float, 4> mat4{ 12.0f, 13.0f,  14.0f, 15.0f,
								 16.0f, 15.0f,  14.0f, 13.0f,
								 38.0f, 75.0f, -15.0f,  5.0f,
								 92.0f, 32.0f,  27.0f,  5.0f };

		float det4 = mat4.Determinant();
		std::cout << "det4 = " << det4 << " (108948)\n";

		MatrixNM<int32_t, 3> mat3{ 1, 2, 3,
								 0, 4, 5,
								 0, 0, 6 };
		int32_t det5 = mat3.Determinant();
		std::cout << "det5 = " << det5 << " (24)\n";


		// LU decomposition
		MatrixNM<double, 5> lowerUpper;
		VectorN<std::size_t, 5> piv;
		int32_t _sign;
		luDecomp.LU(lowerUpper, piv, _sign);

		assert(static_cast<int64_t>(lowerUpper(0, 0)) == 92);
		assert(static_cast<int64_t>(lowerUpper(0, 1)) == 32);
		assert(static_cast<int64_t>(lowerUpper(0, 2)) == 27);
		assert(static_cast<int64_t>(lowerUpper(0, 3)) == 5);
		assert(static_cast<int64_t>(lowerUpper(0, 4)) == -31);
		assert(static_cast<int64_t>(lowerUpper(1, 0) * 1000000) == 413043);
		assert(static_cast<int64_t>(lowerUpper(1, 1) * 1000000) == 61782608);
		assert(static_cast<int64_t>(lowerUpper(1, 2) * 1000000) == -26152173);
		assert(static_cast<int64_t>(lowerUpper(1, 3) * 1000000) == 2934782);
		assert(static_cast<int64_t>(lowerUpper(1, 4) * 1000000) == 23804347);
		assert(static_cast<int64_t>(lowerUpper(2, 0) * 1000000) == 206521);
		assert(static_cast<int64_t>(lowerUpper(2, 1) * 1000000) == 22519);
		assert(static_cast<int64_t>(lowerUpper(2, 2) * 1000000) == 70012843);
		assert(static_cast<int64_t>(lowerUpper(2, 3) * 1000000) == 3901301);
		assert(static_cast<int64_t>(lowerUpper(2, 4) * 1000000) == 60866115);
		assert(static_cast<int64_t>(lowerUpper(3, 0) * 1000000) == 130434);
		assert(static_cast<int64_t>(lowerUpper(3, 1) * 1000000) == 142857);
		assert(static_cast<int64_t>(lowerUpper(3, 2) * 1000000) == 203023);
		assert(static_cast<int64_t>(lowerUpper(3, 3) * 1000000) == 13136513);
		assert(static_cast<int64_t>(lowerUpper(3, 4) * 1000000) == 4285576);
		assert(static_cast<int64_t>(lowerUpper(4, 0) * 1000000) == 173913);
		assert(static_cast<int64_t>(lowerUpper(4, 1) * 1000000) == 152709);
		assert(static_cast<int64_t>(lowerUpper(4, 2) * 1000000) == 189937);
		assert(static_cast<int64_t>(lowerUpper(4, 3) * 1000000) == 832889);
		assert(static_cast<int64_t>(lowerUpper(4, 4) * 1000000) == -1373981);

		// solve a cubic linear system
		VectorN<double, 5> bc{ 1.0f, 2.0f, 3.0f, 4.0f, 5.0f }, x;
		luDecomp.SolveSquare(bc, x);
		assert(static_cast<int64_t>(x[0] * 10000) == -4733);
		assert(static_cast<int64_t>(x[1] * 10000) == 4612);
		assert(static_cast<int64_t>(x[2] * 10000) == 5411);
		assert(static_cast<int64_t>(x[3] * 10000) == 1415);
		assert(static_cast<int64_t>(x[4] * 10000) == -5634);
		
		// QR decomposition
		MatrixNM<double, 3> QR33({ 12, -51,   4,
							 	   6, 167, -68,
								  -4,  24, -41 }),
							q33, r33;
		QR33.QR(q33, r33);
		assert(static_cast<int64_t>(q33(0, 0) * 10000) == 8571);
		assert(static_cast<int64_t>(q33(0, 1) * 10000) == -3942);
		assert(static_cast<int64_t>(q33(0, 2) * 10000) == 3314);
		assert(static_cast<int64_t>(q33(1, 0) * 10000) == 4285);
		assert(static_cast<int64_t>(q33(1, 1) * 10000) == 9028);
		assert(static_cast<int64_t>(q33(1, 2) * 10000) == -342);
		assert(static_cast<int64_t>(q33(2, 0) * 10000) == -2857);
		assert(static_cast<int64_t>(q33(2, 1) * 10000) == 1714);
		assert(static_cast<int64_t>(q33(2, 2) * 10000) == 9428);

		assert(static_cast<int64_t>(std::round(r33(0, 0))) == 14);
		assert(static_cast<int64_t>(std::round(r33(0, 1))) == 21);
		assert(static_cast<int64_t>(std::round(r33(0, 2))) == -14);
		assert(static_cast<int64_t>(std::round(r33(1, 0))) == 0);
		assert(static_cast<int64_t>(std::round(r33(1, 1))) == 175);
		assert(static_cast<int64_t>(std::round(r33(1, 2))) == -70);
		assert(static_cast<int64_t>(std::round(r33(2, 0))) == 0);
		assert(static_cast<int64_t>(std::round(r33(2, 1))) == 0);
		assert(static_cast<int64_t>(std::round(r33(2, 2))) == -35);

		MatrixNM<double, 5, 3> QR53({ 12, -51,   4,
									  6, 167, -68,
									 -4,  24, -41,
									 -1,   1,   0,
									 2,   0,   3 });
		MatrixNM<double, 3, 5> q53;
		MatrixNM<double, 3, 3> r53;
		QR53.QR(q53, r53);

		assert(static_cast<int64_t>(q53(0, 0) * 1000000) == 846414);
		assert(static_cast<int64_t>(q53(0, 1) * 1000000) == 423207);
		assert(static_cast<int64_t>(q53(0, 2) * 1000000) == -282138);
		assert(static_cast<int64_t>(q53(0, 3) * 1000000) == -70534);
		assert(static_cast<int64_t>(q53(0, 4) * 1000000) == 141069);
		assert(static_cast<int64_t>(q53(1, 0) * 1000000) == -391290);
		assert(static_cast<int64_t>(q53(1, 1) * 1000000) == 904087);
		assert(static_cast<int64_t>(q53(1, 2) * 1000000) == 170420);
		assert(static_cast<int64_t>(q53(1, 3) * 1000000) == 14040);
		assert(static_cast<int64_t>(q53(1, 4) * 1000000) == -16655);
		assert(static_cast<int64_t>(q53(2, 0) * 1000000) == -343124);
		assert(static_cast<int64_t>(q53(2, 1) * 1000000) == 29270);
		assert(static_cast<int64_t>(q53(2, 2) * 1000000) == -932855);
		assert(static_cast<int64_t>(q53(2, 3) * 1000000) == 1099);
		assert(static_cast<int64_t>(q53(2, 4) * 1000000) == 105771);

		assert(static_cast<int64_t>(r53(0, 0) * 1000) == 14177);
		assert(static_cast<int64_t>(r53(0, 1) * 1000) == 20666);
		assert(static_cast<int64_t>(r53(0, 2) * 1000) == -13401);
		assert(static_cast<int64_t>(r53(1, 0)) == 0);
		assert(static_cast<int64_t>(r53(1, 1) * 1000) == 175042);
		assert(static_cast<int64_t>(r53(1, 2) * 1000) == -70080);
		assert(static_cast<int64_t>(r53(2, 0)) == 0);
		assert(static_cast<int64_t>(r53(2, 1)) == 0);
		assert(static_cast<int64_t>(r53(2, 2) * 1000) == 35201);

		// rectangular linear system
		MatrixNM<double, 3, 2> A32({ 3.0, -2.0,
									 0.0, 3.0,
									 4.0, 4.0 });
		VectorN<double, 3> B3{ 3, 5, 4 };
		VectorN<double, 2> X2;
		A32.SolveRectangualr(B3, X2);

		assert(static_cast<int64_t>(X2[0] * 100) == 76);
		assert(static_cast<int64_t>(X2[1] * 10) == 6);

		// SVD test (cubic)
		MatrixNM<double, 5> usvd5, uwsvd5, vsvd5;
		VectorN<double, 5> wsvd5;
		bool flag = luDecomp.SVD(uwsvd5, wsvd5, vsvd5);
		assert(flag == true);

		assert(static_cast<int64_t>(std::floor(std::sqrt(wsvd5[0]) * 10000)) == 1301140);
		assert(static_cast<int64_t>(std::floor(std::sqrt(wsvd5[1]) * 10000)) == 931333);
		assert(static_cast<int64_t>(std::floor(std::sqrt(wsvd5[2]) * 10000)) == 621900);
		assert(static_cast<int64_t>(std::floor(std::sqrt(wsvd5[3]) * 10000)) == 158604);
		assert(static_cast<int64_t>(std::floor(std::sqrt(wsvd5[4]) * 10000)) == 6009);

		vsvd5 *= 10000;
		vsvd5 = floor(abs(vsvd5));

		assert(static_cast<int64_t>(vsvd5(0, 0)) == 7578);
		assert(static_cast<int64_t>(vsvd5(0, 1)) == 2182);
		assert(static_cast<int64_t>(vsvd5(0, 2)) == 3706);
		assert(static_cast<int64_t>(vsvd5(0, 3)) == 381);
		assert(static_cast<int64_t>(vsvd5(0, 4)) == 4890);
		assert(static_cast<int64_t>(vsvd5(1, 0)) == 5309);
		assert(static_cast<int64_t>(vsvd5(1, 1)) == 2600);
		assert(static_cast<int64_t>(vsvd5(1, 2)) == 6708);
		assert(static_cast<int64_t>(vsvd5(1, 3)) == 1601);
		assert(static_cast<int64_t>(vsvd5(1, 4)) == 4180);
		assert(static_cast<int64_t>(vsvd5(2, 0)) == 3572);
		assert(static_cast<int64_t>(vsvd5(2, 1)) == 7084);
		assert(static_cast<int64_t>(vsvd5(2, 2)) == 3421);
		assert(static_cast<int64_t>(vsvd5(2, 3)) == 1289);
		assert(static_cast<int64_t>(vsvd5(2, 4)) == 4865);
		assert(static_cast<int64_t>(vsvd5(3, 0)) == 1037);
		assert(static_cast<int64_t>(vsvd5(3, 1)) == 535);
		assert(static_cast<int64_t>(vsvd5(3, 2)) == 741);
		assert(static_cast<int64_t>(vsvd5(3, 3)) == 9778);
		assert(static_cast<int64_t>(vsvd5(3, 4)) == 1570);
		assert(static_cast<int64_t>(vsvd5(4, 0)) == 727);
		assert(static_cast<int64_t>(vsvd5(4, 1)) == 6163);
		assert(static_cast<int64_t>(vsvd5(4, 2)) == 5385);
		assert(static_cast<int64_t>(vsvd5(4, 3)) == 91);
		assert(static_cast<int64_t>(vsvd5(4, 4)) == 5697);

		for (std::size_t i{}; i < 5; ++i) {
			usvd5(0, i) = uwsvd5(0, i) / std::sqrt(wsvd5[i]);
			usvd5(1, i) = uwsvd5(1, i) / std::sqrt(wsvd5[i]);
			usvd5(2, i) = uwsvd5(2, i) / std::sqrt(wsvd5[i]);
			usvd5(3, i) = uwsvd5(3, i) / std::sqrt(wsvd5[i]);
			usvd5(4, i) = uwsvd5(4, i) / std::sqrt(wsvd5[i]);
		}

		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(0, 0)) * 10000)) == 1822);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(0, 1)) * 10000)) == 1566);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(0, 2)) * 10000)) == 1481);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(0, 3)) * 10000)) == 7179);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(0, 4)) * 10000)) == 6363);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(1, 0)) * 10000)) == 2099);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(1, 1)) * 10000)) == 1140);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(1, 2)) * 10000)) == 1088);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(1, 3)) * 10000)) == 5817);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(1, 4)) * 10000)) == 7698);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(2, 0)) * 10000)) == 4963);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(2, 1)) * 10000)) == 3368);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(2, 2)) * 10000)) == 7663);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(2, 3)) * 10000)) == 2290);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(2, 4)) * 10000)) == 207);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(3, 0)) * 10000)) == 7272);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(3, 1)) * 10000)) == 3018);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(3, 2)) * 10000)) == 6141);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(3, 3)) * 10000)) == 307);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(3, 4)) * 10000)) == 435);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(4, 0)) * 10000)) == 3839);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(4, 1)) * 10000)) == 8705);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(4, 2)) * 10000)) == 426);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(4, 3)) * 10000)) == 3046);
		assert(static_cast<int64_t>(std::floor(std::abs(usvd5(4, 4)) * 10000)) == 94);

		// another SVD test (cubic)
		MatrixNM<double, 6> svd6({ 27, 92, 32, 5, -31, 12,
								   14, 12, 13, -15, 16, -8,
								   14, 16, 15, 13, 12 ,11,
								  -15, 38, 75, 5, 11, 5,
								   75, 19, -8, 5, 55, 55,
								   11, -9, 3, 17, 29, 35 }),
			u6, uw6, v6;
		VectorN<double, 6> w6;
		flag = svd6.SVD(uw6, w6, v6);
		assert(flag == true);

		for (std::size_t i = 0; i < 5; ++i) {
			u6(0, i) = uw6(0, i) / std::sqrt(w6[i]);
			u6(1, i) = uw6(1, i) / std::sqrt(w6[i]);
			u6(2, i) = uw6(2, i) / std::sqrt(w6[i]);
			u6(3, i) = uw6(3, i) / std::sqrt(w6[i]);
			u6(4, i) = uw6(4, i) / std::sqrt(w6[i]);
			u6(5, i) = uw6(5, i) / std::sqrt(w6[i]);
		}

		VectorN<double, 6> w66;
		w66 = sqrt(w6);
		w66 *= 10000;
		w66 = floor(w66);

		assert(static_cast<int64_t>(w66[0]) == 1323184);
		assert(static_cast<int64_t>(w66[1]) == 1127663);
		assert(static_cast<int64_t>(w66[2]) == 650989);
		assert(static_cast<int64_t>(w66[3]) == 334712);
		assert(static_cast<int64_t>(w66[4]) == 88533);
		assert(static_cast<int64_t>(w66[5]) == 26458);

		v6 *= 10000;
		v6 = floor(abs(v6));

		assert(static_cast<int64_t>(v6(0, 0)) == 4585);
		assert(static_cast<int64_t>(v6(0, 1)) == 4642);
		assert(static_cast<int64_t>(v6(0, 2)) == 3618);
		assert(static_cast<int64_t>(v6(0, 3)) == 3085);
		assert(static_cast<int64_t>(v6(0, 4)) == 2823);
		assert(static_cast<int64_t>(v6(0, 5)) == 5179);
		assert(static_cast<int64_t>(v6(1, 0)) == 6821);
		assert(static_cast<int64_t>(v6(1, 1)) == 4045);
		assert(static_cast<int64_t>(v6(1, 2)) == 3501);
		assert(static_cast<int64_t>(v6(1, 3)) == 255);
		assert(static_cast<int64_t>(v6(1, 4)) == 63);
		assert(static_cast<int64_t>(v6(1, 5)) == 4977);
		assert(static_cast<int64_t>(v6(2, 0)) == 3942);
		assert(static_cast<int64_t>(v6(2, 1)) == 4443);
		assert(static_cast<int64_t>(v6(2, 2)) == 6639);
		assert(static_cast<int64_t>(v6(2, 3)) == 1118);
		assert(static_cast<int64_t>(v6(2, 4)) == 118);
		assert(static_cast<int64_t>(v6(2, 5)) == 4400);
		assert(static_cast<int64_t>(v6(3, 0)) == 867);
		assert(static_cast<int64_t>(v6(3, 1)) == 467);
		assert(static_cast<int64_t>(v6(3, 2)) == 1050);
		assert(static_cast<int64_t>(v6(3, 3)) == 6751);
		assert(static_cast<int64_t>(v6(3, 4)) == 7225);
		assert(static_cast<int64_t>(v6(3, 5)) == 369);
		assert(static_cast<int64_t>(v6(4, 0)) == 1830);
		assert(static_cast<int64_t>(v6(4, 1)) == 5194);
		assert(static_cast<int64_t>(v6(4, 2)) == 5331);
		assert(static_cast<int64_t>(v6(4, 3)) == 3290);
		assert(static_cast<int64_t>(v6(4, 4)) == 1472);
		assert(static_cast<int64_t>(v6(4, 5)) == 5314);
		assert(static_cast<int64_t>(v6(5, 0)) == 3577);
		assert(static_cast<int64_t>(v6(5, 1)) == 3889);
		assert(static_cast<int64_t>(v6(5, 2)) == 1015);
		assert(static_cast<int64_t>(v6(5, 3)) == 5722);
		assert(static_cast<int64_t>(v6(5, 4)) == 6135);
		assert(static_cast<int64_t>(v6(5, 5)) == 810);

		u6 *= 10000;
		u6 = floor(abs(u6));

		assert(static_cast<int64_t>(u6(0, 0)) == 6560);
		assert(static_cast<int64_t>(u6(0, 1)) == 4443);
		assert(static_cast<int64_t>(u6(0, 2)) == 5457);
		assert(static_cast<int64_t>(u6(0, 3)) == 1846);
		assert(static_cast<int64_t>(u6(0, 4)) == 1008);
		assert(static_cast<int64_t>(u6(0, 5)) == 0);
		assert(static_cast<int64_t>(u6(1, 0)) == 1397);
		assert(static_cast<int64_t>(u6(1, 1)) == 32);
		assert(static_cast<int64_t>(u6(1, 2)) == 845);
		assert(static_cast<int64_t>(u6(1, 3)) == 7783);
		assert(static_cast<int64_t>(u6(1, 4)) == 516);
		assert(static_cast<int64_t>(u6(1, 5)) == 0);
		assert(static_cast<int64_t>(u6(2, 0)) == 2305);
		assert(static_cast<int64_t>(u6(2, 1)) == 397);
		assert(static_cast<int64_t>(u6(2, 2)) == 1255);
		assert(static_cast<int64_t>(u6(2, 3)) == 1409);
		assert(static_cast<int64_t>(u6(2, 4)) == 9533);
		assert(static_cast<int64_t>(u6(2, 5)) == 0);
		assert(static_cast<int64_t>(u6(3, 0)) == 3993);
		assert(static_cast<int64_t>(u6(3, 1)) == 4236);
		assert(static_cast<int64_t>(u6(3, 2)) == 7498);
		assert(static_cast<int64_t>(u6(3, 3)) == 631);
		assert(static_cast<int64_t>(u6(3, 4)) == 1603);
		assert(static_cast<int64_t>(u6(3, 5)) == 0);
		assert(static_cast<int64_t>(u6(4, 0)) == 5621);
		assert(static_cast<int64_t>(u6(4, 1)) == 7172);
		assert(static_cast<int64_t>(u6(4, 2)) == 563);
		assert(static_cast<int64_t>(u6(4, 3)) == 1787);
		assert(static_cast<int64_t>(u6(4, 4)) == 1214);
		assert(static_cast<int64_t>(u6(4, 5)) == 0);
		assert(static_cast<int64_t>(u6(5, 0)) == 1465);
		assert(static_cast<int64_t>(u6(5, 1)) == 3271);
		assert(static_cast<int64_t>(u6(5, 2)) == 3374);
		assert(static_cast<int64_t>(u6(5, 3)) == 5516);
		assert(static_cast<int64_t>(u6(5, 4)) == 1946);
		assert(static_cast<int64_t>(u6(5, 5)) == 0);

		// SVD test (rectangular)
		MatrixNM<double, 7, 5> svd75({ 1, 2, 3, 4, 5,
									  6, 7, 8, 9, 1,
									  1, 2, 3, 4, 5,
									  6, 7, 8, 9, 1,
									  2, 3, 4, 5, 6,
									  7, 8, 9, 1, 2,
									  3, 4, 5, 6, 7 }),
								uw75;
		MatrixNM<double, 5> v75;
		VectorN<double, 5> w75;

		svd75.SVD(uw75, w75, v75);

		w75 = sqrt(w75);
		w75 *= 10000;
		w75 = floor(w75);

		assert(static_cast<int64_t>(w75[0]) == 293662);
		assert(static_cast<int64_t>(w75[1]) == 96911);
		assert(static_cast<int64_t>(w75[2]) == 62797);
		assert(static_cast<int64_t>(w75[3]) == 5207);
		assert(static_cast<int64_t>(w75[4]) == 0);

		v75 *= 10000;
		v75 = floor(abs(v75));

		assert(static_cast<int64_t>(v75(0, 0)) == 3784);
		assert(static_cast<int64_t>(v75(0, 1)) == 3505);
		assert(static_cast<int64_t>(v75(0, 2)) == 1402);
		assert(static_cast<int64_t>(v75(0, 3)) == 7399);
		assert(static_cast<int64_t>(v75(0, 4)) == 4082);
		assert(static_cast<int64_t>(v75(1, 0)) == 4644);
		assert(static_cast<int64_t>(v75(1, 1)) == 2845);
		assert(static_cast<int64_t>(v75(1, 2)) == 1866);
		assert(static_cast<int64_t>(v75(1, 3)) == 427);
		assert(static_cast<int64_t>(v75(1, 4)) == 8164);
		assert(static_cast<int64_t>(v75(2, 0)) == 5504);
		assert(static_cast<int64_t>(v75(2, 1)) == 2184);
		assert(static_cast<int64_t>(v75(2, 2)) == 2330);
		assert(static_cast<int64_t>(v75(2, 3)) == 6544);
		assert(static_cast<int64_t>(v75(2, 4)) == 4082);
		assert(static_cast<int64_t>(v75(3, 0)) == 5070);
		assert(static_cast<int64_t>(v75(3, 1)) == 3049);
		assert(static_cast<int64_t>(v75(3, 2)) == 8053);
		assert(static_cast<int64_t>(v75(3, 3)) == 378);
		assert(static_cast<int64_t>(v75(3, 4)) == 0);
		assert(static_cast<int64_t>(v75(4, 0)) == 2846);
		assert(static_cast<int64_t>(v75(4, 1)) == 8096);
		assert(static_cast<int64_t>(v75(4, 2)) == 4925);
		assert(static_cast<int64_t>(v75(4, 3)) == 1446);
		assert(static_cast<int64_t>(v75(4, 4)) == 0);


		uw75 *= 10000;
		uw75 = floor(abs(uw75));

		assert(static_cast<int64_t>(uw75(0, 0)) == 64098);
		assert(static_cast<int64_t>(uw75(0, 1)) == 36928);
		assert(static_cast<int64_t>(uw75(0, 2)) == 4545);
		assert(static_cast<int64_t>(uw75(0, 3)) == 2637);
		assert(static_cast<int64_t>(uw75(0, 4)) == 0);
		assert(static_cast<int64_t>(uw75(1, 0)) == 147731);
		assert(static_cast<int64_t>(uw75(1, 1)) == 22884);
		assert(static_cast<int64_t>(uw75(1, 2)) == 27417);
		assert(static_cast<int64_t>(uw75(1, 3)) == 121);
		assert(static_cast<int64_t>(uw75(1, 4)) == 0);
		assert(static_cast<int64_t>(uw75(2, 0)) == 64098);
		assert(static_cast<int64_t>(uw75(2, 1)) == 36928);
		assert(static_cast<int64_t>(uw75(2, 2)) == 4545);
		assert(static_cast<int64_t>(uw75(2, 3)) == 2637);
		assert(static_cast<int64_t>(uw75(2, 4)) == 0);
		assert(static_cast<int64_t>(uw75(3, 0)) == 147731);
		assert(static_cast<int64_t>(uw75(3, 1)) == 22884);
		assert(static_cast<int64_t>(uw75(3, 2)) == 27417);
		assert(static_cast<int64_t>(uw75(3, 3)) == 121);
		assert(static_cast<int64_t>(uw75(3, 4)) == 0);
		assert(static_cast<int64_t>(uw75(4, 0)) == 85948);
		assert(static_cast<int64_t>(uw75(4, 1)) == 39538);
		assert(static_cast<int64_t>(uw75(4, 2)) == 7018);
		assert(static_cast<int64_t>(uw75(4, 3)) == 469);
		assert(static_cast<int64_t>(uw75(4, 4)) == 0);
		assert(static_cast<int64_t>(uw75(5, 0)) == 123949);
		assert(static_cast<int64_t>(uw75(5, 1)) == 47717);
		assert(static_cast<int64_t>(uw75(5, 2)) == 47533);
		assert(static_cast<int64_t>(uw75(5, 3)) == 418);
		assert(static_cast<int64_t>(uw75(5, 4)) == 0);
		assert(static_cast<int64_t>(uw75(6, 0)) == 107798);
		assert(static_cast<int64_t>(uw75(6, 1)) == 42149);
		assert(static_cast<int64_t>(uw75(6, 2)) == 9492);
		assert(static_cast<int64_t>(uw75(6, 3)) == 3575);
		assert(static_cast<int64_t>(uw75(6, 4)) == 0);

		// Cholesky decomposition
		MatrixNM<double, 3> spd{{ 4,  12, -16,
								  12,  37, -43,
								 -16, -43,  98 } },
							chol;
		spd.Cholesky(chol);

		assert(chol.IsLowerTriangular() == true);
		assert(static_cast<int64_t>(chol(0, 0)) == 2);
		assert(static_cast<int64_t>(chol(0, 1)) == 0);
		assert(static_cast<int64_t>(chol(0, 2)) == 0);
		assert(static_cast<int64_t>(chol(1, 0)) == 6);
		assert(static_cast<int64_t>(chol(1, 1)) == 1);
		assert(static_cast<int64_t>(chol(1, 2)) == 0);
		assert(static_cast<int64_t>(chol(2, 0)) == -8);
		assert(static_cast<int64_t>(chol(2, 1)) == 5);
		assert(static_cast<int64_t>(chol(2, 2)) == 3);

		// solve a linear system using Cholesky
		VectorN<double, 3> bcc{ 6, 3, 1 }, xc;
		spd.SolveCholesky(bcc, xc);

		assert(static_cast<int64_t>(xc[0] * 100) == 25761);
		assert(static_cast<int64_t>(xc[1] * 100) == -7055);
		assert(static_cast<int64_t>(xc[2] * 100) == 1111);

		// is???

		flag = ar.IsDiagonal();
		assert(flag == false);
		flag = diagonalWise.IsDiagonal();
		assert(flag == true);

		flag = ar.IsLowerTriangular();
		assert(flag == false);
		flag = acLow.IsLowerTriangular();

		flag = acLow.IsSymmetric();
		assert(flag == false);
		flag = acLow.IsDiagonal();
		assert(flag == false);

		flag = acLow.IsUpperTriangular();
		assert(flag == false);

		flag = acUp.IsUpperTriangular();
		assert(flag == true);
		flag = acUp.IsSymmetric();
		assert(flag == false);
		flag = acUp.IsDiagonal();
		assert(flag == false);

		MatrixNM<int32_t, 5> mn({ 1, 2, 3, 4, 5,
								  2, 1, 2, 3, 4,
								  3, 2, 1, 2, 3,
								  4, 3, 2, 1, 2,
								  5, 4, 3, 2, 1 }),
							 dg({ 1, 0, 0, 0, 0,
							 	  0, 1, 0, 0, 0,
								  0, 0, 1, 0, 0,
								  0, 0, 0, 1, 0,
								  0, 0, 0, 0, 1 });
		flag = mn.IsSymmetric();
		assert(flag == true);

		flag = dg.IsDiagonal();
		assert(flag == true);

		MatrixNM<int32_t, 4> tri({ 0, 1, 0, 0,
							 	   0, 0, 0, 1,
								   1, 0, 0, 0,
								   0, 0, 1, 0 });
		flag = dg.IsPermutation();
		assert(flag == true);
		flag = tri.IsPermutation();
		assert(flag == true);
		flag = mn.IsPermutation();
		assert(flag == false);	
	}

	// --- quaternion tests ---
	{
		using namespace Numeric;

		Quaternion<float> q{ 0.0f, 1.0f, 0.0f, 1.0f },
						  r{ 0.5f, 0.5f, 0.75f, 1.0f },
						  mult{ q * r },
						  quat_1,
						  quat_2(1.0f),
						  quat_3(quat_2);

		assert(static_cast<int>(mult.X() * 100) == 125);
		assert(static_cast<int>(mult.Y() * 100) == 150);
		assert(static_cast<int>(mult.Z() * 100) == 25);
		assert(static_cast<int>(mult.W() * 100) == 50);

		q *= r;
		assert(static_cast<int>(q.X() * 100) == 125);
		assert(static_cast<int>(q.Y() * 100) == 150);
		assert(static_cast<int>(q.Z() * 100) == 25);
		assert(static_cast<int>(q.W() * 100) == 50);

		quat_3.Normalize();
		std::cout << "quat_3 = " << quat_3 << "\n";

		VectorN<float, 3> v{ 1.0f, 0.5f, 0.25f };
		RotatePointUsingQuaternion(v, Quaternion<float>{ -0.1314f, 0.6602f, 0.5066f, 0.5387f });
		assert(static_cast<int>(v.X() * 1000) == -600);
		assert(static_cast<int>(v.Y() * 1000) == 801);
		assert(static_cast<int>(v.Z() * 1000) == -557);

		quat_1 = quat_2.Conjugate();
		assert(static_cast<int>(quat_1[0]) == -1);
		assert(static_cast<int>(quat_1[1]) == -1);
		assert(static_cast<int>(quat_1[2]) == -1);
		assert(static_cast<int>(quat_1[3]) == 1);

		Quaternion<float> fromAxisAngle(VectorN<float, 3>{ 0.87287f, 0.43643f, 0.21821f }, 0.17453f);  // 0.17453 ~= 10[deg]
		assert(static_cast<int>(fromAxisAngle.X() * 1000) == 76);
		assert(static_cast<int>(fromAxisAngle.Y() * 1000) == 38);
		assert(static_cast<int>(fromAxisAngle.Z() * 1000) == 19);
		assert(static_cast<int>(fromAxisAngle.W() * 1000) == 996);

		VectorN<float, 3> axis{ fromAxisAngle.Axis() };
		float angle{ fromAxisAngle.Angle() };
		assert(static_cast<int>(axis.X() * 1000) == 872);
		assert(static_cast<int>(axis.Y() * 1000) == 436);
		assert(static_cast<int>(axis.Z() * 1000) == 218);
		assert(static_cast<int>(angle * 1000) == 174);

		Quaternion<float> quatFromDCM(MatrixNM<float, 3>{ -0.38513f, 0.37233f, -0.844442f,
														  -0.71928f, 0.45216f,  0.52743f,
														   0.57819f, 0.81051f,  0.09367f });
		assert(static_cast<int>(std::floor(quatFromDCM.X() * 10000)) == -1314);
		assert(static_cast<int>(std::floor(quatFromDCM.Y() * 10000)) == 6602);
		assert(static_cast<int>(std::floor(quatFromDCM.Z() * 10000)) == 5066);
		assert(static_cast<int>(std::floor(quatFromDCM.W() * 10000)) == 5386);

		MatrixNM<float, 3> dcm(quatFromDCM.DCM());
		assert(static_cast<int>(std::floor(dcm(0, 0) * 10000)) == -3852);
		assert(static_cast<int>(std::floor(dcm(0, 1) * 10000)) == 3723);
		assert(static_cast<int>(std::floor(dcm(0, 2) * 10000)) == -8445);
		assert(static_cast<int>(std::floor(dcm(1, 0) * 10000)) == -7193);
		assert(static_cast<int>(std::floor(dcm(1, 1) * 10000)) == 4521);
		assert(static_cast<int>(std::floor(dcm(1, 2) * 10000)) == 5274);
		assert(static_cast<int>(std::floor(dcm(2, 0) * 10000)) == 5782);
		assert(static_cast<int>(std::floor(dcm(2, 1) * 10000)) == 8105);
		assert(static_cast<int>(std::floor(dcm(2, 2) * 10000)) == 936);

		VectorN<float, 3> vp0{ 0.043477f, 0.036412f, 0.998391f },
						  vp1{ 0.60958f, 0.7354f, 0.29597f };
		vp0.Normalize();
		vp1.Normalize();
		Quaternion<float> FromTwoVectors(vp0, vp1);
		MatrixNM<float, 3, 3> dcm2{ FromTwoVectors.DCM() };
		assert(static_cast<int>(std::floor(dcm2(0, 0) * 10000)) == 7368);
		assert(static_cast<int>(std::floor(dcm2(0, 1) * 10000)) == -3098);
		assert(static_cast<int>(std::floor(dcm2(0, 2) * 10000)) == -6010);
		assert(static_cast<int>(std::floor(dcm2(1, 0) * 10000)) == -3294);
		assert(static_cast<int>(std::floor(dcm2(1, 1) * 10000)) == 6118);
		assert(static_cast<int>(std::floor(dcm2(1, 2) * 10000)) == -7192);
		assert(static_cast<int>(std::floor(dcm2(2, 0) * 10000)) == 5904);
		assert(static_cast<int>(std::floor(dcm2(2, 1) * 10000)) == 7277);
		assert(static_cast<int>(std::floor(dcm2(2, 2) * 10000)) == 3488);
	}

	// --- interval ---
	{
		 using namespace Numeric;

                 Interval<int32_t> i0,
                                   i1(0, 1),
                                   i2{ 1, 2 },
                                   i3(i2),
                                   i4(1, 4);
                 i1 = { 1, 0 };
         
                 assert(i0.Inf() == 0);
                 assert(i0.Sup() == 0);
         
                 assert(i1.Inf() == 1);
                 assert(i1.Sup() == 0);
         
                 assert(i2.Inf() == 1);
                 assert(i2.Sup() == 2);
         
                 assert(i3.Inf() == i2.Inf());
                 assert(i3.Sup() == i2.Sup());
         
                 assert(i4.Inf() == 1);
                 assert(i4.Sup() == 4);
         
                 assert(i2.Inf() == 1);
                 assert(i2.Sup() == 2);
         
                 i2 = { 0, 3 };
                 i1 = { 0, 1 };
                 i3 = Union(i1, i2);
                 assert(i3.Inf() == 0);
                 assert(i3.Sup() == 3);
         
                 i3 = Union(i1, i2, i4);
                 assert(i3.Inf() == 0);
                 assert(i3.Sup() == 4);
         
                 i3 = Interesection(i1, i2);
                 assert(i3.Inf() == 0);
                 assert(i3.Sup() == 1);
         
                 i3 = Interesection(i1, i2, i4);
                 assert(i3.Inf() == 1);
                 assert(i3.Sup() == 1);
         
                 i3 = IntervalMin(i1, i2);
                 assert(i3.Inf() == 0);
                 assert(i3.Sup() == 1);
         
                 i3 = IntervalMin(i1, i2, i4);
                 assert(i3.Inf() == 0);
                 assert(i3.Sup() == 1);
         
                 i3 = IntervalMax(i1, i2);
                 assert(i3.Inf() == 0);
                 assert(i3.Sup() == 3);
         
                 i3 = IntervalMax(i1, i2, i4);
                 assert(i3.Inf() == 1);
                 assert(i3.Sup() == 4);
         
         
                 i1 = { 1, 2 };
                 i2 = { 3, 4 };
                 i4 = Interesection(i1, i2);
                 assert(i4.Inf() == std::numeric_limits<int32_t>::max());
                 assert(i4.Sup() == std::numeric_limits<int32_t>::min());
         
                 i1 = { 1, 20 };
                 i2 = { 0, 1 };
                 i4 = i1 / i2;
                 assert(i4.Inf() == 1);
                 assert(i4.Sup() == std::numeric_limits<int32_t>::max());
         
                 i1 = { -1, 0 };
                 i2 = { 1, std::numeric_limits<int32_t>::max() };
                 i4 = i1 * i2;
                 assert(i4.Inf() == -std::numeric_limits<int32_t>::max());
                 assert(i4.Sup() == 0);
	}

	
	// --- 3D tests ---
	{
		using namespace Numeric;

		// SVD symmetric (2x2)
		MatrixNM<float, 2> svd2{ 1.0f, 2.0f,
								 2.0f, 4.0f },
						   usvd;
		VectorN<float, 2> wsvd;
		SVDsymmetric2x2(svd2, usvd, wsvd);

		assert(static_cast<int>(std::floor(wsvd[0])) == 5);
		assert(static_cast<int>(std::floor(wsvd[1])) == 0);

		assert(static_cast<int>(std::floor(usvd(0, 0) * 10000)) == 4472);
		assert(static_cast<int>(std::floor(usvd(0, 1) * 10000)) == -8945);
		assert(static_cast<int>(std::floor(usvd(1, 0) * 10000)) == 8944);
		assert(static_cast<int>(std::floor(usvd(1, 1) * 10000)) == 4472);

		// polar decomposition
		MatrixNM<float, 2> r22, s22,
						   f22{ 1.3f, -0.375f,
								0.75f, 0.65f };
		PolarDecomposition(f22, r22, s22);

		assert(static_cast<int>(std::floor(r22(0, 0) * 1000)) == 866);
		assert(static_cast<int>(std::floor(r22(0, 1) * 10)) == 4);
		assert(static_cast<int>(std::floor(r22(1, 0) * 10)) == -5);
		assert(static_cast<int>(std::floor(r22(1, 1) * 1000)) == 866);

		assert(static_cast<int>(std::floor(s22(0, 0) * 10)) == 15);
		assert(static_cast<int>(std::floor(s22(0, 1))) == 0);
		assert(static_cast<int>(std::abs(s22(1, 0))) == 0);
		assert(static_cast<int>(std::floor(s22(1, 1) * 100)) == 75);

		MatrixNM<float, 2> A{ 1.0f,  1.0f,
							  1.0f, -1.0f };
		bool flag = IsOrthonormal(A);
		assert(flag == false);
		A /= std::sqrt(2.0f);;
		flag = IsOrthonormal(A);
		assert(flag == true);


		MatrixNM<float, 3> B{ 2.0f, -2.0f,  1.0f,
							  1.0f,  2.0f,  2.0f,
							  2.0f,  1.0f, -2.0f };
		flag = IsOrthonormal(B);
		assert(flag == false);
		B /= 3.0f;
		flag = IsOrthonormal(B);
		assert(flag == true);

		float temp = AngleBetweenThreePoints(VectorN<float, 3>{ 1.0f, 2.0f, 3.0f },
											 VectorN<float, 3>{ 3.0f, 2.0f, 1.0f },
											 VectorN<float, 3>{ 1.0f, 1.0f, 1.0f });
		assert(static_cast<int>(std::floor(temp * 1000)) == 50768);

		// inverse 3x3
		VectorN<float, 3> vp0{ 0.043477f, 0.036412f, 0.998391f },
						  vp1{ 0.60958f, 0.7354f, 0.29597f };
		vp0.Normalize();
		vp1.Normalize();
		Quaternion<float> FromTwoVectors(vp0, vp1);
		MatrixNM<float, 3, 3> dcm{ FromTwoVectors.DCM() },
							  dcmInv{ Inverse(dcm) },
							  dcmTrans{ dcm.Transpose() };
		for (std::size_t i = 0; i < 9; ++i) {
			assert(static_cast<int>(std::floor(dcmInv[i] * 1000)) == static_cast<int>(std::floor(dcmTrans[i] * 1000)));
		}

		// inverse 4x4
		MatrixNM<float, 4> mat4{ 12.0f, 13.0f,  14.0f, 15.0f,
								 16.0f, 15.0f,  14.0f, 13.0f,
								 38.0f, 75.0f, -15.0f,  5.0f,
								 92.0f, 32.0f,  27.0f,  5.0f },
						   mat4inv;

		temp = mat4.Determinant();
		assert(static_cast<int>(temp) == 108948);

		mat4inv = Inverse(mat4);
		mat4inv *= temp;
		assert(static_cast<int>(std::floor(mat4inv(0, 0))) == 26404);
		assert(static_cast<int>(std::floor(mat4inv(1, 0))) == -31835);
		assert(static_cast<int>(std::floor(mat4inv(2, 0))) == 476);
		assert(static_cast<int>(std::floor(mat4inv(3, 0))) == 3080);

		assert(static_cast<int>(std::floor(mat4inv(0, 1))) == -31698);
		assert(static_cast<int>(std::floor(mat4inv(1, 1))) == 37350);
		assert(static_cast<int>(std::floor(mat4inv(2, 1))) == 588);
		assert(static_cast<int>(std::floor(mat4inv(3, 1))) == -2604);

		assert(static_cast<int>(std::floor(mat4inv(0, 2))) == -66402);
		assert(static_cast<int>(std::floor(mat4inv(1, 2))) == 79170);
		assert(static_cast<int>(std::floor(mat4inv(2, 2))) == -2604);
		assert(static_cast<int>(std::floor(mat4inv(3, 2))) == -4032);

		assert(static_cast<int>(std::floor(mat4inv(0, 3))) == 75585);
		assert(static_cast<int>(std::floor(mat4inv(1, 3))) == -80794);
		assert(static_cast<int>(std::floor(mat4inv(2, 3))) == 1540);
		assert(static_cast<int>(std::floor(mat4inv(3, 3))) == 3556);
	}

	return 1;
}
