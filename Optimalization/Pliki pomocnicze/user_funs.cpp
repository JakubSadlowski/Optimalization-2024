#include"user_funs.h"
#define _USE_MATH_DEFINES
#include <iomanip>
#include<math.h>

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * m2d(x)) * pow(exp(1), -pow(0.1 * m2d(x) - 2 * M_PI, 2)) + 0.002 * pow(0.1 * m2d(x), 2);
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) 
{
	// dane
	double a = 0.98; // wsp odpowiadajacy za lepkosc cieczy
	double b = 0.63; // wsp odpowiadajacy za zwezenie strumienia cieczy
	double g = 9.81; // przysp ziemskie
	double Pa = 0.5;
	double Pb = 1.0;
	double FbIn = 0.01;
	double Da = m2d(ud2);
	double Db = 0.00365665;
	double Va = Y(0);
	double Vb = Y(1);
	double Ta0 = 90.0;
	double Tb0 = Y(2);
	double TbIn = 20.0;

	matrix dY(3, 1);

	double FaOut, FbOut;

	if (Y(0) > 0)
	{
		FaOut = a * b * Da * sqrt(2.0 * g * Va / Pa);
	}
	else
	{
		FaOut = 0;
	}


	if (Y(1) > 0)
	{
		FbOut = a * b * Db * sqrt(2.0 * g * Vb / Pb);
	}
	else
	{
		FbOut = 0;
	}

	dY(0) = -1.0 * FaOut;
	dY(1) = FaOut + FbIn - FbOut;
	dY(2) = FbIn / Vb * (TbIn - Tb0) + FaOut / Vb * (Ta0 - Tb0);
	return dY;
}

//matrix ff1R(matrix X, matrix ud1, matrix ud2) 
// {
//	matrix Yi;
//	matrix Y0 = matrix(3, new double[3] {5.0, 1.0, 20.0});
//	// 5 - Va, 1 - Vb, 20 - Tb0
//
//	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, X);
//
//	int n = get_len(Y[0]);
//	double max = Y[1](0, 2);
//	for (int i = 0; i < n; i++)
//	{
//		if (max < Y[1](i, 2)) {
//			max = Y[1](i, 2);
//		}
//	}
//	return abs(max - 50);
//	//return max;
//}

matrix ff1R(matrix X, matrix ud1, matrix ud2) 
{
	matrix Yi;
	matrix Y0 = matrix(3, new double[3] {5.0, 1.0, 20.0});
	// 5 - Va, 1 - Vb, 20 - Tb0

	matrix* lag = solve_ode(df1, 0, 1, 2000, Y0, ud1, 0.00116775);
	matrix* fib = solve_ode(df1, 0, 1, 2000, Y0, ud1, 0.00116724);

	std::ofstream Sout_fibonacci("results_fibonacci4.csv");
	Sout_fibonacci << "t; FVA; LVA; FVB; LVB; FTB; LTB" << std::endl;
	int n = get_len(lag[0]);
	double max = lag[1](0, 2);
	for (int i = 0; i < n; i++)
	{
		Sout_fibonacci << i << "; " << fib[1](i, 0) << "; " << lag[1](i, 0) << "; " << fib[1](i, 1) << "; " << lag[1](i, 1) << "; " << fib[1](i, 2) << "; " << lag[1](i, 2) << "\n";
	}
	Sout_fibonacci.close();
	return abs(max - 50);
	//return max;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) 
{
	matrix y;
	y = pow(m2d(x(0)), 2) + pow(m2d(x(1)), 2) - cos(2.5 * M_PI * m2d(x(0))) - cos(2.5 * M_PI * m2d(x(1))) + 2;
	return y;
}

matrix ff2TTest(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = 2.5 * pow(pow(x(0), 2) - x(1), 2) + pow(1 - x(0), 2);
	return y;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
	//dane
	double l = 1.;
	double mr = 1.;
	double mc = 5.;
	double b = 0.5;
	double I = 1. / 3. * mr * pow(l, 2) + mc * pow(l, 2);
	matrix dY(2, 1);

	dY(0) = Y(1);
	dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1)) / I;

	return dY;
}

//matrix ff2R(matrix X, matrix ud1, matrix ud2) {
//	matrix y = 0;
//	matrix Y0(2, 1);
//	matrix Yref(2, new double[2] { 3.14, 0. });
//	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Yref, X);
//	int n = get_len(Y[0]);
//	for (int i = 0; i < n; ++i) {
//		y = y + 10 * pow(Yref(0) - Y[1](i, 0), 2) + pow(Yref(1) - Y[1](i, 1), 2) + pow(X(0) * (Yref(0) - Y[1](i, 0)) + X(1) * (Yref(1) - Y[1](i, 1)), 2);
//		std::cout << i << std::setprecision(12) << " " << Y[1](i, 0) << " " << Y[1](i, 1) << " " << y(0) << std::endl;
//	}
//	
//	y = y * 0.1;
//	return y;
//}

matrix ff2R(matrix X, matrix ud1, matrix ud2) {
	matrix y = 0;
	matrix Y0(2, 1);
	matrix Yref(2, new double[2] { 3.14, 0. });
	
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Yref, X);
	std::ofstream Sout_simulation("results_simulation_rosen.csv");
	Sout_simulation << "t; alpha; omega" << std::endl;

	int n = get_len(Y[0]);
	for (int i = 0; i < n; ++i) {
		y = y + 10 * pow(Yref(0) - Y[1](i, 0), 2) + pow(Yref(1) - Y[1](i, 1), 2) + pow(X(0) * (Yref(0) - Y[1](i, 0)) + X(1) * (Yref(1) - Y[1](i, 1)), 2);
		Sout_simulation << i << "; " << setprecision(7) << Y[1](i, 0) << "; " << Y[1](i, 1) << "\n";
	}

	y = y * 0.1;
	Sout_simulation.close();
	return y;
}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	//Podstawowa funkcja celu
	y = sin((M_PI * sqrt(pow(x(0)/M_PI, 2) + pow(x(1) / M_PI, 2))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))));
	return y;
}

//Zewnętrzna funkcja kary
matrix fT3a(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = sin((M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))));

	if (-x(0) + 1 > 0) {
		y = y + ud2(0) * pow(-x(0) + 1, 2);
	}

	if (-x(1) + 1 > 0) {
		y = y + ud2(0) * pow(-x(1) + 1, 2);
	}

	if (norm(x) - ud1(0) > 0) {
		y = y + ud2(0) * pow(norm(x) - ud1(0), 2);
	}

	return y;
}

//Wewnętrzna funkcja kary
matrix fT3b(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = sin((M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))));

	if (-x(0) + 1 > 0)
		y = 1 * exp(10);
	else
		y = y - ud2(0) / (-x(0) + 1);

	if(-x(1) + 1 > 0)
		y = 1 * exp(10);
	else
		y = y - ud2(0) / (-x(1) + 1);

	if (norm(x) - ud1(0) > 0)
		y = 1 * exp(10);
	else
		y = y - (ud2(0) / norm(x) - ud1(0));

	return y;
}

matrix ff3Test(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2);
	return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {

	matrix dY(4, 1);

	double C = 0.47;
	double p = 1.2;
	double r = 0.12;
	double S = 3.14 * r * r;
	double m = 0.6;
	double g = 9.81;

	double Dx = 1.0 / 2.0 * C * p * S * Y(1) * abs(Y(1));
	double Fmx = p * Y(3) * ud1(0) * S * r;

	double Dy = 1.0 / 2.0 * C * p * S * Y(3) * abs(Y(3));
	double Fmy = p * Y(1) * ud1(0) * S * r;

	if (Y(2) > 0) {
		dY(0) = Y(1);

		dY(1) = (-Fmx - Dx) / m;

		dY(2) = Y(3);

		dY(3) = (-m * g - Dy - Fmy) / m;
	}
	else {

		dY(0) = 0;

		dY(1) = 0;

		dY(2) = 0;

		dY(3) = 0;

	}
	return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {

	matrix y;
	matrix Y0(4, new double[4] {0, x(0), 100, 0});
	matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, x(1), ud2);
	int n = get_len(Y[0]);
	int i50 = 0, i0 = 0;
	for (int i = 0; i < n; i++) {
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50)) i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2))) i0 = i;
	}
	y = -Y[1](i0, 0);

	if (abs(x(0)) - 10 > 0) y = y + ud2(0) * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 15 > 0) y = y + ud2(0) * pow(abs(x(1)) - 15, 2);
	if (abs(Y[1](i50, 0) - 5) - 0.5 > 0) y = y + ud2(0) * pow(abs(Y[1](i50, 0) - 5) - 0.5, 2);
	return y;
}

matrix ff4SD(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2);
	return y;
}

matrix gradientSD(matrix x, matrix ud1, matrix ud2) {
	matrix g(2, 1);
	g(0) = 2 * x(0);
	g(1) = 2 * x(1);
	return g;
}

matrix ff4Newton(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - 40 / (pow(x(0), 2) + pow(x(1), 2) + 1);
	return y;
}

matrix gradientNewton(matrix x, matrix ud1, matrix ud2) {
	matrix g(2, 1);
	double denominator = pow(pow(x(0), 2) + pow(x(1), 2) + 1, 2);

	// Derivative for x(0)
	g(0) = 2 * x(0) + 80 * x(0) / denominator;

	// Derivative for x(1)
	g(1) = 2 * x(1) + 80 * x(1) / denominator;

	return g;
}

matrix ff4Golden(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = x + 1 / pow(x, 2);
	return y;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x(0) + 2.0 * x(1) - 7.0, 2) + pow(2.0 * x(0) + x(1) - 5.0, 2);
	return y;
}

matrix gradient4T(matrix x, matrix ud1, matrix ud2) {
	matrix g(2);
	g(0) = 2 * (x(0) + 2 * x(1) - 7) + 4 * (2 * x(0) + x(1) - 5);
	g(1) = 4 * (x(0) + 2 * x(1) - 7) + 2 * (2 * x(0) + x(1) - 5);
	return g;
}

matrix gf4T(matrix x, matrix ud1, matrix ud2) {
	double valueX1 = 10 * x(0) + 8 * x(1) - 34;
	double valueX2 = 8 * x(0) + 10 * x(1) - 38;
	matrix y(2, new double[2] {valueX1, valueX2});
	return y;
}

matrix hf4T(matrix x, matrix ud1, matrix ud2) {
	double x0 = x(0);
	double x1 = x(1);

	double h00 = 10;
	double h10 = 8;
	double h01 = 8;
	double h11 = 10;

	matrix y(2, 2);
	y(0, 0) = h00;
	y(1, 0) = h10;
	y(0, 1) = h01;
	y(1, 1) = h11;
	return y;
}