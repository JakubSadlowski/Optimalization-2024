#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2TTest(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);
matrix ff3T(matrix, matrix = NAN, matrix = NAN);
matrix ff3Test(matrix, matrix = NAN, matrix = NAN);
matrix fT3a(matrix, matrix = NAN, matrix = NAN);
matrix fT3b(matrix, matrix = NAN, matrix = NAN);
matrix df3(double, matrix, matrix = NAN, matrix = NAN);
matrix ff3R(matrix, matrix = NAN, matrix = NAN);
matrix ff4T(matrix, matrix = NAN, matrix = NAN);
matrix gradient4T(matrix, matrix = NAN, matrix = NAN);
matrix ff4SD(matrix, matrix = NAN, matrix = NAN);
matrix gradientSD(matrix, matrix = NAN, matrix = NAN);