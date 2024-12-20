﻿#include"opt_alg.h"
#include <cmath>
#include <iomanip>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

solution expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution interval;
		interval.x = matrix(2, 1);
		int i = 0;
		solution X0 = x0;
		solution X1 = X0.x + d;

		if (X1.fit_fun(ff) == X0.fit_fun(ff))
		{
			interval.x(0, 0) = X0.x(0);
			interval.x(1, 0) = X1.x(0);
			return interval;
		}

		if (X1.fit_fun(ff) > X0.fit_fun(ff))
		{
			d = -d;
			X1 = X0.x + d;
			if (X1.fit_fun(ff) >= X0.fit_fun(ff))
			{
				interval.x(0, 0) = X1.x(0);
				interval.x(1, 0) = X0.x(0) - d;
				return interval;
			}
		}

		solution xPre;

		while (true)
		{
			if (i >= Nmax) throw "Maximum number of function calls exceeded";

			i++;
			xPre = X0;
			X0 = X1;
			X1.x(0) = X0.x(0) + pow(alpha, i) * d;

			if (X0.fit_fun(ff) <= X1.fit_fun(ff))
				break;
		}

		if (d > 0)
		{
			interval.x(0, 0) = xPre.x(0);
			interval.x(1, 0) = X1.x(0);
		}
		else
		{
			interval.x(0, 0) = X1.x(0);
			interval.x(1, 0) = xPre.x(0);
		}

		return interval;
	}
	catch (string ex_info)
	{
		throw ("solution expansion(...):\n" + ex_info);
	}
}

double generateFibonacci(int n) {
	return (pow(1. + sqrt(5.), n) - pow(1. - sqrt(5.), n)) / (pow(2., n) * sqrt(5.));;
}

solution fib(matrix(ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int k = 0;

		while (generateFibonacci(k) <= (b - a) / epsilon)
		{
			k++;
		}

		solution A0 = a;
		solution B0 = b;
		solution C0, D0;
		C0.x(0) = B0.x(0) - (generateFibonacci(k - 1) / generateFibonacci(k)) * (B0.x(0) - A0.x(0));
		D0.x(0) = A0.x(0) + B0.x(0) - C0.x(0);

		for (int i = 0; i < k - 3; i++)
		{
			if (C0.fit_fun(ff) < D0.fit_fun(ff))
			{
				B0.x(0) = D0.x(0);
			}
			else
			{
				A0.x(0) = C0.x(0);
			}

			C0.x(0) = B0.x(0) - (generateFibonacci(k - i - 2) / generateFibonacci(k - i - 1)) * (B0.x(0) - A0.x(0));
			D0.x(0) = A0.x(0) + B0.x(0) - C0.x(0);
		}

		Xopt.x = C0.x(0);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}


solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution A0 = a;
		solution B0 = b;
		solution C0 = 0.5 * (A0.x + B0.x);
		solution D0;
		solution D1, A1, B1, C1;
		solution l, m;
		int i = 0;

		while (true)
		{
			if (i > Nmax)
				throw "Maximum number of function calls exceeded";

			l = A0.fit_fun(ff) * (pow(B0.x, 2) - pow(C0.x, 2)) +
				B0.fit_fun(ff) * (pow(C0.x, 2) - pow(A0.x, 2)) +
				C0.fit_fun(ff) * (pow(A0.x, 2) - pow(B0.x, 2));
			m = A0.fit_fun(ff) * (B0.x - C0.x) +
				B0.fit_fun(ff) * (C0.x - A0.x) +
				C0.fit_fun(ff) * (A0.x - B0.x);

			if (m.x <= 0)
			{
				Xopt.x(0) = 666;
				Xopt.y = 666;
				return Xopt;
			}

			D0 = 0.5 * l.x / m.x;

			if (A0.x < D0.x && D0.x < C0.x)
			{
				if ((D0.fit_fun(ff) - 50) < (C0.fit_fun(ff) - 50))
				{
					A1 = A0;
					B1 = C0;
					C1 = D0;
				}
				else
				{
					A1 = D0;
					C1 = C0;
					B1 = B0;
				}
			}
			else if (C0.x < D0.x && D0.x < B0.x)
			{
				if (D0.fit_fun(ff) < C0.fit_fun(ff))
				{
					A1 = C0;
					C1 = D0;
					B1 = B0;
				}
				else
				{
					A1 = A0;
					C1 = C0;
					B1 = D0;
				}
			}
			else
			{
				Xopt.x(0) = 666;
				Xopt.y(0) = 666;
				return Xopt;
			}

			if ((B0.x - A0.x < epsilon) || (fabs(D0.x(0) - D1.x(0)) < gamma))
				break;
			D1.x(0) = D0.x(0);
			A0 = A1;
			B0 = B1;
			C0 = C1;
			D1 = D0;
			i++;
		}

		Xopt = D0;
		Xopt.flag = 0;
		return Xopt;
	}

	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}


solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		/*std::ofstream iterFile("HJ_iterations.csv");
		iterFile << "Iteration; x1; x2\n";*/

		solution Xopt;
		solution x = x0;
		solution xb;
		//int iteration = 0;
	
		do
		{
			xb = x;
			x = HJ_trial(ff, xb, s);
			if (x.fit_fun(ff) < xb.fit_fun(ff))
			{
				while (x.fit_fun(ff) < xb.fit_fun(ff))
				{
					solution tempxb = xb;
					xb = x;
					x = 2 * xb.x - tempxb.x;
					x = HJ_trial(ff, x, s);

					/*iterFile << std::fixed << setprecision(8) << iteration << "; " << xb.x(0) << "; " << xb.x(1) << "\n";
					iteration++;*/

					if (solution::f_calls > Nmax)
						throw "Maximum number of function calls exceeded";
				}
				x = xb;
			}
			else
			{
				s = alpha * s;
			}

			if (solution::f_calls > Nmax)
				throw "Maximum number of function calls exceeded";
		} while (s >= epsilon);

		//iterFile.close();

		Xopt = xb;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_len(XB.x);
		matrix ej = ident_mat(n);

		for (int j = 0; j < n; j++)
		{
			solution xbPlus = XB.x + s * ej[j];
			solution xbMinus = XB.x - s * ej[j];
			if (xbPlus.fit_fun(ff) < XB.fit_fun(ff))
			{
				XB = xbPlus;
			}
			else if (xbMinus.fit_fun(ff) < XB.fit_fun(ff))
			{
				XB = xbMinus;
			}
		}
		
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

double max_abs(const matrix& vec) {
	int n = get_len(vec);
	double max_val = fabs(vec(0, 0));
	for (int i = 1; i < n; i++) {
		max_val = std::max(max_val, fabs(vec(i, 0)));
	}
	return max_val;
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		/*std::ofstream iterFile("Rosen_iterations.csv");
		iterFile << "Iteration; x1; x2\n";*/
		solution Xopt;
		solution xB = x0;          
		matrix s = s0;              
		int n = get_len(x0);
		matrix dj = ident_mat(n);   
		matrix lambda(n, 1);        
		matrix p(n, 1);
		//int iteration = 0;

		do {
			for (int j = 0; j < n; j++) 
			{
				matrix step = s(j) * get_col(dj, j); 
				solution xBNew = xB.x + step;

				if (xBNew.fit_fun(ff) < xB.fit_fun(ff)) 
				{
					xB = xBNew;
					lambda(j) = lambda(j) + s(j);
					s(j) = alpha * s(j);
				}
				else 
				{
					s(j) = -beta * s(j);
					p(j) = p(j) + 1;
				}
			}

			/*iterFile << std::fixed << setprecision(8) << iteration << "; " << xB.x(0) << "; " << xB.x(1) << "\n";
			iteration++;*/

			if (solution::f_calls > Nmax)
				throw "Maximum number of function calls exceeded";

			bool changeDirection = true;
			for (int j = 0; j < n; j++) 
			{
				if (lambda(j) == 0 || p(j) == 0) 
				{
					changeDirection = false;
					break;
				}
			}

			if (changeDirection) 
			{
				matrix Q(n, n);
				matrix vj(n, 1);
				for (int i = 0; i < n; ++i) 
				{
					for (int j = 0; j <= i; ++j) 
					{
						Q(i, j) = lambda(i); 
					}
				}

				Q = dj * Q;
				vj = Q[0] / norm(Q[0]);
				dj.set_col(vj, 0);

				for (int j = 0; j < n; j++) 
				{
					matrix sum(n, 1);
					for (int k = 0; k < j; k++) 
					{
						sum = sum + (trans(Q[j]) * dj[k]) * dj[k];
					}
					vj = Q[j] - sum;
					dj.set_col(vj, j);
				}
				lambda = matrix(n, 1);
				p = matrix(n, 1);
				s = s0;
			}

		} while (max_abs(s) >= epsilon); 

		//iterFile.close();
		Xopt = xB;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try 
	{
		int iterator = 0;
		solution Xopt(x0);
		solution xPrev(x0);
		Xopt.fit_fun(ff, ud1, ud2);

		do {
			xPrev = Xopt;			
			Xopt = sym_NM(ff, Xopt.x, 0.1, 1.0, 0.5, 2.0, 0.5, epsilon, Nmax, ud1, ud2);		
			ud2(0) *= dc;
			iterator++;

			if (solution::f_calls > Nmax)
				throw "Maximum number of function calls exceeded";

		} while (norm(Xopt.x - xPrev.x) >= epsilon);

		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

void sort(solution& min, solution& mid, solution& max) {
	if (min.y(0) > mid.y(0)) {
		std::swap(min, mid);
	}
	if (min.y(0) > max.y(0)) {
		std::swap(min, max);
	}
	if (mid.y(0) > max.y(0)) {
		std::swap(mid, max);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		// Kod dla n = 2
		int iter = 1;
		solution Xopt;
		matrix e1(2, new double[2] { 1.0, 0.0 });
		matrix e2(2, new double[2] { 0.0, 1.0 });
		solution p0 = x0;
		solution p1 = p0.x + s * e1;
		solution p2 = p0.x + s * e2;
		solution minPoint, midPoint, maxPoint;
		solution centerOfGravity;
		solution pReflection;
		solution pExpansion;
		solution pReduction;
		matrix emptyMatrix(2, new double[2] { 0.0, 0.0 });
		
		minPoint = p0;
		midPoint = p1;
		maxPoint = p2;
		minPoint.fit_fun(ff, ud1, ud2);
		midPoint.fit_fun(ff, ud1, ud2);
		maxPoint.fit_fun(ff, ud1, ud2);

		sort(minPoint, midPoint, maxPoint);

		do
		{		
			centerOfGravity.x = emptyMatrix;
			centerOfGravity.x(0) = (minPoint.x(0) + midPoint.x(0)) / 2;
			centerOfGravity.x(1) = (minPoint.x(1) + midPoint.x(1)) / 2;

			//odbicie
			pReflection.x = emptyMatrix;
			pReflection.x(0) = centerOfGravity.x(0) + alpha * (centerOfGravity.x(0) - maxPoint.x(0));
			pReflection.x(1) = centerOfGravity.x(1) + alpha * (centerOfGravity.x(1) - maxPoint.x(1));

			//ekspansja			
			if (pReflection.fit_fun(ff, ud1, ud2) < minPoint.fit_fun(ff, ud1, ud2))
			{
				pExpansion.x = emptyMatrix;
				pExpansion.x(0) = centerOfGravity.x(0) + gamma * (pReflection.x(0) - centerOfGravity.x(0));
				pExpansion.x(1) = centerOfGravity.x(1) + gamma * (pReflection.x(1) - centerOfGravity.x(1));
				if (pExpansion.fit_fun(ff, ud1, ud2) < pReflection.fit_fun(ff, ud1, ud2))
				{
					maxPoint = pExpansion;
				}	
				else 
				{
					maxPoint = pReflection;
				}					
			}
			else 
			{
				if (minPoint.fit_fun(ff, ud1, ud2) <= pReflection.fit_fun(ff, ud1, ud2) && pReflection.fit_fun(ff, ud1, ud2) < maxPoint.fit_fun(ff, ud1, ud2))
					maxPoint = pReflection;
				else 
				{
					//redukcja
					pReduction.x = emptyMatrix;
					pReduction.x(0) = centerOfGravity.x(0) + beta * (maxPoint.x(0) - centerOfGravity.x(0));
					pReduction.x(1) = centerOfGravity.x(1) + beta * (maxPoint.x(1) - centerOfGravity.x(1));

					if (pReduction.fit_fun(ff, ud1, ud2) >= maxPoint.fit_fun(ff, ud1, ud2))
					{
						maxPoint.x(0) = delta * (maxPoint.x(0) + minPoint.x(0));
						maxPoint.x(1) = delta * (maxPoint.x(1) + minPoint.x(1));

						midPoint.x(0) = delta * (midPoint.x(0) + minPoint.x(0));
						midPoint.x(1) = delta * (midPoint.x(1) + minPoint.x(1));
					}
					else maxPoint = pReduction;
				}
			}

			minPoint.fit_fun(ff, ud1, ud2);
			midPoint.fit_fun(ff, ud1, ud2);
			maxPoint.fit_fun(ff, ud1, ud2);

			sort(minPoint, midPoint, maxPoint);

			if (solution::f_calls > Nmax) 
				throw "Maximum number of function calls exceeded";

		} while (sqrt(pow((maxPoint.x(0) - minPoint.x(0)), 2) + pow(maxPoint.x(1) - minPoint.x(1), 2)) > epsilon || sqrt(pow((midPoint.x(0) - minPoint.x(0)), 2) + pow(midPoint.x(1) - minPoint.x(1), 2)) > epsilon);

		Xopt = minPoint;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

double calculateStepSize(bool constH, matrix(*ff)(matrix, matrix, matrix), int Nmax,
	double hPrev, const matrix& dPrev, const matrix& xPrev) {
	if (constH) return hPrev;

	matrix ud(2, 2, 0.0);

	if (get_size(xPrev)[0] > 0 && get_size(xPrev)[1] > 0 &&
		get_size(dPrev)[0] > 0 && get_size(dPrev)[1] > 0) {
		ud(0, 0) = xPrev(0, 0);
		ud(0, 1) = get_size(xPrev)[0] > 1 ? xPrev(1, 0) : 0.0;
		ud(1, 0) = dPrev(0, 0);
		ud(1, 1) = get_size(dPrev)[0] > 1 ? dPrev(1, 0) : 0.0;
	}

	solution exp = expansion(ff, hPrev, 0.5, 1.2, Nmax, matrix(1, 1, 0.0), ud);
	return golden(ff, exp.x(0), exp.x(1), 0.001, Nmax, matrix(1, 1, 0.0), ud).x(0);
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	try {
		solution Xopt;
		solution x, xNext;
		matrix d;
		double h = h0;
		x.x = x0;
		xNext.x = x0;
		int i = 0;
		//d = gf(x0, matrix(1, 1, 0.0), matrix(1, 1, 0.0)); // Pass null matrices as matrix(1,1,0.0)
		//d = -d; // Use matrix negation operator
		bool constH;
		if (h == 0)
			constH = false;
		else 
			constH = true;

		do {
			matrix d = xNext.grad(gf);
			d = -d;
			h = calculateStepSize(constH, ff, Nmax, h, d, xNext.x);
			x = xNext;
			xNext.x = x.x + d * h;

			if (solution::f_calls > Nmax) {
				throw string("Max fcalls");
			}

			i++;
			if (i == 20) {
				Xopt = xNext;
				Xopt.fit_fun(ff);
				return Xopt;
			}
		} while (norm(xNext.x - x.x) >= epsilon);

		Xopt = xNext;
		Xopt.fit_fun(ff);
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	try {
		solution Xopt;
		solution x, xNext;
		matrix d, g0, g1;
		double beta;
		double h = h0;
		x.x = x0;
		xNext.x = x0;
		int i = 0;
		bool constH;
		if (h == 0)
			constH = false;
		else
			constH = true;

		g0 = gf(x0, ud1, ud2);
		d = -g0;

		do {
			h = calculateStepSize(constH, ff, Nmax, h, d, xNext.x);
			x = xNext;
			xNext.x = x.x + d * h;

			if (solution::f_calls > Nmax) {
				throw string("Max fcalls");
			}

			matrix g1 = xNext.grad(gf);
			beta = pow(norm(g1), 2) / pow(norm(g0), 2);  
			d = -g1 + beta * d;
			g0 = g1;

			i++;
			if (i == 20) {
				Xopt = xNext;
				Xopt.fit_fun(ff);
				return Xopt;
			}
		} while (norm(xNext.x - x.x) >= epsilon);

		Xopt = xNext;
		Xopt.fit_fun(ff);
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
	int Nmax, matrix ud1, matrix ud2) {
	try {
		solution Xopt;
		solution x, xNext;
		matrix d;
		double h = h0;
		x.x = x0;
		xNext.x = x0;
		int i = 0;
		bool constH;
		if (h == 0)
			constH = false;
		else
			constH = true;

		do {
			matrix H = xNext.hess(Hf);
			matrix g = xNext.grad(gf);
			d = -inv(H) * g;

			h = calculateStepSize(constH, ff, Nmax, h, d, xNext.x);;
			x = xNext;
			xNext.x = x.x + d * h;

			if (solution::f_calls > Nmax) {
				throw string("Max fcalls");
			}
			i++;
			if (i == 20) {
				Xopt = xNext;
				Xopt.fit_fun(ff);
				return Xopt;
			}
		} while (norm(xNext.x - x.x) >= epsilon);

		Xopt = xNext;
		Xopt.fit_fun(ff);
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double alpha = (sqrt(5.0) - 1.0) / 2.0;
		solution Xa(a);
		solution Xb(b);
		solution Xc(b - alpha * (b - a));
		solution Xd(a + alpha * (b - a));

		matrix fc = Xc.fit_fun(ff, ud1, ud2);
		matrix fd = Xd.fit_fun(ff, ud1, ud2);

		while (true)
		{
			if (solution::f_calls > Nmax)
			{
				Xopt = Xc; 
				Xopt.flag = 0; 
				return Xopt;
			}

			if (abs(Xb.x() - Xa.x()) < epsilon)
			{
				Xopt.x = (Xa.x + Xb.x) / 2.0;
				Xopt.y = Xopt.fit_fun(ff, ud1, ud2); 
				Xopt.flag = 1; 
				return Xopt;
			}

			if (m2d(fc) < m2d(fd))
			{
				Xb = Xd;
				Xd = Xc;
				fd = fc;
				Xc.x = Xb.x - alpha * (Xb.x - Xa.x);
				fc = Xc.fit_fun(ff, ud1, ud2);
			}
			else
			{
				Xa = Xc;
				Xc = Xd;
				fc = fd;
				Xd.x = Xa.x + alpha * (Xb.x - Xa.x);
				fd = Xd.fit_fun(ff, ud1, ud2);
			}
		}
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
