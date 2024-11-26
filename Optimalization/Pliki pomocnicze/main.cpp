/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include <iomanip>

#include"opt_alg.h"
#include <random>

void lab0();
void lab1();
void lab2();
void lab2Table3();
void lab2Iterations();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab3();
		//matrix x0(2, new double[]{ 2.8731934, 4.8817139 });
		/*matrix x0(2, new double[] { 2.8731868, 4.8817043 });
		ff2R(x0);*/
		//lab2Table3();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	double epsilon = 1e-6;
	double gamma = 1e-8;
	int Nmax = 100;
	double x0;
	double d = 0.5;
	double alpha = 1.5;
	double alpha2 = 3;
	double alpha3 = 5.5;
	//double a = 0.0001;
	//double b = 0.01;

	std::ofstream Sout_fibonacci("G:/Programowanie_projekty/C++/Optymalizacja/Optymalizacja/results_fibonacci3.csv");

	Sout_fibonacci << "x; y; f_calls" << std::endl;

	std::mt19937 gen(42);
	std::uniform_real_distribution<double> unif(-100.0, 100.0);

	for (int j = 0; j < 100; ++j)
	{
		x0 = unif(gen);

		solution exp = expansion(ff1T, x0, d, alpha3, Nmax);
		matrix interval = exp.x;
		double a = interval(0, 0);
		double b = interval(1, 0);

		solution fibonacci = fib(ff1T, a, b, epsilon);
		Sout_fibonacci << fibonacci.x(0) << "; " << ff1T(fibonacci.x) << solution::f_calls << std::endl;
		solution::clear_calls();
	}
	/*
	// Metoda Fibonacciego
	solution fibonacci_solution = fib(ff1R, a, b, epsilon);
	cout << "Fibonacci Solution D_A: " << fibonacci_solution.x(0) << "\ny: " << ff1T(fibonacci_solution.x) << "f_calls: " << solution::f_calls << endl;
	solution::clear_calls();

	// Metoda Lagrange'a
	solution lagrange_solution = lag(ff1R, a, b, epsilon, gamma, Nmax);
	cout << "Lagrange Solution D_A: " << lagrange_solution.x(0) << "\ny: " << ff1T(lagrange_solution.x) << "f_calls: " << solution::f_calls << endl;
	solution::clear_calls();
	*/

	Sout_fibonacci.close();
}


void lab2()
{
	double epsilon = 0.00001;
	int Nmax = 100000;
	double s = 0.1;
	double alpha = 0.5;
	double alpha2 = 2.0;
	double beta = 0.5;

	std::ofstream sout("results_combined.csv");
	sout << "D³ugoœæ kroku;Lp.;x1(0);x2(0);x1* (HJ);x2* (HJ);y* (HJ);Liczba wywo³añ funkcji celu (HJ);x1* (Rosen);x2* (Rosen);y* (Rosen);Liczba wywo³añ funkcji celu (Rosen)" << std::endl;

	std::mt19937 gen(42);
	std::uniform_real_distribution<double> unif(-1.0, 1.0);

	for (int j = 1; j <= 100; ++j)
	{
		double initial_values[2] = {unif(gen), unif(gen)};
		double initial_values2[2] = { 0.1, 0.1 };
		matrix x0(2, initial_values);
		matrix s0(2, initial_values2);

		// Hooke-Jeeves
		solution hj = HJ(ff2T, x0, s, alpha, epsilon, Nmax);
		int f_calls_hj = solution::f_calls;
		solution::clear_calls();

		// Rosenbrock
		solution rosen = Rosen(ff2T, x0, s0, alpha2, beta, epsilon, Nmax);
		int f_calls_rosen = solution::f_calls;
		solution::clear_calls();

		sout << s << "; "  << j << "; " << x0(0) << "; " << x0(1) << "; " << hj.x(0) << "; " << hj.x(1) << "; " << ff2T(hj.x) << f_calls_hj << "; " << rosen.x(0) << "; " << rosen.x(1) << "; " << ff2T(rosen.x) << f_calls_rosen << std::endl; 
	}

	sout.close();
}

void lab2Table3()
{
	double epsilon = 0.00001;
	int Nmax = 100000;
	double s = 0.1;
	double alpha = 0.5;
	double alpha2 = 2.0;
	double beta = 0.5;

	std::mt19937 gen(42);
	std::uniform_real_distribution<double> unif(-1.0, 1.0);

	
	double initial_values[2] = { 0.3, 0.3 };
	double initial_values2[2] = { 0.1, 0.1 };
	matrix x0(2, initial_values);
	matrix s0(2, initial_values2);

	// Hooke-Jeeves
	solution hj = HJ(ff2R, x0, s, alpha, epsilon, Nmax);
	int f_calls_hj = solution::f_calls;
	solution::clear_calls();

	// Rosenbrock
	solution rosen = Rosen(ff2R, x0, s0, alpha2, beta, epsilon, Nmax);
	int f_calls_rosen = solution::f_calls;
	solution::clear_calls();

	std::cout << setprecision(8) << "k1_HJ" << "; " << hj.x(0) << "; " << "k2_HJ" << "; " << hj.x(1) << "; " << "Q_HJ" << "; " << ff2R(hj.x) << "fcalls_HJ" << "; " << f_calls_hj << "; "
	<< "k1_Rosen" << "; " << rosen.x(0) << "; " << "k2_Rosen" << "; " << rosen.x(1) << "; " << "Q_Rosen" << "; " << ff2R(rosen.x) << "fcalls_Rosen" << "; " << f_calls_rosen << std::endl;
}

void lab2Iterations()
{
	double epsilon = 0.00001;
	int Nmax = 100000;
	double s = 0.1;
	double alpha = 0.5;
	double alpha2 = 2.0;
	double beta = 0.5;

	double initial_values[2] = { -0.5, 0.5 };
	double initial_values2[2] = { 0.1, 0.1 };
	matrix x0(2, initial_values);
	matrix s0(2, initial_values2);

	// Hooke-Jeeves
	/*solution hj = HJ(ff2T, x0, s, alpha, epsilon, Nmax);
	cout << solution::f_calls << "\n";
	solution::clear_calls();*/

	// Rosenbrock
	solution rosen = Rosen(ff2R, x0, s0, alpha2, beta, epsilon, Nmax);
	int f_calls_rosen = solution::f_calls;
	cout << solution::f_calls << "\n";
	solution::clear_calls();
}

void lab3()
{	
	double s = 1.0;
	double alpha = 1;
	double beta = 0.5;
	double gamma = 2.0;
	double delta = 0.5;
	double epsilon = 1e-3;
	
	double cIntern = 10.0;
	double dcIntern = 0.5;

	double cExtern = 0.5;
	double dcExtern = 2.0;
	
	matrix ud1(5);
	matrix ud2Intern(cIntern);
	matrix ud2Extern(cExtern);
	int Nmax = 10000;
	solution symplexNelder;
	solution penIn, penOut;
	std::ofstream soutPenIn("PenIn.csv");
	std::ofstream soutPenOut("PenOut.csv");

	std::mt19937 gen(42);
	std::uniform_real_distribution<double> unif(-1.0, 1.0);

	for (int j = 1; j <= 100; ++j)
	{
		double initial_values[2] = { unif(gen), unif(gen) };
		matrix X0(2, initial_values);		
		
		penIn = pen(fT3a, X0, cIntern, dcIntern, epsilon, Nmax, ud1, ud2Intern);
		int f_calls_penin = solution::f_calls;
		solution::clear_calls();
		soutPenIn << std::setprecision(10) << penIn.x(0) << ";" << penIn.x(1) << ";" << penIn.y(0) << ";" << f_calls_penin << '\n';

		penOut = pen(fT3b, X0, cExtern, dcExtern, epsilon, Nmax, ud1, ud2Extern);
		int f_calls_penOut = solution::f_calls;
		solution::clear_calls();
		soutPenOut << std::setprecision(10) << penOut.x(0) << ";" << penOut.x(1) << ";" << penOut.y(0) << ";" << f_calls_penOut << '\n';
	}

	soutPenIn.close();
	soutPenOut.close();
}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
