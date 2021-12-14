#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdlib>   
#include <ctime> 
#include <map>
#include <assert.h>

//using namespace std;
const double PI = 3.141592653589793238;
template<typename T>
//void pop_front(std::vector<T>& vec)
//{
//	assert(!vec.empty());
//	vec.front() = std::move(vec.back());
//	vec.pop_back();
//}

const double eps = 1e-7;
int N = 10, l = 2;
template <class T> const T& max(const T& a, const T& b) {
	return (a < b) ? b : a;
}
template <class T> const T& min(const T& a, const T& b) {
	return (a < b) ? a : b;
}
#pragma region Matrix
class Matrix
{
public:
	Matrix();
	Matrix(int rows, int cols, double init);
	Matrix(int N, double init);

	Matrix& operator=(const Matrix& m);
	Matrix(std::vector<std::vector<double>> A);
	~Matrix();
	Matrix transpose();
	double dotProduct(Matrix& a, Matrix& b);
	std::vector<std::vector<double>> m;
	int cols_;
	int rows_;
private:
};

Matrix::Matrix()
{
}

Matrix::Matrix(int rows, int cols, double init) {
	m.resize(rows);
	for (int i = 0; i < rows; i++) {
		m[i].resize(cols, init);
	}
	rows_ = m.size();
	cols_ = m[0].size();
}

Matrix::Matrix(int N, double init)
{
	m.resize(N);
	for (int i = 0; i < N; i++) {
		m[i].resize(N, init);
	}
	rows_ = m.size();
	cols_ = m[0].size();
}

Matrix& Matrix::operator=(const Matrix& m)
{
	this->m = m.m;
	this->cols_ = m.cols_;
	this->rows_ = m.rows_;
	return *this;
}

Matrix::Matrix(std::vector<std::vector<double>> A)
{
	m = A;
	rows_ = m.size();
	cols_ = m[0].size();
}
Matrix Matrix::transpose()
{
	Matrix ret(rows_, cols_, 0);
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			ret.m[j][i] = m[i][j];
		}
	}
	return ret;
}
double Matrix::dotProduct(Matrix& a, Matrix& b)
{
	double sum = 0;
	for (int i = 0; i < a.rows_; ++i) {
		sum += (a.m[i][0] * b.m[i][0]);
	}
	return sum;
}
Matrix::~Matrix()
{
}
Matrix operator+(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, const Matrix& b);
Matrix operator-(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, double b);
Matrix operator*(double b, const Matrix& a);
Matrix operator/(const Matrix& a, double b);

#pragma endregion
#pragma region Основые методы
//Генерирует случайное число в заданном диапазоне
double random_in_range(double a, double b);
//Генерирует матрицу по заданным условиям
void generateMatrix(Matrix& A, double q);
//Генерирует вектор по заданным условиям
//bool compare_float(double a, double b)
//{
//	return fabs(a - b) < eps;
//}

void generateVector(Matrix& x);


//std::vector <double> GaussianCoefsX = { 0 ,			-0.5773503,		0.5773503,	-0.7745967,
//										0, 			0.7745967,	-0.8611363,		-0.3399810,
//										0.3399810,	0.8611363,	-0.9061798,
//										-0.5384693, 0, 0.5384693,0.9061798,
//										-0.9324700, -0.6612094, -0.2386142,
//										0.2386142, 0.6612094, 0.9324700 };
//std::vector <double> GaussianCoefsC = { 2 ,1, 1,0.5555556,0.8888889,
//										0.5555556, 0.3478548, 0.6521451,
//										0.6521451,0.3478548, 0.4786287,
//										0.2369269, 0.5688888, 0.2369269,0.4786287,
//										0.1713245, 0.3607616,0.4679140,
//										0.4679140, 0.3607616,  0.1713245 };
#pragma region TASK1

double* GaussianMethod(int n, double** a, double* b, double* x)
{
	double d, s;
	for (int k = 1; k <= n; k++) // прямой ход
	{
		for (int j = k + 1; j <= n; j++)
		{
			d = a[j][k] / a[k][k]; // формула (1)
			for (int i = k; i <= n; i++)
			{
				a[j][i] = a[j][i] - d * a[k][i]; // формула (2)
			}
			b[j] = b[j] - d * b[k]; // формула (3)
		}

	}
	for (int k = n; k >= 1; k--) // обратный ход
	{
		d = 0;
		for (int j = k + 1; j <= n; j++)
		{
			s = a[k][j] * x[j]; // формула (4)
			d = d + s; // формула (4)
		}
		x[k] = (b[k] - d) / a[k][k]; // формула (4)
	}
	std::cout << "Korni sistemy: " << std::endl;
	for (int i = 1; i <= n; i++)
		std::cout << "x[" << i << "]=" << x[i] << " " << std::endl;
	return x;
}
double Interpolation(double* coeffs, double x) {
	return (coeffs[3] * x * x + coeffs[2] * x + coeffs[1]) / (coeffs[6] * x * x * x + coeffs[5] * x * x + coeffs[4] * x + 1);
}
void task_1() {
	std::ifstream in("input.txt");
	std::ofstream out("task_1.txt");
	int n = 6, i, j, k;

	double** a = new double* [n];
	for (i = 0; i <= n; i++)
		a[i] = new double[n];
	double** a1 = new double* [n];
	for (i = 0; i <= n; i++)
		a1[i] = new double[n];
	double* b = new double[n];
	double* x = new double[n];

	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			in >> a[i][j];
			a1[i][j] = a[i][j];
		}
		in >> b[i];
	}
	double* PolinomialCoeffs = GaussianMethod(n, a, b, x);
	double N = 10000;
	double h = (1 + 2) / N;
	for (double i = -2; i <= 1; i += h) {
		double r = Interpolation(PolinomialCoeffs, i);
		if (r > 19 || r < -19) {

			out << std::fixed << std::setprecision(12) << std::endl << std::endl;
		}
		out << std::fixed << std::setprecision(12) << r << std::endl;
	}
	out.close();
}
#pragma endregion
#pragma region TASK2
double Interpolation2(double* coeffs2, double x) {
	return (coeffs2[2] * x * x + coeffs2[1] * x + coeffs2[0]) / (coeffs2[5] * x * x * x + coeffs2[4] * x * x + coeffs2[3] * x + 1);
}
double* SolveUsingLU(std::vector<std::vector<double> > matrix, std::vector<double> rightPart, int n)
{
	// decomposition of matrix
	std::vector<std::vector<double>> lu(n);
	for (int i = 0; i < n; i++)
		lu[i].resize(n);
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += lu[i][k] * lu[k][j];
			lu[i][j] = matrix[i][j] - sum;
		}
		for (int j = i + 1; j < n; j++)
		{
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += lu[j][k] * lu[k][i];
			lu[j][i] = (1 / lu[i][i]) * (matrix[j][i] - sum);
		}
	}

	// lu = L+U-I
	// find solution of Ly = b
	std::vector<double> y(n);
	for (int i = 0; i < n; i++)
	{
		sum = 0;
		for (int k = 0; k < i; k++)
			sum += lu[i][k] * y[k];
		y[i] = rightPart[i] - sum;
	}
	// find solution of Ux = y
	double* x = new double[n];
	for (int i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (int k = i + 1; k < n; k++)
			sum += lu[i][k] * x[k];
		x[i] = (1 / lu[i][i]) * (y[i] - sum);
	}
	return x;
}
void task_2() {
	std::ifstream in("input.txt");
	std::ofstream out2("task_2.txt");

	int n = 6;
	std::vector<std::vector<double>> LES(n);
	std::vector<double> B(n);
	for (int i = 0; i < n; i++) {
		LES[i].resize(n);
		for (int k = 0; k < n; k++) {
			in >> LES[i][k];
		}
		in >> B[i];
	}
	double* PolinomialCoeffs2 = SolveUsingLU(LES, B, n);

	for (int i = 0; i < n; i++)
		std::cout << "x[" << i << "]=" << PolinomialCoeffs2[i] << " " << std::endl;
	double N = 10000;
	double h = (1 + 2) / N;
	for (double i = -2; i <= 1; i += h) {
		double r = Interpolation2(PolinomialCoeffs2, i);
		if (r > 3 || r < -3) {
			out2 << std::fixed << std::setprecision(12) << std::endl << std::endl;
		}
		else
			out2 << std::fixed << std::setprecision(12) << r << std::endl;
	}
	out2.close();
}


#pragma endregion
#pragma region TASK3
std::vector<std::vector<double>> mul(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B) {
	std::vector<std::vector<double>> temp(A.size());
	for (int i = 0; i < A.size(); i++) {
		temp[i].resize(B[1].size());
	}
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < temp[1].size(); ++j) {
			for (int k = 0; k < A[1].size(); ++k) {
				temp[i][j] += (A[i][k] * B[k][j]);
			}
		}
	}
	return temp;
}
std::vector<std::vector<double>> trans(std::vector<std::vector<double>> A) {
	std::vector<std::vector<double>> ret(A[1].size());
	for (int i = 0; i < A[1].size(); i++) {
		ret[i].resize(A.size());
	}
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A[1].size(); ++j) {
			ret[j][i] = A[i][j];
		}
	}
	return ret;
}
std::vector<std::vector<double>> CholeskyDecompos(std::vector<std::vector<double>> A, int n) {
	std::vector<std::vector<double>> L(n);

	for (int i = 0; i < n; i++)
	{
		L[i].resize(n); //L - треугольная матрица, поэтому в i-ой строке i+1 элементов

		double temp;
		//Сначала вычисляем значения элементов слева от диагонального элемента,
		//так как эти значения используются при вычислении диагонального элемента.
		for (int j = 0; j < i; j++)
		{
			temp = 0;
			for (int k = 0; k < j; k++)
			{
				temp += L[i][k] * L[j][k];
			}
			L[i][j] = (A[i][j] - temp) / L[j][j];
		}

		//Находим значение диагонального элемента
		temp = A[i][i];
		for (int k = 0; k < i; k++)
		{
			temp -= L[i][k] * L[i][k];
		}
		L[i][i] = sqrt(temp);
	}
	return L;
}
int SLAU(std::vector<std::vector<double>>& matrica_a, int n, std::vector<double>& massiv_b, std::vector<double>& x)

{
	int i, j, k, r;
	double c, M, max, s, ** a, * b; a = new double* [n];
	for (i = 0; i < n; i++) a[i] = new double[n];
	b = new double[n];
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) a[i][j] = matrica_a[i][j];
	for (i = 0; i < n; i++) b[i] = massiv_b[i];
	for (k = 0; k < n; k++) {
		max = fabs(a[k][k]);
		r = k;
		for (i = k + 1; i < n; i++)
			if (fabs(a[i][k]) > max)
			{
				max = fabs(a[i][k]);
				r = i;
			}
		for (j = 0; j < n; j++)
		{
			c = a[k][j];
			a[k][j] = a[r][j];
			a[r][j] = c;
		}
		c = b[k]; b[k] = b[r]; b[r] = c; for (i = k + 1; i < n; i++)
		{
			for (M = a[i][k] / a[k][k], j = k; j < n; j++)
				a[i][j] -= M * a[k][j]; b[i] -= M * b[k];
		}
	}
	if (a[n - 1][n - 1] == 0) if (b[n - 1] == 0) return -1; else return -2;
	else {
		for (i = n - 1; i >= 0; i--)
		{
			for (s = 0, j = i + 1; j < n; j++)
				s += a[i][j] * x[j]; x[i] = (b[i] - s) / a[i][i];
		}return 0;
	}
}
std::vector<std::vector<double>>  INVERSE(std::vector<std::vector<double>>& a, int n)
{
	std::vector<std::vector<double>> y(n);
	for (int i = 0; i < n; i++) {
		y[i].resize(n);
	}
	int i, j, res;
	std::vector<double> b(n);
	std::vector<double> x(n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			if (j == i) b[j] = 1; else b[j] = 0;
		res = SLAU(a, n, b, x);
		if (res != 0) break;
		else
			for (j = 0; j < n; j++)
				y[j][i] = x[j];
	}
	return y;

}
std::vector<std::vector<double>> CholeskySolve(std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& b) {
	std::vector<std::vector<double>> Li = INVERSE(L, L[1].size());
	std::vector<std::vector<double>> y = mul(Li, b);
	std::vector<std::vector<double>> LT = trans(L);
	std::vector<std::vector<double>> LTi = INVERSE(LT, L[1].size());
	std::vector<std::vector<double>> x = mul(LTi, y);
	return x;
}
double function(double x) {
	return sin(4 * pow(x, 2));
}
double aproximation(std::vector<std::vector<double>>& b, double x) {
	double Rx = 0;
	for (int i = 0; i < b.size(); i++) {
		Rx += b[i][0] * pow(x, i);
	}
	return Rx;
}
void task_3() {
	std::ofstream out3("task_3.txt");
	std::ofstream out4("task_3.1.txt");
	for (int n = 2; n < 13; n++) {
		std::vector<std::vector<double>> matrix(n + 1);
		std::vector<std::vector<double>> b(n + 1);
		for (int i = 0; i <= n; i++) {
			matrix[i].resize(n + 1);
			b[i].resize(1);
		}
		double Count = 20;
		double step = (1 - 0) / Count;
		std::vector<std::vector<double>> xvector(Count);
		std::vector<std::vector<double>> yvector(Count);
		double temp = step;
		for (int i = 0; i < Count; i++) {
			xvector[i].resize(1);
			yvector[i].resize(1);
			xvector[i][0] = temp;
			yvector[i][0] = function(xvector[i][0]);
			temp += step;


		}
		/****************************
		* Заполняем матрицу и вектор*
		*****************************/
		for (int i = 0; i <= n; i++) {

			for (int j = 0; j <= n; j++) {
				if (i == 0 && j == 0) {
					matrix[i][j] = Count;
				}
				double sum = 0;
				for (int k = 0; k < Count; k++) {
					sum += pow(xvector[k][0], i + j);
				}
				matrix[i][j] = sum;
			}

			double sum2 = 0;
			for (int p = 0; p < Count; p++) {
				sum2 += yvector[p][0] * pow(xvector[p][0], i);
			}
			b[i][0] = sum2;
		}
		std::cout << "System of equations for  n = " << n << std::endl;
		for (auto a : matrix) {
			for (auto k : a) {
				std::cout << std::fixed << std::setprecision(2) << k << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::vector<std::vector<double>> L = CholeskyDecompos(matrix, n);
		std::vector<std::vector<double>> output2 = CholeskySolve(L, b);

		std::cout << "System solution for n = " << n << std::endl;
		for (auto a : output2) {
			for (auto k : a) {
				std::cout << std::fixed << std::setprecision(2) << k << "\t";
			}
			std::cout << std::endl;
		}
		double N = 100000;
		double stepx = (1) / N;
		double x = stepx;
		double sum = 0;
		for (int i = 0; i < N; i++) {
			double a = aproximation(output2, x);
			double b = function(x);

			sum += pow(b - a, 2);

			out3 << a << std::endl;
			out4 << b << std::endl;
			x += stepx;
		}
		std::cout << std::endl << "The value of the sum of squared deviations for  n = " << n << "-->" << sum << std::endl << std::endl;
	}
	out3.close();
	out4.close();
}
#pragma endregion
#pragma region TASK4


void task_4(Matrix& A, Matrix& ff) {
#pragma region INIT
	std::vector<double> a(N + 1);
	std::vector<double> b(N + 1);
	std::vector<double> c(N + 1);
	std::vector<double> d(N + 1);
	std::vector<double> e(N + 1);
	std::vector<double> f(N + 1);
	for (int i = 0; i <= N; i++) {
		if (i >= 2) {
			a[i] = A.m[i][i - 2];
		}
		if (i >= 1) {
			b[i] = A.m[i][i - 1];
			if (b[i - 1] > 0) b[i - 1] = b[i - 1] * -1;
		}
		if (i < N - 1) {
			e[i] = A.m[i][i + 2];
		}
		if (i < N) {
			d[i] = A.m[i][i + 1];
			if (d[i] > 0) d[i] *= -1;
		}
		c[i] = A.m[i][i];
		f[i] = ff.m[i][0];
	}

#pragma endregion
	std::vector<double> alpha(N + 1), betha(N), gamma(N + 2);


	//i = 0 
	alpha[1] = d[0] / c[0];
	betha[1] = e[0] / c[0];
	gamma[1] = f[0] / c[0];
	//i = 1
	alpha[2] = (d[1] - b[1] * betha[1]) / (c[1] - b[1] * alpha[1]);
	betha[2] = e[1] / (c[1] - b[1] * alpha[1]);
	gamma[2] = (f[1] + b[1] * gamma[1]) / (c[1] - b[1] * alpha[1]);
	//i = 2 ... N - 2
	for (int i = 2; i <= N - 2; i++) {
		double delta = c[i] - a[i] * betha[i - 1] + alpha[i] * (a[i] * alpha[i - 1] - b[i]);
		alpha[i + 1] = 1 / delta * (d[i] + betha[i] * (a[i] * alpha[i - 1] - b[i]));
		betha[i + 1] = e[i] / delta;
		gamma[i + 1] = 1 / delta * (f[i] - a[i] * gamma[i - 1] - gamma[i] * (a[i] * alpha[i - 1] - b[i]));
	}
	int i = N - 1;
	double delta = c[i] - a[i] * betha[i - 1] + alpha[i] * (a[i] * alpha[i - 1] - b[i]);
	alpha[i + 1] = 1 / delta * (d[i] + betha[i] * (a[i] * alpha[i - 1] - b[i]));
	gamma[i + 1] = 1 / delta * (f[i] - a[i] * gamma[i - 1] - gamma[i] * (a[i] * alpha[i - 1] - b[i]));

	i = N;
	delta = c[i] - a[i] * betha[i - 1] + alpha[i] * (a[i] * alpha[i - 1] - b[i]);
	gamma[i + 1] = 1 / delta * (f[i] - a[i] * gamma[i - 1] - gamma[i] * (a[i] * alpha[i - 1] - b[i]));

	//Восстанавливаем ответ
	std::vector <double> y(N + 1);
	Matrix answer(N + 1, 1, 0);
	y[N] = gamma[N + 1];
	y[N - 1] = alpha[N] * y[N] + gamma[N];
	for (int i = N - 2; i >= 0; i--) {
		y[i] = alpha[i + 1] * y[i + 1] - betha[i + 1] * y[i + 2] + gamma[i + 1];
	}
	int k = 0;
	for (auto a : y) {
		answer.m[k][0] = a;
		k++;
	}
	Matrix m = A * answer;
	std::cout << "q = 10000" << std::endl;
	std::cout << "The resulting solution  \tFirst part \tResulted right part \tError margin" << std::endl;
	for (int j = 0; j <= N; j++) {
		std::cout << std::fixed << std::setprecision(8) << y[j] << "\t\t\t" << f[j] << "\t" << m.m[j][0] << "\t\t" << abs(f[j] - m.m[j][0]) << std::endl;
	}
}

#pragma endregion
int main() {

	srand(time(NULL));
	/*******************
	* Начальные условия*
	*	  Задание 4    *
	********************/
	//Исходные матрицы
	Matrix A;
	N++;
	Matrix x = Matrix(N, 1, 0.0);

	double q = 10'000;

	generateMatrix(A, q);

	generateVector(x);


	for (int i = 0; i < N; i++) {
		for (int k = 0; k < N; k++) {
			std::cout << std::fixed << std::setprecision(2) << A.m[i][k] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	for (int i = 0; i < N; i++) {
		std::cout << std::fixed << std::setprecision(5) << x.m[i][0] << std::endl;
	}

	std::cout << std::endl;
	setlocale(LC_ALL, "Russian");
	while (1) {
		std::cout << "enter number for task(1,2,3,4), 0 for exit" << std::endl;
		int taskN;
		std::cin >> taskN;
		switch (taskN)
		{
		case 1:
			task_1();
			break;
		case 2:

			task_2();
			break;
		case 3:

			task_3();
			break;
		case 4:
			N--;
			task_4(A, x);
			break;
		case 0:
			return 0;
			break;
		}
	}

}



Matrix mul(Matrix& a, Matrix& b)
{
	Matrix temp(a.rows_, b.cols_, 0);
	for (int i = 0; i < temp.rows_; ++i) {
		for (int j = 0; j < temp.cols_; ++j) {
			for (int k = 0; k < a.cols_; ++k) {
				temp.m[i][j] += (a.m[i][k] * b.m[k][j]);
			}
		}
	}
	return temp;
}
void matrixOut(Matrix& A) {
	for (auto a : A.m) {
		std::cout << std::fixed << std::setprecision(4) << a[0] << "\t";

		std::cout << std::endl;
	}
}
Matrix operator+(const Matrix& a, const Matrix& b)
{
	Matrix temp(a.rows_, b.cols_, 0);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] + b.m[i][j];
		}
	}
	return temp;
}
Matrix operator*(const Matrix& a, const Matrix& b)
{
	Matrix temp(a.rows_, b.cols_, 0);
	for (int i = 0; i < a.rows_; i++)
	{
		for (int j = 0; j < b.cols_; j++)
		{
			temp.m[i][j] = 0;
			for (int k = 0; k < b.rows_; k++)
			{
				temp.m[i][j] += a.m[i][k] * b.m[k][j];
			}
		}
	}
	return temp;
}
Matrix operator-(const Matrix& a, const Matrix& b)
{
	Matrix temp(a.rows_, b.cols_, 0);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] - b.m[i][j];
		}
	}
	return temp;
}
Matrix operator*(const Matrix& a, double b)
{
	Matrix temp(a.rows_, a.cols_, 0);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] * b;
		}
	}
	return temp;
}
Matrix operator*(double b, const Matrix& a)
{
	Matrix temp(a.rows_, a.cols_, 0);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] * b;
		}
	}
	return temp;
}
Matrix operator/(const Matrix& a, double b)
{
	Matrix temp(a.rows_, a.cols_);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] / b;
		}
	}
	return temp;
}
double random_in_range(double a, double b)
{
	double f = (double)rand() / RAND_MAX;
	return a + f * (b - a);
}
void generateMatrix(Matrix& A, double q) {
	Matrix TEMP(N, 0);
	for (signed i = 0; i < signed(N); ++i)
		for (signed j = 0; j < signed(N); j++) {
			if (j >= max(0, i - l) && j <= min(i + l, N)) {
				TEMP.m[i][j] = random_in_range(-1, 1);
			}
			else {
				TEMP.m[i][j] = 0;
			}
		}
	for (signed j = 0; j < signed(N); j++) {
		double sum = 0;

		for (signed k = 0; k < signed(N); k++) {
			if (k != j)
				sum += abs(TEMP.m[k][j]);
		}
		TEMP.m[j][j] = q * sum;
	}

	A = TEMP;
}

void generateVector(Matrix& x)
{
	for (int i = 0; i < N; i++) {
		x.m[i][0] = (random_in_range(-1, 1));
	}
}

