#include<iostream>
#include<cmath>
#include<windows.h>
#include<algorithm>
#define SIZE_VAR 2
using namespace std;

// function квадратическая  = 30*pow(x,2)+56*x+1+4*pow(y,2)+5*y+178*x*y;
// funcio розенброка = 110*pow((y-pow(x,2)),2)+pow((3-x),2)
double t = 0, t1 = 0, t2 = 0, t3 = 0;
// t = x[0][0]   t1 = p[0][0]  t2 = x[1][0]  t3 = p[1][0]
double funcP(double x)
{
	// квадратическая функция
	/*return 6 * pow((t + t1 * x), 2) + 4 * sqrt(5) * (t + t1 * x) + 22 + 
		3 * pow((t2 + t3 * x), 2)
		+ 8 * sqrt(5) * (t2 + t3 * x) - 4 * (t + t1 * x)*(t2 + t3 * x);*/

	// функция розенброка
	return 110 * pow(((t2 + t3 * x) - pow((t + t1 * x), 2)), 2) + 
		pow((3 - (t + t1 * x)), 2);
}
long double funcBFGS(double x1, double x2)
{
	// квадратическая функция
	//return 6 * pow(x1, 2) + 4 * sqrt(5) * x1 + 
		22 + 3 * pow(x2, 2) + 8 * sqrt(5) * x2 - 4 * x1*x2;

	// функция розенброка
	return 110 * pow((x2 - pow(x1, 2)), 2) + pow((3 - x1), 2);
}
double** gradFuncBFGS(double x, double y)
{
	double** arr = new double*[SIZE_VAR];
	for (int i = 0; i < SIZE_VAR; i++)
		arr[i] = new double[1];
	/*
	аналитический ввод градиента
	*/
	// квадратическая функция
	//arr[0][0] = 12 * x + 4 * sqrt(5) - 4 * y;
	//arr[1][0] = 6 * y + 8 * sqrt(5) - 4 * x;

	// функция розенброка
	arr[0][0] = -1 * 440 * x*y + 440 * pow(x, 3) - 6 + 2 * x;
	arr[1][0] = 220 * y - 220 * pow(x, 2);
	return arr;
}

struct Matrix
{
	int height;
	int width;
	double **arr;
	Matrix() { ; }
	Matrix(int _height, int _width)
	{
		height = _height;
		width = _width;
		arr = new double*[height];
		for (int i = 0; i < height; i++)
		{
			arr[i] = new double[width];
			for (int j = 0; j < width; j++)
			{
				arr[i][j] = 0;
			}
		}
	}
	Matrix(int _height, int _width,double **_arr)
	{
		height = _height;
		width = _width;
		arr = new double*[height];
		for (int i = 0; i < height; i++)
		{
			arr[i] = new double[width];
			for (int j = 0; j < width; j++)
			{
				arr[i][j] = _arr[i][j];
			}
		}
	}
	void clear()
	{
		for (int i = 0; i < height; i++)
		{
			arr[i] = new double[width];
			for (int j = 0; j < width; j++)
			{
				arr[i][j] = 0;
			}
		}
	}
};
class MethodPowell
{
private:
	double x1, x2, x3; // точки
	double h; // величина шага
	double e1,e2; // точность
	double fMin, xMin; // мин.значения функции и точки ей соответствующей
	double xMinPol; // минимум интерп.полинома
	double fMinPol; // значение функции минимум интерп.полинома
	double answer; // минимум функции (ответ)
	bool tp = true; // вместо goto
public:
	MethodPowell(double _x1, double _h, double _e1, double _e2)
		:x1(_x1), h(_h), e1(_e1), e2(_e2) { ; }
	double* sort(double &_x1, double &_x2, double &_x3, double _xMinPol)
	{
		double temp = 0;
		double* arr = new double[4];
		arr[0] = _x1;
		arr[1] = _x2;
		arr[2] = _x3;
		arr[3] = _xMinPol;
		for (int i = 0; i < 4 - 1; i++) {
			for (int j = 0; j < 4 - i - 1; j++) {
				if (arr[j] > arr[j + 1]) {
					temp = arr[j];
					arr[j] = arr[j + 1];
					arr[j + 1] = temp;
				}
			}
		}
		if (_x3 < _x2) swap(_x3, _x2);
		if (_x2 < _x1) swap(_x2, _x1);
		return arr;
	}
	double algPow()
	{
		/*cout << "----------------------------------------\n";
		cout << "! Начался поиск шага\n";
		cout << "Квадратичная интерполяция\n";
		cout << "----------------------------------------\n";
		cout << "нач x = " << x1 << endl;
		cout << "h = " << h << endl;
		cout << "e1 = " << e1 << endl;
		cout << "e2 = " << e2 << endl;
		cout << "--------------------\n";*/
		int index = 0;
		while (true)
		{
			if (tp == true)
			{
				x2 = x1 + h; // шаг 1

				if (funcP(x1) > funcP(x2)) // шаг 2
					x3 = x1 + 2 * h;
				else
					x3 = x1 - h;
			}
			tp = true;
			fMin = min(funcP(x1), min(funcP(x2), funcP(x3))); // шаг 3
			if (fMin == funcP(x1))
				xMin = x1;
			else if (fMin == funcP(x2))
				xMin = x2;
			else
				xMin = x3;
			// шаг 4
			double numerator = (pow(x2, 2) - pow(x3, 2))*funcP(x1) + (pow(x3, 2) -
				pow(x1, 2))*funcP(x2) + (pow(x1, 2) - pow(x2, 2))*funcP(x3);
			double denominator = (x2 - x3)*funcP(x1) + (x3 - x1)*funcP(x2) + (x1 - x2)*funcP(x3);
			if (denominator == 0)
			{
				x1 = xMin;
				/*cout << "\niteration [" << index << "]" << endl;
				cout << "x1= " << x1 << endl;
				cout << "x2= " << x2 << endl;
				cout << "x3= " << x3 << endl;
				cout << "f1= " << funcP(x1) << endl;
				cout << "f2= " << funcP(x2) << endl;
				cout << "f3= " << funcP(x3) << endl;
				cout << "Знаменатель = 0. К шагу 1\n";*/
				index++;
				continue; // переход к 1 шагу
			}
			xMinPol = 0.5*(numerator / denominator);
			fMinPol = funcP(xMinPol);

			// Шаг 5
			double tmp1 = fabs((fMin - fMinPol) / fMinPol);
			double tmp2 = fabs((xMin - xMinPol) / xMinPol);
			/*cout << "\niteration [" << index << "]" << endl;
			cout << "x1= " << x1 << endl;
			cout << "x2= " << x2 << endl;
			cout << "x3= " << x3 << endl;
			cout << "f1= " << funcP(x1) << endl;
			cout << "f2= " << funcP(x2) << endl;
			cout << "f3= " << funcP(x3) << endl;
			cout << "xMin= " << xMin << endl;
			cout << "xMinPol= " << xMinPol << endl;
			cout << "funcP(xMinPol)= " << funcP(xMinPol) << endl;
			cout << "|(fMin - fMinPol) / fMinPol| = " << tmp1 << endl;
			cout << "|(xMin - xMinPol) / xMinPol| = " << tmp2 << endl;*/
			index++;
			if (tmp1 < e1 && tmp2 < e2)
			{
				answer = xMinPol;
				/*cout << "----------------------------------------\n";
				cout << "Условие остановки выполнено\n";
				cout << "Ответ: " << answer << endl;
				cout << "----------------------------------------\n";*/
				return answer;
			}
			else
			{
				double* arr = sort(x1, x2, x3, xMinPol);
				if (xMinPol >= x1 && xMinPol <= x3)
				{
					tp = false;
					double goodX = min(funcP(xMin), funcP(xMinPol));
					if (goodX == funcP(xMin))
						goodX = xMin;
					else
						goodX = xMinPol;
					bool tmp = true;
					if (goodX == arr[0])
					{
						x1 = arr[0] - h;
						x2 = arr[0];
						x3 = arr[1];
						tmp = false;
					}
					else if (goodX == arr[3])
					{
						x1 = arr[2];
						x2 = arr[3];
						x3 = arr[3] + h;
						tmp = false;;
					}
					if (tmp == true)
					{
						for (int i = 1; i < 3; i++)
						{
							if (goodX >= arr[i])
							{
								x1 = arr[i - 1];
								x2 = arr[i];
								x3 = arr[i + 1];
							}
						}
					}
				}
				else if (xMinPol < x1 || xMinPol > x3)
				{
					x1 = xMinPol;
				}
			}
		}
	}
};

class MethodBFGS
{
private:
	int index = 0; // иттерация
	int n = SIZE_VAR; // колич.переменных
	Matrix p; // направление поиска
	Matrix x0; // Xk
	Matrix x1; // Xk+1
	Matrix s0; // шаг алгоритма на итеррации Sk
	Matrix y0; // изменение градиента на итерации Yk
 	Matrix H0; // гессиан функция Hk
	double e; // точность
	double normaGrad = 0; // норма градиента для проверки условия окончания поиска
	double a; // значения соответст мин.фун при условии Вольфе

public:
	MethodBFGS(Matrix _x0, double _e) :x0(_x0), e(_e) { ; }
	void algBFGS()
	{
		Matrix gradFunc0(n, 1); // градиент функции k
		Matrix gradFunc1(n, 1); // градиент функции k+1
		Matrix iMatrix(n, n); // единичная матрица
		for (int i = 0; i < iMatrix.height; i++) // при первой иттерации единичная гессиан функция
		{
			for (int j = 0; j < iMatrix.width; j++)
			{
				if (i == j)
					iMatrix.arr[i][j] = 1;
				else
					iMatrix.arr[i][j] = 0;
			}
		}
		H0 = iMatrix; // в 0 иттерацию H=единичная матрица
		cout << "---------------Метод BFGS---------------\n----------------------------------------" << endl;
		cout << "x0(первоначальная точка)" << endl;
		outputMatrix(x0);
		cout << "e(эпсило) = " << e << endl;
		cout << "----------------------------------------\n";
		
		while (true)
		{
			cout << "iteration = " << index << endl;
			cout << "X[" << index << "]: \n";
			outputMatrix(x0);
			cout << "Значение функции[" << index << "]: " << funcBFGS(x0.arr[0][0], x0.arr[1][0]) << endl;
			cout << "H(Матрица гессе)[" << index << "]: \n";
			outputMatrix(H0);
			// шаг 1
			// нахождение градиента
			gradFunc0.clear();
			gradFunc0.arr = gradFuncBFGS(x0.arr[0][0], x0.arr[1][0]);
			cout << "Градиент[" << index << "]: \n";
			outputMatrix(gradFunc0);
			// находение нормы градиента
			normaGrad = 0;
			for (int i = 0; i < gradFunc0.height; i++)
				normaGrad += pow(gradFunc0.arr[i][0], 2);
			normaGrad = sqrt(normaGrad);
			cout << "Норма градиента[" << index << "]: " << normaGrad << endl;
			// проверка на условие окончания поиска
			if (normaGrad < e)
			{
				cout << "----------------------------------------\n";
				cout << "Условие остановки выполнено\n";
				cout << "----------------------------------------\n";
				cout << "(" << x0.arr[0][0] << " ; " << x0.arr[1][0] << ")" << endl;
				return;
			}
			// шаг2
			H0 = multipNumber(H0, -1);
			p = multipMatrix(H0, gradFunc0);
			H0 = multipNumber(H0, -1);
			cout << "P(направление поиска)[" << index << "]: \n";
			outputMatrix(p);
			// шаг 3
			t1 = p.arr[0][0];
			t3 = p.arr[1][0];
			t = x0.arr[0][0];
			t2 = x0.arr[1][0];
			MethodPowell mtPow(1, 1, 0.1, 0.1);
			a = mtPow.algPow();
			cout << "a(шаг поиска)[" << index << "]: " << a << endl;
			Matrix center; // середина выражения
			Matrix left; // Левая часть
			Matrix right; // правая часть
			Matrix y_0; // Транспонированная y0
			Matrix s_0; // Транспонированная s0
			Matrix tmp; // для подсчетов

			tmp = multipNumber(p, a);
			x1 = summaMatrix(x0, tmp);
			// шаг 4
			s0 = difMatrix(x1, x0);// разность точек
			/*cout << "S0(разность точек)[" << index << "]\n";
			outputMatrix(s0);*/
			gradFunc1.arr = gradFuncBFGS(x1.arr[0][0], x1.arr[1][0]); 
			y0 = difMatrix(gradFunc1, gradFunc0); // разность градиентов
			//cout << "Y0(разность градиентов)[" << index << "]\n";
			//outputMatrix(y0);
			// шаг 5
			//y_0 = transpose(y0);
			//s_0 = transpose(s0);

			//tmp = multipMatrix(y_0, s0);

			//left = multipMatrix(s0, y_0);
			//left = divisionNumber(left, tmp.arr[0][0]);
			//left = difMatrix(iMatrix, left);
			//left = multipMatrix(left, H0);
			//H0 = left;

			//center = multipMatrix(y0, s_0);
			//center = divisionNumber(center, tmp.arr[0][0]);
			//center = difMatrix(iMatrix, center);
			//H0 = multipMatrix(H0, center);

			//right = multipMatrix(s0, s_0);
			//right = divisionNumber(right, tmp.arr[0][0]);

			//H0 = summaMatrix(H0, right);
			//x0 = x1;

			y_0 = transpose(y0);
			s_0 = transpose(s0);

			tmp = multipMatrix(y_0, s0);
			right = multipMatrix(s0, s_0);
			right = divisionNumber(right, tmp.arr[0][0]);

			left = multipMatrix(s0, y_0);
			left = divisionNumber(left, tmp.arr[0][0]);
			left = difMatrix(iMatrix, left);

			center = multipMatrix(y0, s_0);
			center = divisionNumber(center, tmp.arr[0][0]);
			center = difMatrix(iMatrix, center);

			H0 = multipMatrix(left, H0);
			H0 = multipMatrix(H0, center);
			H0 = summaMatrix(H0, right);
			x0 = x1; 
			cout << "----------------------------------------\n";
			cout << index << " иттерация завершена\n";
			cout << "----------------------------------------\n";
			index++;
		}
	}
	Matrix outputMatrix(Matrix A)
	{
		for (int i = 0; i < A.height; i++)
		{
			for (int j = 0; j < A.width; j++)
			{
				cout << "[" << i << "]" << "[" << j << "]" << A.arr[i][j] << "\t";
			}
			cout << endl;
		}
		return A;
	}
	Matrix divisionNumber(Matrix A, double number)
	{
		Matrix myArray(A.height, A.width);
		for (int i = 0; i < myArray.height; i++)
		{
			for (int j = 0; j < myArray.width; j++)
			{
				myArray.arr[i][j] = A.arr[i][j] / number;
			}
		}
		return myArray;
	}
	Matrix multipMatrix(Matrix A, Matrix B)
	{
		Matrix myArray(A.height, B.width);
		for (int i = 0; i < A.height; i++)
		{
			for (int j = 0; j < B.width; j++)
			{
				for (int t = 0; t < A.width; t++)
				{
					myArray.arr[i][j] += A.arr[i][t] * B.arr[t][j];
				}
			}
		}
		return myArray;
	}
	Matrix multipNumber(Matrix A, double number)
	{
		Matrix myArray(A.height, A.width);
		for (int i = 0; i < myArray.height; i++)
		{
			for (int j = 0; j < myArray.width; j++)
			{
				myArray.arr[i][j] = A.arr[i][j] * number;
			}
		}
		return myArray;
	}
	Matrix summaMatrix(Matrix A, Matrix B)
	{
		Matrix myArray(A.height, A.width);
		for (int i = 0; i < A.height; i++)
		{
			for (int j = 0; j < B.width; j++)
			{
				myArray.arr[i][j] = A.arr[i][j] + B.arr[i][j];
			}
		}
		return myArray;
	}
	Matrix difMatrix(Matrix A, Matrix B)
	{
		Matrix myArray(A.height, A.width);
		for (int i = 0; i < A.height; i++)
		{
			for (int j = 0; j < B.width; j++)
			{
				myArray.arr[i][j] = A.arr[i][j] - B.arr[i][j];
			}
		}
		return myArray;
	}
	Matrix transpose(Matrix A)
	{
		Matrix myArray(A.width, A.height);
		for (int i = 0; i < myArray.height; i++)
		{
			for (int j = 0; j < myArray.width; j++)
			{
				myArray.arr[i][j] = A.arr[j][i];
			}
		}
		return myArray;
	}
};

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	Matrix t(2, 1);
	t.arr[0][0] = 1;
	t.arr[1][0] = 1;
	MethodBFGS gf(t, 0.1);
	gf.algBFGS();
	system("pause");
	return 0;
}
