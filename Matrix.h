#pragma once
#include <vector>

template<typename T>
class Matrix
{

public:
	size_t ncols, nrows;
	std::vector<T> internal_vector;

	Matrix() : ncols(0), nrows(0) {};
	Matrix(const size_t rows, const size_t cols) : ncols(cols), nrows(rows), internal_vector(std::vector<T>(cols*rows)) {};
	Matrix(const size_t rows, const size_t cols, const std::vector<T>& vec) : ncols(cols), nrows(rows), internal_vector(vec) {};
	~Matrix() {};

	T* operator[] (const size_t row) { return &internal_vector[row*ncols]; }
	const T* operator[] (const size_t row) const { return &internal_vector[row*ncols]; }

	void resize(const size_t rows, const size_t cols)
	{
		ncols = cols;
		nrows = rows;
		internal_vector.resize(ncols*nrows);
	};

	void clear()
	{
		ncols = 0;
		nrows = 0;
		internal_vector.clear();
	}

	bool square() const { return ncols == nrows ? true : false; }

	void print() {
		std::cout << "Printing matrix with (rows x cols) " << nrows << " x " << ncols << std::endl;
		for (size_t m = 0; m < nrows; m++)
		{
			for (size_t n = 0; n < ncols; n++)
			{
				std::cout << internal_vector[m*ncols + n] << "  ";
			}
			std::cout << std::endl;
		}
	};

	void swap(Matrix<T>& rhs)
	{
		std::swap(internal_vector, rhs.internal_vector);
		std::swap(nrows, rhs.nrows);
		std::swap(ncols, rhs.ncols);
	}

	// Assign
	Matrix<T>& operator=(Matrix<T>& rhs)
	{
		Matrix<T> temp(rhs);
		swap(temp);
		return *this;
	}

	// Scalar operators
	Matrix<T>& operator+=(T& a);
	Matrix<T>& operator-=(T& a);
	Matrix<T>& operator*=(T a);

	Matrix<T> operator-() const;
	bool operator==(const Matrix& rhs) const;


	// Matrix multiplication
	Matrix<T> operator*(Matrix<T>& b);

	void transpose(); //also implemented is transpose(mat) which makes Matrix
	void inverse();   //also implement is function inverse(mat) which returns newmat
};

template<typename T>
inline bool Matrix<T>::operator==(const Matrix & rhs) const
{
	if (nrows != rhs.nrows || ncols != rhs.ncols) return false;
	for (size_t i = 0; i < nrows; ++i)
	{
		const double* ai = operator[](i);
		const double* bi = rhs[i];
		for (size_t j = 0; j < ncols; ++j)
		{
			if (fabs(ai[j] - bi[j]) > 1.0e-12) return false;
		}
	}
	return true;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator-() const
{
	Matrix<T> tmp(nrows, ncols);
	for (size_t i = 0; i < ncols*nrows; i++)
	{
		tmp.internal_vector[i] = -internal_vector[i];
	}
	return tmp;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator+=(T& a)
{
	for (size_t i = 0; i < ncols*nrows; i++)
	{
		internal_vector[i] = internal_vector[i] + a;
	}
	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator-=(T& a)
{
	for (size_t i = 0; i < ncols*nrows; i++)
	{
		internal_vector[i] = internal_vector[i] - a;
	}
	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator*=(T a)
{
	for (size_t i = 0; i < ncols*nrows; i++)
	{
		internal_vector[i] = internal_vector[i] * a;
	}
	return *this;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator*(Matrix<T>& b)
{
	if (ncols != b.nrows) {
		std::cout << "DIMENSIONALTY PROBLEMS MATRIX MULT stop alt" << std::endl;
	}

	const size_t rows = nrows, cols = b.ncols, n = ncols;

	Matrix<T> c(rows, cols);
	//  outermost loop on result rows
	for (size_t i = 0; i < rows; ++i)
	{
		//  loop on result columns
		for (size_t j = 0; j < cols; ++j)
		{
			//  compute dot product
			T res = 0.0;
			for (size_t k = 0; k < n; ++k)
			{
				res = res + internal_vector[i*ncols + k] * b[k][j];
			}
			c[i][j] = res;
		}   //  columns
	}   //  rows
	return c;
}

template<typename T>
inline void Matrix<T>::transpose()
{
	std::vector<T> tmp(ncols*nrows);
	for (size_t i = 0; i < ncols; i++)
	{
		for (size_t j = 0; j < nrows; j++)
		{
			tmp[nrows*i + j] = internal_vector[ncols*j + i];
		}
	}
	std::swap(ncols, nrows);
	internal_vector = tmp;
}

template<typename T>
Matrix<T> transpose(const Matrix<T>& mat) {
	Matrix<T> tmp = mat;
	tmp.transpose();
	return tmp;
}

template<typename T>
Matrix<T> identity(size_t n) {
	Matrix<T> tmp(n, n);

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (i == j) tmp[i][i] = 1.0;
			else tmp[i][j] = 0.0;
		}
	}

	return tmp;
}



template<typename T>
inline void Matrix<T>::inverse()
{
	Matrix<T> a(nrows, ncols);
	a.internal_vector = internal_vector;


	size_t n = ncols;
	size_t m = n;
	if (n != m)
		std::cout << "ERROR in .inverse(): not square matrix" << std::endl;

	size_t icol, irow;
	T big, dum, pivinv;

	Matrix<T> b = identity<T>(n);

	std::vector<size_t> indxc(n);
	std::vector<size_t> indxr(n);
	std::vector<size_t> ipiv(n);

	// Iterators
	size_t i, j, k, l, ll;

	for (j = 0; j < n; j++)
	{
		ipiv[j] = 0;
	}

	for (i = 0; i < n; i++)
	{
		big = 0.0;
		for (j = 0; j < n; j++)
		{
			if (ipiv[j] != 1)
				for (k = 0; k < n; k++)
				{
					if (ipiv[k] == 0) {
						if (abs(a[j][k]) >= big) {
							big = abs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}
		}
		ipiv[icol] = ipiv[icol] + 1;
		if (irow != icol)
		{
			for (l = 0; l < n; l++) { std::swap(a[irow][l], a[icol][l]); }
			for (l = 0; l < m; l++) { std::swap(b[irow][l], b[icol][l]); }
		}
		indxr[i] = irow;
		indxc[i] = icol;


		if (abs(a[icol][icol]) < 10e-7) {
			std::cout << "SINGULAR" << std::endl;
		}
		pivinv = 1.0 / a[icol][icol];
		for (l = 0; l < n; l++) { a[icol][l] = a[icol][l] * pivinv; }
		for (l = 0; l < m; l++) { b[icol][l] = b[icol][l] * pivinv; }
		for (ll = 0; ll < n; ll++)
		{
			if (ll != icol)
			{
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 0; l < n; l++) { a[ll][l] = a[ll][l] - a[icol][l] * dum; }
				for (l = 0; l < m; l++) { b[ll][l] = b[ll][l] - b[icol][l] * dum; }
			}
		}
	}
	swap(b);
}

template<typename T>
Matrix<T> inverse(const Matrix<T>& mat) {
	Matrix<T> tmp = mat;
	tmp.inverse();
	return tmp;
}

template<typename T>
Matrix<T> diag(const std::vector<T>& vec, const size_t rows, const size_t cols) {
	Matrix<T> tmp(rows, cols);
	for (size_t i = 0; i < rows*cols; i++) tmp.internal_vector[i] = 0;

	for (size_t i = 0; i < vec.size(); i++)
	{
		tmp[i][i] = vec[i];
	}
	return tmp;
}

template<typename T>
void diaginverse(Matrix<T>& m) {
	for (size_t i = 0; i < m.ncols; i++)
	{
		if (abs(m[i][i]) > 10e-4) {
			m[i][i] = 1 / m[i][i];
		}
		else
		{
			m[i][i] = 0.0;
		}
	}
}

template<typename T>
Matrix<double> as_double(Matrix<T>& mat) {
	Matrix<double> out(mat.nrows, mat.ncols);

	out.internal_vector.resize(mat.internal_vector.size());
	for (size_t i = 0; i < mat.internal_vector.size(); i++)
	{
		out[i] = value(mat.internal_vector[i]);
	}
	return out;
}

template<typename T>
std::vector<double> as_double(std::vector<T>& vec) {
	std::vector<double> out(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		out[i] = value(vec[i]);
	}
	return out;
}

template<typename T>
T pythag(T a, T b) {
	T absa, absb;
	absa = abs(a);
	absb = abs(b);
	if (absa > absb) return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb)));
}

template<typename T>
T SIGN(T& a, T& b) {
	return b >= 0.0 ? abs(a) : -abs(b);
}

template<typename T>
int svd(Matrix<T>& a, std::vector<T>& w, Matrix<T>& v)
{
	// Takes matrix a, and decomposes it to vector w, and matrix v while turning a into u
	// a = u * diag(w) * v^T

	//http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c
	int flag, i, its, j, jj, k, l, nm;
	T c, f, h, s, x, y, z;
	T anorm = 0.0, g = 0.0, scale = 0.0;

	int m = (int)a.nrows;
	int n = (int)a.ncols;

	v.resize(n, n);
	w.resize(n);

	if (m < n)
	{
		std::cout << "#rows must be > #cols " << std::endl;
		return(0);
	}

	std::vector<T> rv1(n);


	/* Householder reduction to bidiagonal form */
	for (i = 0; i < n; i++)
	{
		/* left-hand reduction */
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i < m)
		{
			for (k = i; k < m; k++)
				scale = scale + abs((T)a[k][i]);
			if (value(scale))
			{
				for (k = i; k < m; k++)
				{
					a[k][i] = (T)((T)a[k][i] / scale);
					s = s + ((T)a[k][i] * (T)a[k][i]);
				}
				f = (T)a[i][i];


				T tmp_s = sqrt(s) > 0.0 ? sqrt(s) : -sqrt(s);
				g = value(f) >= 0.0 ? tmp_s : -tmp_s;
				g = -g;

				h = f * g - s;
				a[i][i] = (T)(f - g);
				if (i != n - 1)
				{
					for (j = l; j < n; j++)
					{
						for (s = 0.0, k = i; k < m; k++)
							s = s + ((T)a[k][i] * (T)a[k][j]);
						f = s / h;
						for (k = i; k < m; k++)
							a[k][j] = a[k][j] + (T)(f * (T)a[k][i]);
					}
				}
				for (k = i; k < m; k++)
					a[k][i] = (T)((T)a[k][i] * scale);
			}
		}
		w[i] = (T)(scale * g);

		/* right-hand reduction */
		g = s = scale = 0.0;
		if (i < m && i != n - 1)
		{
			for (k = l; k < n; k++)
				scale = scale + abs((T)a[i][k]);
			if (value(scale))
			{
				for (k = l; k < n; k++)
				{
					a[i][k] = (T)((T)a[i][k] / scale);
					s = s + ((T)a[i][k] * (T)a[i][k]);
				}
				f = (T)a[i][l];

				T tmp_s = sqrt(s) > 0.0 ? sqrt(s) : -sqrt(s);
				g = value(f) >= 0.0 ? tmp_s : -tmp_s;
				g = -g;


				h = f * g - s;
				a[i][l] = (T)(f - g);
				for (k = l; k < n; k++)
					rv1[k] = (T)a[i][k] / h;
				if (i != m - 1)
				{
					for (j = l; j < m; j++)
					{
						for (s = 0.0, k = l; k < n; k++)
							s = s + ((T)a[j][k] * (T)a[i][k]);
						for (k = l; k < n; k++)
							a[j][k] = a[j][k] + (T)(s * rv1[k]);
					}
				}
				for (k = l; k < n; k++)
					a[i][k] = (T)((T)a[i][k] * scale);
			}
		}
		anorm = value(anorm) > (abs((T)w[i]) + abs(rv1[i])) ? anorm : (abs((T)w[i]) + abs(rv1[i]));
	}

	/* accumulate the right-hand transformation */
	for (i = n - 1; i >= 0; i--)
	{
		if (i < n - 1)
		{
			if (value(g))
			{
				for (j = l; j < n; j++)
					v[j][i] = (T)(((T)a[i][j] / (T)a[i][l]) / g);
				/* T division to avoid underflow */
				for (j = l; j < n; j++)
				{
					for (s = 0.0, k = l; k < n; k++)
						s = s + ((T)a[i][k] * (T)v[k][j]);
					for (k = l; k < n; k++)
						v[k][j] = v[k][j] + (T)(s * (T)v[k][i]);
				}
			}
			for (j = l; j < n; j++)
				v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}

	/* accumulate the left-hand transformation */
	for (i = n - 1; i >= 0; i--)
	{
		l = i + 1;
		g = (T)w[i];
		if (i < n - 1)
			for (j = l; j < n; j++)
				a[i][j] = 0.0;
		if (value(g))
		{
			g = 1.0 / g;
			if (i != n - 1)
			{
				for (j = l; j < n; j++)
				{
					for (s = 0.0, k = l; k < m; k++)
						s = s + ((T)a[k][i] * (T)a[k][j]);
					f = (s / (T)a[i][i]) * g;
					for (k = i; k < m; k++)
						a[k][j] = a[k][j] + (T)(f * (T)a[k][i]);
				}
			}
			for (j = i; j < m; j++)
				a[j][i] = (T)((T)a[j][i] * g);
		}
		else
		{
			for (j = i; j < m; j++)
				a[j][i] = 0.0;
		}
		a[i][i] = a[i][i] + 1;
	}

	/* diagonalize the bidiagonal form */
	for (k = n - 1; k >= 0; k--)
	{                             /* loop over singular values */
		for (its = 0; its < 30; its++)
		{                         /* loop over allowed iterations */
			flag = 1;
			for (l = k; l >= 0; l--)
			{                     /* test for splitting */
				nm = l - 1;
				if (abs(rv1[l]) + anorm == anorm)
				{
					flag = 0;
					break;
				}
				if (abs((T)w[nm]) + anorm == anorm)
					break;
			}
			if (flag)
			{
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++)
				{
					f = s * rv1[i];
					if (abs(f) + anorm != anorm)
					{
						g = (T)w[i];
						h = pythag(f, g);
						w[i] = (T)h;
						h = 1.0 / h;
						c = g * h;
						s = (-f * h);
						for (j = 0; j < m; j++)
						{
							y = (T)a[j][nm];
							z = (T)a[j][i];
							a[j][nm] = (T)(y * c + z * s);
							a[j][i] = (T)(z * c - y * s);
						}
					}
				}
			}
			z = (T)w[k];
			if (l == k)
			{                  /* convergence */
				if (z < 0.0)
				{              /* make singular value nonnegative */
					w[k] = (T)(-z);
					for (j = 0; j < n; j++)
						v[j][k] = (-v[j][k]);
				}
				break;
			}
			if (its >= 30) {
				fprintf(stderr, "No convergence after 30,000! iterations \n");
				return(0);
			}

			/* shift from bottom 2 x 2 minor */
			x = (T)w[l];
			nm = k - 1;
			y = (T)w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);

			T one = 1.0;

			g = pythag(f, one);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

			/* next QR transformation */
			c = s = 1.0;
			for (j = l; j <= nm; j++)
			{
				i = j + 1;
				g = rv1[i];
				y = (T)w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for (jj = 0; jj < n; jj++)
				{
					x = (T)v[jj][j];
					z = (T)v[jj][i];
					v[jj][j] = (T)(x * c + z * s);
					v[jj][i] = (T)(z * c - x * s);
				}
				z = pythag(f, h);
				w[j] = (T)z;
				if (value(z))
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = (c * g) + (s * y);
				x = (c * y) - (s * g);
				for (jj = 0; jj < m; jj++)
				{
					y = (T)a[jj][j];
					z = (T)a[jj][i];
					a[jj][j] = (T)(y * c + z * s);
					a[jj][i] = (T)(z * c - y * s);
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = (T)x;
		}
	}
	return 1;
}

template<typename T>
void matrix_to_txt(std::string filename, Matrix<T>& matrix) {

	// This function takes in a filename (location) and a matrix and returns a .txt file. 
	std::ofstream myfile;
	myfile.open(filename, std::fstream::out);

	for (unsigned int i = 0; i < matrix.nrows; i++)
	{
		for (unsigned int j = 0; j < matrix.ncols; j++)
		{
			myfile << matrix[i][j] << "\t";
		}
		myfile << std::endl;
	}
	myfile.close();
}