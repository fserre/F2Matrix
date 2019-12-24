#pragma once
#include <assert.h>
#include <fstream>
#include <immintrin.h>
#include <vector>
#include <iostream>

#define T unsigned long long
#ifdef __AVX2__
#define TZCNT(x) (_tzcnt_u64(x))
#else
#define TZCNT(x) ((x)?(_tzcnt_u64(x)):64)
#endif
using namespace std;


template <int m, int n>
struct Matrix {
	T value;
	Matrix(const T value) :value(value) {}
	Matrix() :value(0){}
	void print() const
	{
		for (int i = 0; i < m; i++)
		{
			cout << "(";
			for (int j = 0; j < n; j++)
				if ((value >> (8 * i + j)) & 1uLL)
					cout << "1";
				else
					cout << " ";
			cout << ")" << endl;
		}
		cout << endl;
	}
	void printTex(ofstream &of) const
	{
		of << "\\begin{pmatrix}" << endl;
		for (int i = 0; i < m; i++)
		{

			for (int j = 0; j < n; j++)
			{
				if ((value >> (8 * i + j)) & 1uLL)
					of << "1";
				else
					of << " ";
				if (j != m - 1)
					of << " & ";
			}
			of << "\\\\" << endl;
		}
		of<<"\\end{pmatrix}"<< endl;
	}
	template <int i, int j,int h,int w>
	Matrix <h,w> getSubMatrix() const
	{
		T mask = (
			((1uLL << w) - 1) * 0x0101010101010101uLL
			) & (
			(1uLL<<(h*8))-1
				);
		return Matrix<h, w>((value >> (j + i * 8))&mask);
	}
	Matrix<n, m> transpose() const
	{
		T res = value; T t;
		t = (res ^ (res >> 7)) & 0x00AA00AA00AA00AAuLL; res ^= t ^ (t << 7);
		t = (res ^ (res >> 14)) & 0x0000CCCC0000CCCCuLL; res ^= t ^ (t << 14);
		t = (res ^ (res >> 28)) & 0x00000000F0F0F0F0uLL; res ^= t ^ (t << 28);
		return Matrix<n, m>(res);
	}
	Matrix inverse() const
	{
		assert(m == n);
#define inverseSwap(mat,m,k) t = ((mat >> k) ^ mat) & m; mat ^= t ^ (t << k);

		T t; T pos;
		T mat = value;
		T res = identity().value;
		pos = TZCNT(mat & 0x0101010101010101uLL);
		inverseSwap(mat, 0xFF, pos);
		inverseSwap(res, 0xFF, pos);
		t = mat & 0x0101010101010100uLL;
		mat ^= t*(mat & 0xFF);
		res ^= t*(res & 0xFF);
		if (n > 1)
		{
			pos = TZCNT(mat & 0x0202020202020200uLL) - 1 - 8;
			inverseSwap(mat, 0xFF00, pos);
			inverseSwap(res, 0xFF00, pos);
			t = (mat & 0x0202020202020002uLL) >> 1;
			mat ^= t *((mat & 0xFF00) >> 8);
			res ^= t *((res & 0xFF00) >> 8);
		}
		if (n > 2)
		{
			pos = TZCNT(mat & 0x0404040404040000uLL) - 2 - 16;
			inverseSwap(mat, 0xFF0000, pos);
			inverseSwap(res, 0xFF0000, pos);
			t = (mat & 0x0404040404000404uLL) >> 2;
			mat ^= t *((mat & 0xFF0000) >> 16);
			res ^= t *((res & 0xFF0000) >> 16);
		}
		if (n > 3)
		{
			pos = TZCNT(mat & 0x0808080808000000uLL) - 3 - 24;
			inverseSwap(mat, 0xFF000000, pos);
			inverseSwap(res, 0xFF000000, pos);
			t = (mat & 0x0808080800080808uLL) >> 3;
			mat ^= t *((mat & 0xFF000000) >> 24);
			res ^= t *((res & 0xFF000000) >> 24);
		}
		if (n > 4)
		{
			pos = TZCNT(mat & 0x1010101000000000uLL) - 4 - 32;
			inverseSwap(mat, 0xFF00000000, pos);
			inverseSwap(res, 0xFF00000000, pos);
			t = (mat & 0x1010100010101010uLL) >> 4;
			mat ^= t *((mat & 0xFF00000000) >> 32);
			res ^= t *((res & 0xFF00000000) >> 32);
		}
		if (n > 5)
		{
			pos = TZCNT(mat & 0x2020200000000000uLL) - 5 - 40;
			inverseSwap(mat, 0xFF0000000000, pos);
			inverseSwap(res, 0xFF0000000000, pos);
			t = (mat & 0x2020002020202020uLL) >> 5;
			mat ^= t*((mat & 0xFF0000000000) >> 40);
			res ^= t*((res & 0xFF0000000000) >> 40);
		}
		if (n > 6)
		{
			pos = TZCNT(mat & 0x4040000000000000uLL) - 6 - 48;
			inverseSwap(mat, 0xFF000000000000, pos);
			inverseSwap(res, 0xFF000000000000, pos);
			t = (mat & 0x4000404040404040uLL) >> 6;
			mat ^= t *((mat & 0xFF000000000000) >> 48);
			res ^= t *((res & 0xFF000000000000) >> 48);
		}
		if (n == 8)
		{
			t = (mat & 0x0080808080808080uLL) >> 7;
			mat ^= t*(mat >> 56);
			res ^= t*(res >> 56);
		}
		return Matrix(mat == identity().value ? res : 0);
	}
	int rk() const
	{
		if (m < n)
			return transpose().rk();
		
#define inverseSwap(mat,m,k) t = ((mat >> k) ^ mat) & m; mat ^= t ^ (t << k);
		
		int res = 0;
		T t; T pos;
		T mat = value;


		t = mat & 0x0101010101010101uLL;
		if (t == 0)
			mat<<= 8;
		else
		{
			res++;
			pos = TZCNT(t);
			inverseSwap(mat, 0xFF, pos);
			t = mat & 0x0101010101010100uLL;
			mat ^= t*(mat & 0xFF);
		}
		if (n > 1)
		{
			t = mat & 0x0202020202020200uLL;
			if (t == 0)
				mat<<= 8;
			else
			{
				res++;
				pos = TZCNT(t) - 1 - 8;
				inverseSwap(mat, 0xFF00, pos);
				t = (mat & 0x0202020202020002uLL) >> 1;
				mat ^= t *((mat & 0xFF00) >> 8);
			}
		}

		if (n > 2)
		{
			t = mat & 0x0404040404040000uLL;
			if (t == 0)
				mat<<= 8;
			else
			{
				res++;
				pos = TZCNT(t) - 2 - 16;
				inverseSwap(mat, 0xFF0000, pos);
				t = (mat & 0x0404040404000404uLL) >> 2;
				mat ^= t *((mat & 0xFF0000) >> 16);
			}
		}

		if (n > 3)
		{
			t = mat & 0x0808080808000000uLL;
			if (t == 0)
				mat<<= 8;
			else
			{
				res++;
				pos = TZCNT(t) - 3 - 24;
				inverseSwap(mat, 0xFF000000, pos);
				t = (mat & 0x0808080800080808uLL) >> 3;
				mat ^= t *((mat & 0xFF000000) >> 24);
			}
		}

		if (n > 4)
		{
			t = mat & 0x1010101000000000uLL;
			if (t == 0)
				mat<<= 8;
			else
			{
				res++;
				pos = TZCNT(t) - 4 - 32;
				inverseSwap(mat, 0xFF00000000, pos);
				t = (mat & 0x1010100010101010uLL) >> 4;
				mat ^= t *((mat & 0xFF00000000) >> 32);
			}
		}
		if (n > 5)
		{
			t=mat & 0x2020200000000000uLL;
			if (t == 0)
				mat <<= 8;
			else
			{
				res++;
				pos = TZCNT(t) - 5 - 40;
				inverseSwap(mat, 0xFF0000000000, pos);
				t = (mat & 0x2020002020202020uLL) >> 5;
				mat ^= t*((mat & 0xFF0000000000) >> 40);
			}
		}
		if (n > 6)
		{
			t = mat & 0x4040000000000000uLL;
			if (t == 0)
				mat <<= 8;
			else
			{
				res++;
				pos = TZCNT(t) - 6 - 48;
				inverseSwap(mat, 0xFF000000000000, pos);
				t = (mat & 0x4000404040404040uLL) >> 6;
				mat ^= t *((mat & 0xFF000000000000) >> 48);
			}
		}
		if (n == 8)
		{
			res += ((mat & 0x8000000000000000uL)!=0);
		}
		return  res;
	}

	template <int k>
	Matrix<m, k> operator*(Matrix <n, k> rhs) const
	{
		T res = ((value & 0x0101010101010101uLL) * 0xFFuLL) & ((rhs.value & 0xFFuLL) * 0x0101010101010101uLL);
		if (n > 1)
			res ^= (((value >> 1) & 0x0101010101010101uLL) * 0xFFuLL) & (((rhs.value >> 8) & 0xFFuLL) * 0x0101010101010101uLL);
		if (n > 2)
			res ^= (((value >> 2) & 0x0101010101010101uLL) * 0xFFuLL) & (((rhs.value >> 16) & 0xFFuLL) * 0x0101010101010101uLL);
		if (n > 3)
			res ^= (((value >> 3) & 0x0101010101010101uLL) * 0xFFuLL) & (((rhs.value >> 24) & 0xFFuLL) * 0x0101010101010101uLL);
		if (n > 4)
			res ^= (((value >> 4) & 0x0101010101010101uLL) * 0xFFuLL) & (((rhs.value >> 32) & 0xFFuLL) * 0x0101010101010101uLL);
		if (n > 5)
			res ^= (((value >> 5) & 0x0101010101010101uLL) * 0xFFuLL) & (((rhs.value >> 40) & 0xFFuLL) * 0x0101010101010101uLL);
		if (n > 6)
			res ^= (((value >> 6) & 0x0101010101010101uLL) * 0xFFuLL) & (((rhs.value >> 48) & 0xFFuLL) * 0x0101010101010101uLL);
		if (n == 8)
			res ^= (((value >> 7) & 0x0101010101010101uLL) * 0xFFuLL) & ((rhs.value >> 56) * 0x0101010101010101uLL);
		return Matrix<m, k>(res);
	}
	Matrix<m, n> operator+(Matrix <m, n> rhs) const
	{
		return Matrix<m, n>(value^rhs.value);
	}
	template <int p,int q>
	Matrix<m+p, n+q> oplus(Matrix <p, q> rhs) const
	{
		return Matrix<m+p, n+q>(value|(rhs.value<<(m*8+n)));
	}
	static Matrix<n, n> identity()
	{
		assert(m == n);
		if (n == 8)
			return Matrix<n, n>(0x8040201008040201uLL);
		return Matrix<n, n>(0x8040201008040201uLL & ((1uLL << (n * 8)) - 1));
	}

	static void preCompRanks()
	{
		ofstream file("ranks.hpp", fstream::app);
		cout << endl << "Computing ranks for m=" <<m<<" n="<<n<< endl;
		if (m == 1 || n == 1)
		{
			file << "template<>" << endl << "int Matrix<" << m << ", " << n << ">::rk() const" << endl << "{" << endl;
			file << "return value!=0;" << endl << "}" << endl;
		}
		else
		{
			file << "template<>" << endl << "int Matrix<" << m << ", " << n << ">::rk() const" << endl << "{" << endl;
			file << "static const char ranks[] = {";
			for (T cur = 0; cur < (1uLL << (m*n)); cur++)
			{
				Matrix < m, n > Cur = getMatrix<m, n>(cur);
				file << Cur.rk();
				if (cur < ((1uLL << (m*n)) - 1))
					file << ",";
			}
			file << "};" << endl << "return ranks[getId()];" << endl << "}" << endl;
		}
	}
};
template <int m,int n>
T getMatrix(T id)
{
	T res = 0;
	for (int j = 0; j < m; j++)
		res |= (((((1uLL << n) - 1) << (j*n)) & id) << (j*(8 - n)));

	return res;
}
void preCompRanks()
{
	Matrix<1, 1>::preCompRanks();
	Matrix<2, 1>::preCompRanks();
	Matrix<1, 2>::preCompRanks();
	Matrix<3, 1>::preCompRanks();
	Matrix<2, 2>::preCompRanks();
	Matrix<1, 3>::preCompRanks();
	Matrix<4, 1>::preCompRanks();
	Matrix<3, 2>::preCompRanks();
	Matrix<2, 3>::preCompRanks();
	Matrix<1, 4>::preCompRanks();
	Matrix<5, 1>::preCompRanks();
	Matrix<4, 2>::preCompRanks();
	Matrix<3, 3>::preCompRanks();
	Matrix<2, 4>::preCompRanks();
	Matrix<1, 5>::preCompRanks();
	Matrix<6, 1>::preCompRanks();
	Matrix<5, 2>::preCompRanks();
	Matrix<4, 3>::preCompRanks();
	Matrix<3, 4>::preCompRanks();
	Matrix<2, 5>::preCompRanks();
	Matrix<1, 6>::preCompRanks();
	Matrix<7, 1>::preCompRanks();
	Matrix<6, 2>::preCompRanks();
	Matrix<5, 3>::preCompRanks();
	Matrix<4, 4>::preCompRanks();
	Matrix<3, 5>::preCompRanks();
	Matrix<2, 6>::preCompRanks();
	Matrix<1, 7>::preCompRanks();
	Matrix<8, 1>::preCompRanks();
	Matrix<7, 2>::preCompRanks();
	Matrix<6, 3>::preCompRanks();
	Matrix<5, 4>::preCompRanks();
	Matrix<4, 5>::preCompRanks();
	Matrix<3, 6>::preCompRanks();
	Matrix<2, 7>::preCompRanks();
	Matrix<1, 8>::preCompRanks();
	Matrix<8, 2>::preCompRanks();
	Matrix<2, 8>::preCompRanks();
}

//#include "ranks.hpp"