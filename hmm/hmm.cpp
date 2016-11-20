#include<string>
#include<iostream>
#include<random>
#include"hmm.h"
#include"matrix.h"

using std::string;
using std::cout;
using std::endl;
using std::random_device;
using std::default_random_engine;
using std::uniform_int_distribution;

HMM::HMM(fstream &f)
{
	string skip;
	f >> skip >> M;	// skip = "M= "

	f >> skip >> N;	// skip = "N= "

	f >> skip;	// skip = "A: "
	A = matrix<double>(0, N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			f >> A[i][j];
		}
	}
	
	f >> skip;	// skip = "B: "
	B = matrix<double>(0, N, M);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			f >> B[i][j];
		}
	}

	f >> skip;	// skip = "pi: "
	pi = new double[N];
	for (int i = 0; i < N; ++i)
	{
		f >> pi[i];
	}
}

HMM::HMM(int n, int m)
{
	N = n;
	M = m;

	random_device rd;
	default_random_engine e(rd());
	uniform_int_distribution<double> u(0, 1);

	A = matrix<double>(0, N, N);
	for (int i = 0; i < N; ++i)
	{
		double sum = 0;
		for (int j = 0; j < N; ++j)
		{
			A[i][j] = u(e);
			sum += A[i][j];
		}
		for (int j = 0; j < N; ++j)
		{
			A[i][j] /= sum;
		}
	}

	B = matrix<double>(0, N, M);
	for (int i = 0; i < N; ++i)
	{
		double sum = 0;
		for (int j = 0; j < M; ++j)
		{
			B[i][j] = u(e);
			sum += B[i][j];
		}
		for (int j = 0; j < M; ++j)
		{
			B[i][j] /= sum;
		}
	}

	pi = new double[N];
	for (int i = 0; i < N; ++i)
	{
		pi[i] = u(e);
	}
}

HMM::~HMM()
{
	freeMatrix<double>(A, N);
	freeMatrix<double>(B, N);
	delete[] pi;
}

void HMM::print()
{
	cout << "M= " << M << endl;
	cout << "N= " << N << endl;
	cout << "A: " << endl;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			cout << A[i][j] << "";
		}
		cout << endl;
	}

	cout << "B: " << endl;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			cout << B[i][j] << " ";
		}
		cout << endl;
	}

	cout << "pi: " << endl;
	for (int i = 0; i < N; ++i)
	{
		cout << pi[i] << " ";
	}
	cout << endl;
}