#include"baumwelch.h"
#include"matrix.h"

void computeAlpha(HMM& h, Seq& s, double **alpha)
{
	/* 1. Initialization */
	for (int i = 0; i < h.N; ++i)
	{
		alpha[i][0] = h.pi[0] * h.B[i][s.O[0]];
	}

	/* 2. Induction */
	for (int t = 1; t < s.T; ++t)
	{
		for (int j = 0; j < h.N; ++j)
		{
			double sum = 0;
			for (int i = 0; i < h.N; ++i)
			{
				sum += alpha[i][t - 1] * h.A[i][j];
			}
			alpha[j][t] = sum * h.B[j][s.O[t]];
		}
	}
}

void computeBeta(HMM& h, Seq& s, double **beta)
{
	/* 1. Initialization */
	for (int i = 0; i < h.N; ++i)
	{
		beta[i][s.T - 1] = 1;
	}

	/* 2. Induction */
	for (int t = s.T - 2; t >= 0; --t)
	{
		for (int i = 0; i < h.N; ++i)
		{
			double sum = 0;
			for (int j = 0; j < h.N; ++j)
			{
				sum += h.A[i][j] * h.B[j][s.O[t + 1]] * beta[j][t + 1];
			}
			beta[i][t] = sum;
		}
	}
}

void computeXi(HMM &h, Seq &s, double **alpha, double **beta, double ***xi)
{
	for (int t = 0; t < s.T; ++t)
	{
		double sum = 0;
		for (int i = 0; i < h.N; ++i)
		{
			for (int j = 0; j < h.N; ++j)
			{
				xi[t][i][j] = alpha[i][t] * beta[j][t + 1] * h.A[i][j] * h.B[j][s.O[t + 1]];
				sum += xi[t][i][j];
			}
		}
		for (int i = 0; i < h.N; ++i)
		{
			for (int j = 0; j < h.N; ++j)
			{
				xi[t][i][j] /= sum;
			}
		}
	}
}

void computeGamma(HMM &h, Seq &s, double **alpha, double **beta, double **gamma)
{
	for (int t = 0; t < s.T; ++t)
	{
		double sum = 0;
		for (int i = 0; i < h.N; ++i)
		{
			gamma[i][t] = alpha[i][t] * beta[i][t];
			sum += gamma[i][t];
		}
		for (int i = 0; i < h.N; ++i)
		{
			gamma[i][t] /= sum;
		}
	}
}

void baumWelch(HMM &h, Seq &s, double DELTA)
{
	double **gamma = matrix<double>(0, h.N, s.T);
	double **alpha = matrix<double>(0, h.N, s.T);
	double **beta = matrix<double>(0, h.N, s.T);
	double ***xi = new double**[s.T];
	for (int i = 0; i < s.T; ++i)
	{
		xi[i] = new double*[h.N];
		for (int j = 0; j < h.N; ++j)
		{
			xi[i][j] = new double[h.N];
		}
	}

	computeAlpha(h, s, alpha);
	computeBeta(h, s, beta);
	computeGamma(h, s, alpha, beta, gamma);
	computeXi(h, s, alpha, beta, xi);

	/* free memeroy */
	for (int i = 0; i < s.T; ++i)
	{
		for (int j = 0; j < h.N; ++j)
		{
			delete[] xi[i][j];
		}
		delete[] xi[i];
	}
	delete[] xi;
	freeMatrix(gamma, h.N);
	freeMatrix(alpha, h.N);
	freeMatrix(beta, h.N);
}