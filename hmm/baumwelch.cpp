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

double computeProb(double **alpha, int N, int T)
{
	double sum = 0;
	for (int i = 0; i < N; ++i)
	{
		sum += alpha[i][T - 1];
	}
	return sum;
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

	double probprev = computeProb(alpha, h.N, s.T);
	double prob = 0;
	double delta = 0;

	do
	{
		/* update A */
		for (int i = 0; i < h.N; ++i)
		{
			for (int j = 0; j < h.N; ++j)
			{
				double numerator = 0;
				double denominator = 0;
				for (int t = 0; t < s.T - 1; ++t)
				{
					numerator += xi[t][i][j];
					denominator += gamma[i][t];
				}
				h.A[i][j] = numerator / denominator;
			}
		}

		/* update B */
		for (int j = 0; j < h.N; ++j)
		{
			for (int k = 0; k < h.M; ++k)
			{
				double numerator = 0;
				double denominator = 0;
				for (int t = 0; t < s.T; ++t)
				{
					if (s.O[t] == k)
					{
						numerator += gamma[j][t];
					}
					denominator += gamma[j][t];
				}
				h.B[j][k] = numerator / denominator;
			}
		}

		/* update pi */
		for (int i = 0; i < h.N; ++i)
		{
			h.pi[i] = gamma[i][1];
		}

		computeAlpha(h, s, alpha);
		computeBeta(h, s, beta);
		computeGamma(h, s, alpha, beta, gamma);
		computeXi(h, s, alpha, beta, xi);

		prob = computeProb(alpha, h.N, s.T);
		delta = prob - probprev;
		probprev = prob;

	} while (delta > DELTA);

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