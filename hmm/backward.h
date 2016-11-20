#pragma once

#include"hmm.h"
#include"seq.h"
#include"matrix.h"

double backward(HMM &h, Seq &s)
{
	double **beta = matrix<double>(0, h.N, s.T);

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

	/* 3. Termination*/
	double sum = 0;
	for (int i = 0; i < h.N; ++i)
	{
		sum += h.pi[i] * h.B[i][s.O[0]] * beta[i][0];
	}

	freeMatrix<double>(beta, h.N);
	return sum;
}