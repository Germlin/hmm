#pragma once

#include"hmm.h"
#include"seq.h"
#include"matrix.h"

double forward(HMM& h, Seq& s)
{
	double **alpha = matrix<double>(0, h.N, s.T);

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
				sum += alpha[i][t-1] * h.A[i][j];
			}
			alpha[j][t] = sum * h.B[j][s.O[t]];
		}
	}

	/* 3. Termination */
	double sum = 0;
	for (int i = 0; i < h.N; ++i)
	{
		sum += alpha[i][s.T-1];
	}
	
	freeMatrix<double>(alpha, h.N);
	return sum;
}