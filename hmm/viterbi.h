#pragma once

#include"hmm.h"
#include"seq.h"
#include"matrix.h"

int* viterbi(HMM& h, Seq& s)
{
	double **delta = matrix<double>(0, h.N, s.T);
	int **psi = matrix<int>(0, h.N, s.T);
	int *path = new int[s.T];

	/* 1. Initialization */
	for (int i = 0; i < h.N; ++i)
	{
		delta[i][0] = h.pi[i] * h.B[i][s.O[0]];
		psi[i][0] = 0;
	}

	/* 2. Recursion */
	for (int t = 1; t < s.T; ++t)
	{
		for (int j = 0; j < h.N; ++j)
		{
			double max = 0;
			int max_node = 0;
			for (int i = 0; i < h.N; ++i)
			{
				if (max < delta[i][t - 1] * h.A[i][j])
				{
					max = delta[i][t - 1] * h.A[i][j];
					max_node = i;
				}
			}
			delta[j][t] = max*h.B[j][s.O[t]];
			psi[j][t] = max_node;
		}
	}

	/* 3. Termination */
	double P = 0;
	for (int i = 0; i < h.N; ++i)
	{
		if (P < delta[i][s.T - 1])
		{
			P = delta[i][s.T - 1];
			path[s.T - 1] = i;
		}
	}

	/* 4. Path*/
	for (int t = s.T - 2; t >= 0; --t)
	{
		path[t] = psi[path[t + 1]][t + 1];
	}

	freeMatrix<double>(delta, h.N);
	freeMatrix<int>(psi, h.N);
	return path;
}