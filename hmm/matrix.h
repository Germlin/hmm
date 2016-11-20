#pragma once

template<typename T>
T** matrix(T u, int n, int m)
{
	T **mt = new T*[n];
	for (int i = 0; i < n; ++i)
	{
		mt[i] = new T[m];
		for (int j = 0; j < m; ++j)
		{
			mt[i][j] = u;
		}
	}
	return mt;
}

template<typename T>
void freeMatrix(T** mt, int n)
{
	for (int i = 0; i < n; ++i)
	{
		delete[] mt[i];
	}
	delete[] mt;
}