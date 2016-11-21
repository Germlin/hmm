#pragma once

#include<fstream>
#include"seq.h"

using std::fstream;

class HMM
{
public:
	HMM(fstream&);	// initialize the hmm model from a file
	HMM(int, int);	// initialize a random hmm model by specifying N and M
	~HMM();
	void print();
	friend double forward(HMM&, Seq&);
	friend double backward(HMM&, Seq&);
	friend int* viterbi(HMM&, Seq&);
	/* following function is used for Baum-Welch algorithm */
	friend void computeAlpha(HMM&, Seq&, double**);
	friend void computeBeta(HMM&, Seq&, double**);
	friend void computeXi(HMM&, Seq&, double**, double**, double***);
	friend void computeGamma(HMM&, Seq&, double**, double**, double**);
	friend void baumWelch(HMM&, Seq&, double);
private:
	int N;	// number of states; Q={1,2,...,N}
	int M;	// number of observation symbols; V={1,2,...,M}
	double **A;	/* A[1..N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j
			   at time t+1 */
	double **B;	/* B[1..N][1..M]. b[j][k] is the probability of
			   of observing symbol k in state j */
	double *pi;	/* pi[1..N] pi[i] is the initial state distribution. */
};