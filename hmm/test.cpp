#include"hmm.h"
#include"seq.h"
#include"forward.h"
#include"backward.h"
#include"viterbi.h"

#include<iostream>

using std::cout;
using std::endl;

int main()
{
	/* Test HMM */
	fstream hf("test.hmm");
	HMM h(hf);
	h.print();

	/* Test Seq */
	fstream sf("test.seq");
	Seq s(sf);
	s.print();

	/* Test forward algorithm */
	cout << forward(h, s) << endl;

	/* Test backward algorithm */
	cout << backward(h, s) << endl;

	/* Test Viterbi algorithm */
	int *path = viterbi(h, s);
	for (int t = 0; t < s.T; ++t)
	{
		cout << path[t] << " ";
	}
	cout << endl;
	return 0;
}