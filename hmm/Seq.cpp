#include"seq.h"
#include<string>
#include<iostream>

using std::cout;
using std::endl;
using std::string;

Seq::Seq(fstream& f)
{
	string skip;
	f >> skip >> T;	// skip = "T= "
	O = new int[T];
	for (int i = 0; i < T; ++i)
	{
		f >> O[i];
	}
}

Seq::~Seq()
{
	delete[] O;
}

void Seq::print()
{
	cout << "T= " << T << endl;
	for (int i = 0; i < T; ++i)
	{
		cout << O[i] << " ";
	}
	cout << endl;
}