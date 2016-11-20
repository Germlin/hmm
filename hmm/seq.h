#pragma once

#include<fstream>

using std::fstream;

class Seq
{
public:
	Seq(fstream&);
	~Seq();
	void print();
	int T;
	int* O;
};