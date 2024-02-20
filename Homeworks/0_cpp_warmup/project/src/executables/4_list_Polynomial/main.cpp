#include "PolynomialList.h"

#include <list>
#include <vector>
#include "PolynomialList.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include <cassert>

using namespace std;

int main(int argc, char** argv) {
	PolynomialList p1("../data/P4.txt");
	PolynomialList p2("../data/P2.txt");
	PolynomialList p3;
	p1.Print();
	p2.Print();

	p3 = p1 + p2;
	p3.Print();
	p3 = p1 - p2;
	p3.Print();

	p3 = p1 * p2;
	p3.Print();	
	return 0;
}