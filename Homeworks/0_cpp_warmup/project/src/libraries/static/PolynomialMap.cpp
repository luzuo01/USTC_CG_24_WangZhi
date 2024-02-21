#include "PolynomialMap.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>
#define EPSILON 1.0e-13

using namespace std;

PolynomialMap::PolynomialMap(const PolynomialMap& other): m_Polynomial(other.m_Polynomial){}

PolynomialMap::PolynomialMap(const string& file) {
    ReadFromFile(file);
}

PolynomialMap::PolynomialMap(const double* cof, const int* deg, int n) {
	for (int i = 0; i < n; i++)
    {
        m_Polynomial.insert({deg[i], cof[i]});
    }
}

PolynomialMap::PolynomialMap(const vector<int>& deg, const vector<double>& cof) {
	assert(deg.size() == cof.size());
	for (int i = 0; i < deg.size(); i++)
    {
        m_Polynomial.insert({deg[i], cof[i]});
    }
}

double PolynomialMap::coff(int i) const {
	assert (i >= 0 && i < m_Polynomial.size());
	return m_Polynomial.at(i); // you should return a correct value
}

double& PolynomialMap::coff(int i) {
	assert (i >= 0 && i < m_Polynomial.size());
	return m_Polynomial.at(i);
}

void PolynomialMap::compress() {
	map<int, double> temp = m_Polynomial;
    m_Polynomial.clear();
    for (auto term : temp) {
        if (fabs(term.second) > EPSILON) {
            temp[term.first] = term.second; 
        } 
    }
    m_Polynomial = temp;
}

PolynomialMap PolynomialMap::operator+(const PolynomialMap& right) const {
	PolynomialMap temp = right;
    for (auto it : m_Polynomial)
    {
        temp.m_Polynomial[it.first] += it.second;
    }
    temp.compress();
    return temp;
}

PolynomialMap PolynomialMap::operator-(const PolynomialMap& right) const {
    PolynomialMap temp = *this;
    for (auto it : right.m_Polynomial)
    {
        temp.m_Polynomial[it.first] -= it.second;
    }
    temp.compress();
    return temp;
}
PolynomialMap PolynomialMap::operator*(const PolynomialMap& right) const {
	PolynomialMap result;
    for (auto term1 : m_Polynomial)
    {
        for (auto term2 : right.m_Polynomial)
        {
            int deg = term1.first + term2.first;
            double cof = term1.second * term2.second;
            result.m_Polynomial[deg] = cof;
        }
    }
    result.compress();
    return result;
}

PolynomialMap& PolynomialMap::operator=(const PolynomialMap& right) {
	if (this != &right) {
        m_Polynomial = right.m_Polynomial;
    }
    return *this;
}

void PolynomialMap::Print() const {
	std::cout << "f(x) = ";
    auto it = m_Polynomial.begin();
    while (it != m_Polynomial.end()) {
        if (std::next(it) != m_Polynomial.end() && next(it)->second > 0){
            if (it -> first == 0)
            {
                std::cout << it -> second << " + ";
            } else {
                std::cout << it -> second << "x^" << it -> first << " + ";
            }
        } else if (std::next(it) != m_Polynomial.end() && next(it)-> second <= 0) {
            if (it -> first == 0)
            {
                std::cout << it -> second << ' ';
            } else {
                std::cout << it -> second << "x^" << it -> first << ' ';
            }
        } else if (std::next(it) == m_Polynomial.end() && it -> second != 0){
            std::cout << it -> second << "x^" << it -> first;
        }
        it++;

    }
    std::cout << endl;
}

bool PolynomialMap::ReadFromFile(const string& file) {
    vector<int> deg_vec;
    vector<double> cof_vec;
    ifstream infile(file);
    if (!infile) {
        return false;
    }

    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        char type;
        int length;
        int deg;
        double cof;

        // Read the first character (type) from the line
        if (!(iss >> type)) {
            return false;
        }

        if (type == 'P') {
            continue;
        } else {
            iss.clear();
            iss.seekg(0);
            if (!(iss >> deg >> cof)) {
                return false;
            }
            bool marker = true;
            for (int i = 0; i < deg_vec.size(); i++) {
                if (deg_vec[i] == deg)
                {
                    cof_vec[i] += cof;
                    marker = false;
                } 
            }
            if (marker) {
                deg_vec.push_back(deg);
                cof_vec.push_back(cof);
            }
        }
    }
 
    PolynomialMap result = PolynomialMap(deg_vec, cof_vec);
    m_Polynomial = result.m_Polynomial;
    infile.close();
    return true;
}
