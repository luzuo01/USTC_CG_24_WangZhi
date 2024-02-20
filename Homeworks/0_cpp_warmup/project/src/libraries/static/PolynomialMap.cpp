#include "PolynomialMap.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>

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
	// TODO
}

PolynomialMap PolynomialMap::operator+(const PolynomialMap& right) const {
	auto it1 = m_Polynomial.begin();
    auto it2 = right.m_Polynomial.begin();
    vector<int> deg;
    vector<double> cof;
    while (it1 != m_Polynomial.end() && it2 != right.m_Polynomial.end()) {
        if (it1->first == it2->first) {
            deg.push_back(it1 -> first);
            cof.push_back(it1->second + it2->second);
            ++it1;
            ++it2;
        } else if (it1->first > it2->first) {
            deg.push_back(it2 -> first);
            cof.push_back(it2-> second);
            ++it2;
        } else {
            deg.push_back(it1 -> first);
            cof.push_back(it1->second);
            ++it1;
        }
    }

    while (it1 != m_Polynomial.end()) {
        deg.push_back(it1 -> first);
        cof.push_back(it1-> second);
        ++it1;
    }

    while (it2 != right.m_Polynomial.end()) {
        deg.push_back(it2 -> first);
        cof.push_back(it2->second);
        ++it2;
    }
	PolynomialMap res = PolynomialMap(deg, cof);
	return res;
}

PolynomialMap PolynomialMap::operator-(const PolynomialMap& right) const {
	auto it1 = m_Polynomial.begin();
    auto it2 = right.m_Polynomial.begin();
    vector<int> deg;
    vector<double> cof;
    while (it1 != m_Polynomial.end() && it2 != right.m_Polynomial.end()) {
        if (it1->first == it2->first) {
            deg.push_back(it1 -> first);
            cof.push_back(it1->second - it2->second);
            ++it1;
            ++it2;
        } else if (it1->first > it2->first) {
            deg.push_back(it2 -> first);
            cof.push_back(- it2-> second);
            ++it2;
        } else {
            deg.push_back(it1 -> first);
            cof.push_back(it1->second);
            ++it1;
        }
    }

    while (it1 != m_Polynomial.end()) {
        deg.push_back(it1 -> first);
        cof.push_back(it1-> second);
        ++it1;
    }

    while (it2 != right.m_Polynomial.end()) {
        deg.push_back(it2 -> first);
        cof.push_back(- it2->second);
        ++it2;
    }
	PolynomialMap res = PolynomialMap(deg, cof);
	return res;
}

PolynomialMap PolynomialMap::operator*(const PolynomialMap& right) const {
	PolynomialMap presult;
    vector<int> deg;
    vector<double> cof;

    for (auto term1 : m_Polynomial) {
        for (auto term2 : right.m_Polynomial) {
            int newDegree = term1.first + term2.first;
            double newCoefficient = term1.second * term2.second;

            bool found = false;
            for (auto& resultTerm : presult.m_Polynomial) {
                if (resultTerm.first == newDegree) {
                    resultTerm.second += newCoefficient;
                    found = true;
                    break;
                }
            }

            if (!found) {
                presult.m_Polynomial.insert({newDegree, newCoefficient});
            }
        }
    }
    for (auto it : presult.m_Polynomial) {
        deg.push_back(it.first);
        cof.push_back(it.second);
    }
	PolynomialMap result = PolynomialMap(deg, cof);
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
