#include "PolynomialList.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <cassert>

using namespace std;

PolynomialList::PolynomialList(const PolynomialList& other): m_Polynomial(other.m_Polynomial)
{}

PolynomialList::PolynomialList(const string& file) {
    ReadFromFile(file);
}

PolynomialList::PolynomialList(const double* cof, const int* deg, int n) {
    for (int i = 0; i < n; i++)
    {
        m_Polynomial.push_back(Term(deg[i], cof[i]));
    }
}

PolynomialList::PolynomialList(const vector<int>& deg, const vector<double>& cof) {
    for (int i = 0; i < deg.size(); i++)
    {
        m_Polynomial.push_back(Term(deg[i], cof[i]));
    }
}

double PolynomialList::coff(int i) const {
    auto seer = m_Polynomial.begin();
    advance(seer, i);
    return seer -> cof; // you should return a correct value
}

double& PolynomialList::coff(int i) {
    auto seer = m_Polynomial.begin();
    advance(seer, i);
    return seer -> cof;
}

void PolynomialList::compress() {
    // TODO
}

PolynomialList PolynomialList::operator+(const PolynomialList& right) const {
    auto it1 = m_Polynomial.begin();
    auto it2 = right.m_Polynomial.begin();
    vector<int> deg;
    vector<double> cof;
    while (it1 != m_Polynomial.end() && it2 != right.m_Polynomial.end()) {
        if (it1->deg == it2->deg) {
            deg.push_back(it1 -> deg);
            cof.push_back(it1->cof + it2->cof);
            ++it1;
            ++it2;
        } else if (it1->deg > it2->deg) {
            deg.push_back(it2 -> deg);
            cof.push_back(it2-> cof);
            ++it2;
        } else {
            deg.push_back(it1 -> deg);
            cof.push_back(it1->cof);
            ++it1;
        }
    }

    while (it1 != m_Polynomial.end()) {
        deg.push_back(it1 -> deg);
        cof.push_back(it1-> cof);
        ++it1;
    }

    while (it2 != right.m_Polynomial.end()) {
        deg.push_back(it2 -> deg);
        cof.push_back(it2->cof);
        ++it2;
    }

    int n = deg.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (deg[j] > deg[j + 1]) {
                int temp = deg[j];
                deg[j] = deg[j + 1];
                deg[j + 1] = temp;

                int temp2 = cof[j];
                cof[j] = cof[j + 1];
                cof[j + 1] = temp2;
            }
        }
    }

    PolynomialList result = PolynomialList(deg, cof);
    return result;
}

PolynomialList PolynomialList::operator-(const PolynomialList& right) const {
    auto it1 = m_Polynomial.begin();
    auto it2 = right.m_Polynomial.begin();
    vector<int> deg;
    vector<double> cof;
    while (it1 != m_Polynomial.end() && it2 != right.m_Polynomial.end()) {
        if (it1->deg == it2->deg) {
            deg.push_back(it1 -> deg);
            cof.push_back(it1->cof - it2->cof);
            ++it1;
            ++it2;
        } else if (it1->deg > it2->deg) {
            deg.push_back(it2 -> deg);
            cof.push_back(- it2-> cof);
            ++it2;
        } else {
            deg.push_back(it1 -> deg);
            cof.push_back(it1->cof);
            ++it1;
        }
    }

    while (it1 != m_Polynomial.end()) {
        deg.push_back(it1 -> deg);
        cof.push_back(it1-> cof);
        ++it1;
    }

    while (it2 != right.m_Polynomial.end()) {
        deg.push_back(it2 -> deg);
        cof.push_back(- it2->cof);
        ++it2;
    }

    int n = deg.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (deg[j] > deg[j + 1]) {
                int temp = deg[j];
                deg[j] = deg[j + 1];
                deg[j + 1] = temp;

                int temp2 = cof[j];
                cof[j] = cof[j + 1];
                cof[j + 1] = temp2;
            }
        }
    }

    PolynomialList result = PolynomialList(deg, cof);
    return result;
}

PolynomialList PolynomialList::operator*(const PolynomialList& right) const {
    PolynomialList presult;
    vector<int> deg;
    vector<double> cof;

    for (auto term1 : m_Polynomial) {
        for (auto term2 : right.m_Polynomial) {
            int newDegree = term1.deg + term2.deg;
            double newCoefficient = term1.cof * term2.cof;

            bool found = false;
            for (auto& resultTerm : presult.m_Polynomial) {
                if (resultTerm.deg == newDegree) {
                    resultTerm.cof += newCoefficient;
                    found = true;
                    break;
                }
            }

            if (!found) {
                presult.m_Polynomial.push_back({newDegree, newCoefficient});
            }
        }
    }
    for (auto it : presult.m_Polynomial) {
        deg.push_back(it.deg);
        cof.push_back(it.cof);
    }

    int n = deg.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (deg[j] > deg[j + 1]) {
                int temp = deg[j];
                deg[j] = deg[j + 1];
                deg[j + 1] = temp;

                int temp2 = cof[j];
                cof[j] = cof[j + 1];
                cof[j + 1] = temp2;
            }
        }
    }
    PolynomialList result = PolynomialList(deg, cof);
    return result;
}

PolynomialList& PolynomialList::operator=(const PolynomialList& right) {
    if (this != &right) {
        m_Polynomial = right.m_Polynomial;
    }
    return *this;
}

void PolynomialList::Print() const {
    std::cout << "f(x) = ";
    auto it = m_Polynomial.begin();
    while (it != m_Polynomial.end()) {
        if (std::next(it) != m_Polynomial.end() && next(it)->cof > 0){
            if (it -> deg == 0)
            {
                std::cout << it -> cof << " + ";
            } else {
                std::cout << it -> cof << "x^" << it -> deg << " + ";
            }
        } else if (std::next(it) != m_Polynomial.end() && next(it)->cof <= 0) {
            if (it -> deg == 0)
            {
                std::cout << it -> cof << ' ';
            } else {
                std::cout << it -> cof << "x^" << it -> deg << ' ';
            }
        } else if (std::next(it) == m_Polynomial.end() && it -> cof != 0){
            std::cout << it -> cof << "x^" << it -> deg;
        }
        it++;

    }
    std::cout << endl;
}

bool PolynomialList::ReadFromFile(const string& file) {
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
    //let's just use a bubble sort to sort these 2 at once.
    int n = deg_vec.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (deg_vec[j] > deg_vec[j + 1]) {
                int temp = deg_vec[j];
                deg_vec[j] = deg_vec[j + 1];
                deg_vec[j + 1] = temp;

                int temp2 = cof_vec[j];
                cof_vec[j] = cof_vec[j + 1];
                cof_vec[j + 1] = temp2;
            }
        }
    }
    PolynomialList result = PolynomialList(deg_vec, cof_vec);
    m_Polynomial = result.m_Polynomial;
    infile.close();
    return true;
}

PolynomialList::Term& PolynomialList::AddOneTerm(const Term& term) {
    m_Polynomial.push_back(term);
    return m_Polynomial.back();
}
