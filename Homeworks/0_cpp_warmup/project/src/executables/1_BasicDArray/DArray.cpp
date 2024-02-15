// implementation of class DArray
#include <iostream>

#include "DArray.h"

#include <cassert>

using namespace std;

// default constructor
DArray::DArray() {
	Init();
}

// set an array with default values
DArray::DArray(int nSize, double dValue)
	: m_pData(new double[nSize]), m_nSize(nSize)
{
	for (int i = 0; i < m_nSize; i++) {
		m_pData[i] = dValue;
	}
}

DArray::DArray(const DArray& arr)
	: m_pData(new double[arr.m_nSize]), m_nSize(arr.m_nSize)
{
	for (int i = 0; i < m_nSize; i++) {
		m_pData[i] = arr[i];
	}
}

// deconstructor
DArray::~DArray() {
	Free();
}

// display the elements of the array
void DArray::Print() const {
	cout << "Length: " << m_nSize << " = [";
	for (int i = 0; i < m_nSize; i++)
	{
		cout << GetAt(i) << " ";
	}
	cout << "]" << endl;
}

// initilize the array
void DArray::Init() {
	m_nSize = 0;
	m_pData = nullptr;
}

// free the array
void DArray::Free() {
	delete[] m_pData;
	m_pData = nullptr;
	m_nSize = 0;
}

// get the size of the array
int DArray::GetSize() const {
	return m_nSize;
}

// set the size of the array
void DArray::SetSize(int nSize) {
	if (nSize == m_nSize) {
		return;
	}
	double* temp = new double[nSize];
	int min = nSize < m_nSize ? nSize : m_nSize;
	for (int i = 0; i < min; i++) {
		temp[i] = m_pData[i];
	}
	for (int i = min; i < nSize; i++) {
		temp[i] = 0.0;
	}
	delete[] m_pData;
	m_pData = temp;
	m_nSize = nSize;
}

// get an element at an index
const double& DArray::GetAt(int nIndex) const {
	assert (nIndex >= 0 && nIndex < m_nSize);
	return m_pData[nIndex];
}

// set the value of an element 
void DArray::SetAt(int nIndex, double dValue) {
	assert (nIndex >= 0 && nIndex < m_nSize);
	m_pData[nIndex] = dValue;
}

double& DArray::operator[](int nIndex) {
	assert(nIndex >= 0 && nIndex < m_nSize);
	return m_pData[nIndex];
}

// overload operator '[]'
const double& DArray::operator[](int nIndex) const {
	assert(nIndex >= 0 && nIndex < m_nSize);
	return m_pData[nIndex];
}

// add a new element at the end of the array
void DArray::PushBack(double dValue) {
	double* temp = new double[m_nSize + 1];
	for (int i = 0; i < m_nSize; i++) {
		temp[i] = m_pData[i];
	}
	temp[m_nSize] = dValue;
	delete[] m_pData;
	m_pData = temp;
	m_nSize++;
}

// delete an element at some index
void DArray::DeleteAt(int nIndex) {
	assert(nIndex >= 0 && nIndex < m_nSize);

	double* temp = new double[(m_nSize) - 1];

	for (int i = 0; i < nIndex; i++)
		temp[i] = m_pData[i];

	for (int i = nIndex; i < m_nSize - 1; i++)
		temp[i] = m_pData[i + 1];

	delete[] m_pData;
	m_pData = temp;
	m_nSize--;
}

// insert a new element at some index
void DArray::InsertAt(int nIndex, double dValue) {
	assert(nIndex >= 0 && nIndex <= m_nSize);

	double* temp = new double[(m_nSize) + 1];

	for (int i = 0; i < nIndex; i++)
		temp[i] = m_pData[i];

	temp[nIndex] = dValue;

	for (int i = nIndex + 1; i < m_nSize + 1; i++)
		temp[i] = m_pData[i - 1];

	delete[] m_pData;
	m_pData = temp;
	m_nSize++;
}

// overload operator '='
DArray& DArray::operator = (const DArray& arr) {
	delete[] m_pData;
	m_nSize = arr.m_nSize;
	m_pData = new double[m_nSize];
	for (int i = 0; i < m_nSize; i++){
		m_pData[i] = arr[i];
	}
	return *this;
}
