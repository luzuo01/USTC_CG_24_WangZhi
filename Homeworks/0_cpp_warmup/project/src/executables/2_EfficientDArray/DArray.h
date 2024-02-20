#pragma once

#include <iostream>
#include <assert.h>
using namespace std;


template <class T>
class DArray {
public:
	DArray(); // default constructor
	DArray(int nSize, T dValue = 0); // set an array with default values
	DArray(const DArray& arr); // copy constructor
	~DArray(); // deconstructor

	void Print() const; // print the elements of the array

	int GetSize() const; // get the size of the array
	void SetSize(int nSize); // set the size of the array

	const T& GetAt(int nIndex) const; // get an element at an index
	void SetAt(int nIndex, T dValue); // set the value of an element

	T& operator[](int nIndex); // overload operator '[]'
	const T& operator[](int nIndex) const; // overload operator '[]'

	void PushBack(T dValue); // add a new element at the end of the array
	void DeleteAt(int nIndex); // delete an element at some index
	void InsertAt(int nIndex, T dValue); // insert a new element at some index

	DArray& operator = (const DArray& arr); //overload operator '='

private:
	T* m_pData; // the pointer to the array memory
	int m_nSize; // the size of the array
	int m_nMax;

private:
	void Init(); // initilize the array
	void Free(); // free the array
	void Reserve(int nSize); // allocate enough memory
};

template<typename T>
// default constructor
DArray<T>::DArray() {
	Init();
}

template<typename T>
// set an array with default values
DArray<T>::DArray(int nSize, T dValue) : 
m_pData(new T [nSize]), m_nMax(nSize), m_nSize(nSize)
{
	for (int i = 0; i < nSize; i++) {
		m_pData[i] = dValue;
	}
}

template<typename T>
DArray<T>::DArray(const DArray& arr)
	: m_pData(new T [arr.m_nSize]), m_nSize(arr.m_nSize), m_nMax(arr.m_nSize)
{
	for (int i = 0; i < m_nSize; i++)
		m_pData[i] = arr.m_pData[i];
}

template<typename T>
// deconstructor
DArray<T>::~DArray() {
	Free();
}

template<typename T>
// display the elements of the array
void DArray<T>::Print() const {
	cout << "Length: " << m_nSize << " = [";
	for (int i = 0; i < m_nSize; i++)
	{
		cout << GetAt(i) << " ";
	}
	cout << "]" << endl;
}

template<typename T>
// initilize the array
void DArray<T>::Init() {
	m_nSize = 0;
	m_nMax = 0;
	m_pData = nullptr;
}

template<typename T>
// free the array
void DArray<T>::Free() {
	delete[] m_pData;
	m_pData = nullptr;
	m_nSize = 0;
	m_nMax = 0;
}

template<typename T>
// get the size of the array
int DArray<T>::GetSize() const {
	return m_nSize;
}


template<typename T>
void DArray<T>::Reserve(int nSize) {
	if (nSize <= m_nMax) {
		return;
	}
	if (m_nMax == 0) {
		m_nMax = 1;
	}
	while (m_nMax < nSize)
	{
		m_nMax *= 2;
	}
	T* temp = new T[m_nMax];
	for (int i = 0; i < m_nSize; i++) {
		temp[i] = m_pData[i];
	}
	for (int i = m_nSize; i < m_nMax; i++) {
		temp[i] = 0;
	}
	delete[] m_pData;
	m_pData = temp;
}

template<typename T>
// set the size of the array
void DArray<T>::SetSize(int nSize) {
	if (nSize <= m_nSize) {
		return;
	}
	Reserve(nSize);
	m_nSize = nSize;
}

template<typename T>
// get an element at an index
const T& DArray<T>::GetAt(int nIndex) const {
	assert(nIndex >= 0 && nIndex < m_nSize);
	return m_pData[nIndex];
}

template<typename T>
// set the value of an element 
void DArray<T>::SetAt(int nIndex, T dValue) {
	assert(nIndex >= 0 && nIndex < m_nSize);
	m_pData[nIndex] = dValue;
}

template<typename T>
// overload operator '[]'
T& DArray<T>::operator[](int nIndex) {
	assert(nIndex >= 0 && nIndex < m_nSize);
	return m_pData[nIndex];
}

template<typename T>
// overload operator '[]'
const T& DArray<T>::operator[](int nIndex) const {
	assert(nIndex >= 0 && nIndex < m_nSize);
	return m_pData[nIndex];
}

template<typename T>
// add a new element at the end of the array
void DArray<T>::PushBack(T dValue) {
	Reserve(m_nSize + 1);
	m_pData[m_nSize] = dValue;
	m_nSize++;
}

template<typename T>
// delete an element at some index
void DArray<T>::DeleteAt(int nIndex) {
	assert(nIndex >= 0 && nIndex < m_nSize);
	for (int i = nIndex; i < m_nSize - 1; i++)
	{
		m_pData[i] = m_pData[i + 1];
	}
	m_nSize--;
}

template<typename T>
// insert a new element at some index
void DArray<T>::InsertAt(int nIndex, T dValue) {
	assert(nIndex >= 0 && nIndex <= m_nSize);
	Reserve(m_nSize + 1);

	for (int i = m_nSize; i > nIndex ; i--)
	{
		m_pData[i] = m_pData[i - 1];
	}
	m_pData[nIndex] = dValue;
	m_nSize++;
}

template<typename T>
// overload operator '='
DArray<T>& DArray<T>::operator = (const DArray& arr) {
	delete[] m_pData;
	m_nSize = arr.m_nSize;
	m_pData = new T[m_nSize];
	Reserve(arr.m_nSize);
	for (int i = 0; i < m_nSize; i++){
		m_pData[i] = arr[i];
	}
	return *this;
}