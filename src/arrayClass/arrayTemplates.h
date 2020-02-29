#include "arrayClass.h"

//------------------------------ CONSTRUCTORS ----------------------------//
template <class T>
array<T>::array(unsigned int _size): sizeArray(_size)
{
  this->tab = new T[_size]; // POINTER
  for (unsigned int i = 0; i < _size; ++i)
  {
    this->tab[i] = 0.; // assignement of the value 0
  }

};

template <class T>
array<T>::array(unsigned int _size, const T& value): sizeArray(_size) // constructor with value given
{
  this->tab = new T[_size]; // POINTER
  for (unsigned int i = 0; i < _size; ++i)
  {
    this->tab[i] = value; // assignement of the value 'value'
  }
};

template <class T>
array<T>::array(const array<T> &object): sizeArray(object.sizeArray) // copy constructor
{
  this->tab = new T[object.sizeArray];
  for (unsigned int i = 0; i < object.sizeArray; ++i)
  {
    this->tab[i] = object.tab[i]; // copy the right value in the new array
  }
};

// ---------------------------- DESTRUCTOR -------------------------------//
template <class T>
array<T>::~array()
{
  delete[] this->tab; // delete the pointer
};

// ---------------------------- FUNCTIONS --------------------------------//
template <class T>
unsigned int array<T>::size() // get the size of the data table
{
  return this->sizeArray;
};

template <class T>
array<T>&  array<T>::init(const T& value) // initialization with th value
{
  for (unsigned int i = 0; i < this->sizeArray; ++i)
  {
    this->tab[i] = value;
  }
  return *this;
};

template <class T>
array<T>& array<T>::resize(unsigned int _size) // resize and put 0 for the data value
{
  delete[] this->tab;;
  this->tab = new T[_size];
  this->sizeArray = _size;
  for (unsigned int i = 0; i < _size; ++i)
  {
    this->tab[i] = 0;
  }
  return *this;
};

template <class T>
array<T>& array<T>::resize(unsigned int _size, const T& value) // resize and put 'value' for the data value
{
  delete[] this->tab; // delete
  this->tab = new T[_size]; // build a new
  this->sizeArray = _size;
  for (unsigned int i = 0; i < _size; ++i)
  {
    this->tab[i] = value;
  }
  return *this;
};

// ----------------------- OPERATORS -------------------------------------//
template <class T>
T& array<T>::operator[](unsigned int i) // overload operator [] to access to an element
{
  assert(i < this->sizeArray && "INDEX_OUT_OF_BOUNDS");
  return this->tab[i];
};

template <class T>
array<T>& array<T>::operator=(const array<T>& object)
{
  delete[] this->tab;
  this->tab = new T[object.sizeArray];
  this->sizeArray = object.sizeArray;

  for (unsigned int i = 0; i < object.sizeArray; ++i)
  {
    this->tab[i] = object.tab[i];
  }
  return *this;
};

// ----------------------- FRIEND OPERATORS --------------------------------//
template <class T>
array<T> operator+(const array<T>& object1, const array<T>& object2)
{
  assert(object1.sizeArray == object2.sizeArray && "INCOMPATIBLE_SIZE_OPE+");

  array<T> RTN(object1.sizeArray);
  for (unsigned int i = 0; i < RTN.size(); ++i)
  {
    RTN[i] = object1.tab[i] + object2.tab[i];
  }
  return RTN;
};

template <class T>
array<T> operator-(const array<T>& object1, const array<T>& object2)
{
  assert(object1.sizeArray == object2.sizeArray && "INCOMPATIBLE_SIZE_OPE+");

  array<T> RTN(object1.sizeArray);
  for (unsigned int i = 0; i < RTN.size(); ++i)
  {
    RTN[i] = object1.tab[i] - object2.tab[i];
  }
  return RTN;
};

template <class T>
T operator|(const array<T>& object1, const array<T>& object2)
{
  assert(object1.sizeArray == object2.sizeArray && "INCOMPATIBLE_SIZE_OPE+");

  T RTN(0.);
  for (unsigned int i = 0; i < object1.sizeArray; ++i)
  {
    RTN += object1.tab[i] * object2.tab[i];
  }
  return RTN;
};

template <class T, class U>
array<T> operator*(const U& alpha, const array<T>& object)
{
  array<T> RTN(object.sizeArray);

  for (unsigned int i = 0; i < object.sizeArray; ++i)
  {
    RTN[i] = object.tab[i] * alpha;
  }
  return RTN;
};

template <class T>
ostream& operator<<(ostream& flux, const array<T> &object)
{
  for (unsigned int i = 0; i < object.sizeArray; ++i)
  {
    flux << object.tab[i] << " ";
  }
  flux << endl;
  return flux;
}
