/* C'est ma propre classe 'array', elle utilise des tableaux à la C comme c'était le cas dans le programme
   de départ. Constructeur et destructeur de pointeurs obligent, j'ai passé mon programme avec Valgrind
   pour vérifier s'il n'existerait pas de fuites de mémoire, normalement il ne devrait pas y avoir de problème
   de ce côté ci */


#ifndef ARRAY_CLASS
#define ARRAY_CLASS

#include <cassert> // for assert :: error of using
#include <iostream> // stream :: output and input
#include <stdlib.h> // for new[] and delete[]

using namespace std;

/* class array template :
  Be carefull matrix or array more than 1 dimension are not supported */

template <class T>
class array
{
private:
  unsigned int sizeArray;
  T* tab; // data table
public:
  array(unsigned int _size); // constructor
  array(unsigned int _size, const T& value); // constructor with initialization by a value (type template obviously)
  array(const array<T> &object); // copy constructor
  ~array(); // destructor
  unsigned int size(); // get the size of data array
  array<T>&  init(const T& value); // initialization by a value of the data array
  array<T>& resize(unsigned int _size); // resize and put 0 for data value
  array<T>& resize(unsigned int _size, const T& value); // resize and put 'value' for data value

  T& operator[](unsigned int i); // overload operator [], access to an element
  array<T>& operator=(const array<T>& object); // overload operator = for the assignement, the one created by the compiler doesn't work

  // FRIEND
  template <class U>
  friend array<U> operator+(const array<U>& object1, const array<U>& object2); // sum of two class-object
  template <class U>
  friend array<U> operator-(const array<U>& object1, const array<U>& object2); // diff of two class-object
  template <class U>
  friend U operator|(const array<U>& object1, const array<U>& object2); // dot product of two class-object
  template <class U, class I>
  friend array<U> operator*(const I& alpha, const array<U>& object2); //  multiplication by a scalar - wise operation

  // DISPLAY
  template <class U>
  friend ostream& operator<<(ostream& flux, const array<U> &object); // overload the display stream (and file output stream)
};

#endif // ARRAY_CLASS
