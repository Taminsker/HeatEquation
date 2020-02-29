#ifndef CLASS_HEAT_STEADY
#define CLASS_HEAT_STEADY

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>

#include "array.h"
#include "transformFunction.h"

using namespace std;

class heatSteady{
private:
    double deltaX; // step in the X dimension
    double deltaY; // step in the Y dimension
    double dt; // step of time
    unsigned int dimX; // dimension in X
    unsigned int dimY; // dimension in Y
    array<double> x_n; // x_n solution
    array<double> x_np1; // x_{n+1} solution
    array<double> f; // array of the funtion contribution

public:
    heatSteady(double _deltaX, double _deltaY, unsigned int _dimX, unsigned int _dimY, double _dt); // constructor
    ~heatSteady(); // destructpor

    int resolve(unsigned int maxIt, double eps); //resolve the steady problem
    void init(double valueOnTheEdge); // put 'valueOnTheEdge' on the edge
    void saveSolVTK(string nameFile); // export solution in a vtk file
    array<double> productMatrixArray(array<double> matrix); // do the A*matrix product, with matrix an array and A the matrix problem
    double iteration(); // do an iteration
    array<double> getSolution(); // get the solution x_n
};

#endif // CLASS_HEAT_STEADY
