#ifndef CLASS_HEAT_UNSTEADY_2D
#define CLASS_HEAT_UNSTEADY_2D

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

extern double centerCoalStove_coorX;
extern double centerCoalStove_coorY;

class heatUnsteady_2D{
private:
    double deltaX; // step in the X dimension
    double deltaY; // step in the Y dimension
    double dt; // step of time
    double valueOnTheEdge; //  value on the edge
    double error;
    unsigned int dimX; // dimension in X
    unsigned int dimY; // dimension in Y
    array<double> x_n; // x_n solution
    array<double> x_np1; // x_{n+1} solution
    array<double> f; // array of the funtion contribution

public:
    heatUnsteady_2D(double _deltaX, double _deltaY, unsigned int _dimX, unsigned int _dimY, double _dt); // constructor
    ~heatUnsteady_2D(); //  destructor
    void functionContribution(); // compute the function contribution

    int resolve(unsigned int maxIt, double eps); // resolve the Unsteady problem
    void init(double valueOnTheEdge); // put 'valueOnTheEdge' on the edge, and initialize solution array and call #functionContribution
    void saveSolVTK(string nameFile); // save the solution array in a vtk file
    array<double> productMatrixArray(array<double> matrix); // do the
    array<double> iteration(); // do an iteration
    array<double> getSolution(); // get the solution x_n
    array<double> conjuguateGradient(); // conjuguate gradient function

    friend array<double> function_1(array<double> point); // transform function number 1 :: identity
};

#endif // CLASS_HEAT_UNSTEADY_2D
