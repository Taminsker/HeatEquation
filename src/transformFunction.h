#ifndef TRANSFORM_FUNCTION
#define TRANSFORM_FUNCTION

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string>

#include "array.h"

extern double centerCoalStove_coorX;
extern double centerCoalStove_coorY;
extern double centerCoalStove_coorZ;

array<double> function_1(array<double> point);
// array<double> function_2(array<double> point);

array<double> toCircle_XY(array<double> point, double radius);
array<double> toCircle_XZ(array<double> point, double radius);
array<double> toCircle_YZ(array<double> point, double radius);
array<double> toSphere(array<double> point, double radius);
array<double> rotation(array<double> point, double theta, double gamma, double lambda);
array<double> translation(array<double> point, double xDir, double yDir, double zDir);

#endif // TRANSFORM_FUNCTION
