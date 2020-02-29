#include "transformFunction.h"

// array<double> function_1(array<double> point) // identity
// {
//   point[0] = point[0];
//   point[1] = point[1];
//   point[2] = point[2];
//   return point;
// }

array<double> function_1(array<double> point) // (x,y,z) --> (x, y, 6z)
{
  point[0] = point[0];
  point[1] = point[1];
  point[2] = 6*point[2];
  return point;
}

// array<double> function_1(array<double> point) // essai rate de la casserole
// {
//   if (point[0] <= 0.2 || point[1] <= 0.2 || point[0] > 0.7 || point[1] > 0.7)
//   {
//     point[0] = point[0];
//     point[1] = point[1];
//     point[2] = point[2];
//   } else
//   {
//     point[0] = point[0];
//     point[1] = point[1];
//     point[2] = 0.25 * point[2];
//   }
//
//   return point;
// }

// array<double> function_1(array<double> point) // mapOnCircle(x, y) radius = 1 && z = 6z
// {
//   point = toCircle_XY(point, 1.);
//   // point = translation(point, point[2], point[2], 0);
//   // point = rotation(point, 0.5 * point[2] * M_PI, 0.5 * point[2] * M_PI, 0);
//   point[2] *= 3;
//   return point;
// }

// --------------- toCircle ---------------//
array<double> toCircle_XY(array<double> point, double radius)
{
  double x = 2 * (point[0] - centerCoalStove_coorX);
  double y = 2 * (point[1] - centerCoalStove_coorY);

  double xprime = x * sqrt(1 - y * y / 2);
  double yprime = y * sqrt(1 - x * x / 2);

  point[0] = radius * (xprime / 2 + centerCoalStove_coorX);
  point[1] = radius * (yprime / 2 + centerCoalStove_coorY);
  point[2] = point[2];

  return point;
};
array<double> toCircle_XZ(array<double> point, double radius)
{
  double x = 2 * (point[0] - centerCoalStove_coorX);
  double z = 2 * (point[2] - centerCoalStove_coorZ);

  double xprime = x * sqrt(1 - z * z / 2);
  double zprime = z * sqrt(1 - x * x / 2);

  point[0] = radius * (xprime / 2 + centerCoalStove_coorX);
  point[1] = point[1];
  point[2] = radius * (zprime / 2 + centerCoalStove_coorZ);

  return point;
};
array<double> toCircle_YZ(array<double> point, double radius)
{
  double y = 2 * (point[1] - centerCoalStove_coorY);
  double z = 2 * (point[2] - centerCoalStove_coorZ);

  double yprime = y * sqrt(1 - z * z / 2);
  double zprime = z * sqrt(1 - y * y / 2);

  point[0] = point[0];
  point[1] = radius * (yprime / 2 + centerCoalStove_coorY);
  point[2] = radius * (zprime / 2 + centerCoalStove_coorZ);

  return point;
};

// ---------- sphere --------------------//
array<double> toSphere(array<double> point, double radius)
{
  double x = 2 * (point[0] - centerCoalStove_coorX);
  double y = 2 * (point[1] - centerCoalStove_coorY);
  double z = 2 * (point[2] - centerCoalStove_coorZ);

  double xprime = x * sqrt(1 - y * y / 2 - z * z / 2);
  double yprime = y * sqrt(1 - x * x / 2 - z * z / 2);
  double zprime = z * sqrt(1 - x * x / 2 - y * y / 2);

  point[0] = radius * (xprime / 2 + centerCoalStove_coorX);
  point[1] = radius * (yprime / 2 + centerCoalStove_coorY);
  point[2] = radius * (zprime / 2 + centerCoalStove_coorZ);

  return point;
};


// --------------- rotation ---------------//
array<double> rotation(array<double> point, double theta, double gamma, double lambda)
{
  double xR0 = point[0];
  double yR0 = point[1];
  double zR0 = point[2];

  // rotation autour de x //
  double xR1 = xR0;
  double yR1 = yR0 * cos(theta) - zR0 * sin(theta);
  double zR1 = yR0 * sin(theta) + zR0 * cos(theta);

  // rotation autour de y //
  double xR2 = xR1 * cos(gamma) + zR1 * sin(gamma);
  double yR2 = yR1;
  double zR2 = zR1 * cos(gamma) - xR1 * sin(gamma);

  // rotation autour de z //
  double xR3 = xR2 * cos(lambda) - zR2 * sin(lambda);
  double yR3 = xR2 * sin(lambda) + zR2 * cos(lambda);
  double zR3 = yR2;

  point[0] = xR3;
  point[1] = yR3;
  point[2] = zR3;

  return point;
};

// --------------- translation ---------------//
array<double> translation(array<double> point, double xDir, double yDir, double zDir)
{
  double x = point[0] + xDir;
  double y = point[1] + yDir;
  double z = point[2] + zDir;

  point[0] = x;
  point[1] = y;
  point[2] = z;

  return point;
};
