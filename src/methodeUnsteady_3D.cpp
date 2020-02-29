#include "heatClassUnsteady_3D.h"

heatUnsteady_3D::heatUnsteady_3D( double _deltaX,
                                  double _deltaY,
                                  double _deltaZ,
                                  unsigned int _dimX,
                                  unsigned int _dimY,
                                  unsigned int _dimZ,
                                  double _dt):
                                  deltaX(_deltaX),
        /* constructor */         deltaY(_deltaY),
                                  deltaZ(_deltaZ),
                                  dt(_dt),
                                  dimX(_dimX),
                                  dimY(_dimY),
                                  dimZ(_dimZ),
                                  x_n(_dimX * _dimY * _dimZ),
                                  x_np1(_dimX * _dimY * _dimZ),
                                  f(_dimX * _dimY * _dimZ)
                                  {};


void heatUnsteady_3D::init(double value) // put 'valueOnTheEdge' on the edge, and initialize solution array and call #functionContribution
{
  this->valueOnTheEdge = value;
  this->error = 100;
  this->x_n.init(0.);
  this->x_np1.init(0.);
  this->functionContribution();
};

heatUnsteady_3D::~heatUnsteady_3D(){}; // destructor

void heatUnsteady_3D::functionContribution() // put one in the center
{
  for (unsigned int i = 0; i < this->dimX; i++)
  for (unsigned int j = 0; j < this->dimX; j++)
  for (unsigned int k = 0; k < this->dimZ; k++)
  {
    double distance = sqrt(pow(i * this->deltaX - centerCoalStove_coorX, 2) + pow(j * this->deltaY - centerCoalStove_coorY, 2) + pow(k * this->deltaZ - centerCoalStove_coorZ, 2));

    if (distance < 0.1)
    {
      // this->f[(i * this->dimY + j) * this->dimZ + k] = sin(i * this->deltaX + j * this->deltaY + k * this->deltaZ);
      // this->f[(i * this->dimY + j) * this->dimZ + k] = 1.;
    };
  };
};

void heatUnsteady_3D::saveSolVTK(string nameFile) // save the solution array in a vtk file
{
  ofstream file;
  file.open(nameFile.c_str());

  file << "# vtk DataFile Version 2.0" << '\n';
  file << nameFile << '\n';
  file << "ASCII" << '\n';
  file << "DATASET UNSTRUCTURED_GRID" << '\n';

  file << "POINTS " << (this->dimX * this->dimY * this->dimZ) << " double" << '\n';
  for (unsigned int i = 0; i < this->dimX; ++i)
  for (unsigned int j = 0; j < this->dimY; ++j)
  for (unsigned int k = 0; k < this->dimZ; k++)
  {
    array<double> point(3);
    point[0] = i * this->deltaX;
    point[1] = j * this->deltaY;
    point[2] = k * this->deltaZ;
    file << function_1(point); // transformFunction
  };

  file << "CELLS " << ((this->dimX-1) * (this->dimY-1) * (this->dimZ-1)) << " " << (8+1) * ((this->dimX-1) * (this->dimY-1) * (this->dimZ-1)) << '\n';
  for (unsigned int i = 0; i < this->dimX-1; ++i)
  for (unsigned int j = 0; j < this->dimY-1; ++j)
  for (unsigned int k = 0; k < this->dimZ-1; ++k)
  {
    file << "8 " << (i * this->dimY + j) * this->dimZ + k  << " "  /* node 1 */
                 << ((i + 1) * this->dimY + j) * this->dimZ + k << " "  /* node 2 */
                 << ((i + 1) * this->dimY + (j + 1)) * this->dimZ + k << " "  /* node 3 */
                 << (i * this->dimY + (j + 1)) * this->dimZ + k << " "  /* node 4 */
                 << (i * this->dimY + j) * this->dimZ + (k + 1)  << " "  /* node 5 */
                 << ((i + 1) * this->dimY + j) * this->dimZ + (k + 1) << " "  /* node 6 */
                 << ((i + 1) * this->dimY + (j + 1)) * this->dimZ + (k + 1) << " "  /* node 7 */
                 << (i * this->dimY + (j + 1)) * this->dimZ + (k + 1) << '\n';  /* node 8 */
  };

  file << "CELL_TYPES " << ((this->dimX-1) * (this->dimY-1) * (this->dimZ-1)) << '\n';
  for (unsigned int i = 0; i < this->dimX-1; ++i)
  for (unsigned int j = 0; j < this->dimY-1; ++j)
  for (unsigned int k = 0; k < this->dimZ-1; ++k)
  {
    file << "12" << '\n';
  };

  file << "POINT_DATA " << (this->dimX * this->dimY * this->dimZ) << '\n';
  file << "SCALARS cell double" << '\n';
  file << "LOOKUP_TABLE default" << '\n' << flush;
  for(unsigned int i = 0; i < this->dimX; i++)
  {
    for(unsigned int j = 0; j < this->dimY; j++)
    {
      for(unsigned int k = 0; k < this->dimZ; k++)
      {
        file << this->x_n[(i * this->dimY + j) * this->dimZ + k] << " ";
      };
      file << '\n';
    };
    file << '\n' << '\n';
  };

  file << '\n' << flush;
  file.close();
}


array<double> heatUnsteady_3D::getSolution()
{
  return this->x_n;
}

array<double> heatUnsteady_3D::productMatrixArray(array<double> vec)
{
  array<double> res(this->dimX * this->dimY * this->dimZ, 0.);

  double w_x = this->dt / (this->deltaX * this->deltaX);
  double w_y = this->dt / (this->deltaY * this->deltaY);
  double w_z = this->dt / (this->deltaZ * this->deltaZ);

  double diagonalFactor = 1. + 2. * w_x + 2. * w_y + 2. * w_z;

  unsigned int nodeNumber = 0;
  unsigned sizeTo = this->dimX * this->dimY * this->dimZ;

  for (unsigned int i = 0; i < this->dimX; ++i)
  for (unsigned int j = 0; j < this->dimY; ++j)
  for (unsigned int k = 0; k < this->dimZ; ++k)
  {
    unsigned int itsTheNode = (i * this->dimY + j) * this->dimZ + k;
    res[itsTheNode] = diagonalFactor * vec[itsTheNode];

    nodeNumber = ((i-1) * this->dimY + j) * this->dimZ + k; // node before on x axis : i-1
    if (nodeNumber < sizeTo && i != 0) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_x * vec[nodeNumber];
    };
    nodeNumber = ((i+1) * this->dimY + j) * this->dimZ + k; // node after on x axis : i+1
    if (nodeNumber < sizeTo && i != (this->dimX - 1)) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_x * vec[nodeNumber];
    };

    nodeNumber = (i * this->dimY + (j-1)) * this->dimZ + k; // node before on y axis : j-1
    if (nodeNumber < sizeTo && j != 0) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_y * vec[nodeNumber];
    };
    nodeNumber = (i * this->dimY + (j+1)) * this->dimZ + k; // node after on y axis : j+1
    if (nodeNumber < sizeTo && j != (this->dimY - 1)) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_y * vec[nodeNumber];
    };

    nodeNumber = (i * this->dimY + j) * this->dimZ + (k-1); // node before on z axis : k-1
    if (nodeNumber < sizeTo && k != 0) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_z * vec[nodeNumber];
    };
    nodeNumber = (i * this->dimY + j) * this->dimZ + (k+1); // node after on z axis : k+1
    if (nodeNumber < sizeTo && k != (this->dimZ - 1)) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_z * vec[nodeNumber];
    };
  };

  return res;
}
array<double> heatUnsteady_3D::conjuguateGradient()
{
  // ----- build cond-vector for the edge contribution ------- //

  double w_x = this->dt / (this->deltaX * this->deltaX);
  double w_y = this->dt / (this->deltaY * this->deltaY);
  double w_z = this->dt / (this->deltaZ * this->deltaZ);

  array<double> cond(this->dimX * this->dimY * this->dimZ, 0.);

  /*Impose value on the whole edge */

  // for (unsigned int j = 0; j < this->dimY; ++j)
  // for (unsigned int k = 0; k < this->dimZ; ++k)
  // {
  //   unsigned int i = 0;
  //   cond[(i * this->dimY + j) * this->dimZ + k] += this->valueOnTheEdge * w_x;
  //   i = this->dimX - 1;
  //   cond[(i * this->dimY + j) * this->dimZ + k] += this->valueOnTheEdge * w_x;
  // };
  //
  // for (unsigned int i = 0; i < this->dimX; ++i)
  // for (unsigned int k = 0; k < this->dimZ; ++k)
  // {
  //   unsigned int j = 0;
  //   cond[(i * this->dimY + j) * this->dimZ + k] += this->valueOnTheEdge * w_y;
  //   j = this->dimY - 1;
  //   cond[(i * this->dimY + j) * this->dimZ + k] += this->valueOnTheEdge * w_y;
  // };
  //
  // for (unsigned int i = 0; i < this->dimX; ++i)
  // for (unsigned int j = 0; j < this->dimY; ++j)
  // {
  //   unsigned int k = 0;
  //   cond[(i * this->dimY + j) * this->dimZ + k] += this->valueOnTheEdge * w_z;
  //   k = this->dimZ - 1;
  //   cond[(i * this->dimY + j) * this->dimZ + k] += this->valueOnTheEdge * w_z;
  // };

  /* Impose value on the a part of the edge
     Impose the other part with the last compute value x_{n-1} at (i, j, k) */

  for (unsigned int j = 0; j < this->dimY; ++j)
  for (unsigned int k = 0; k < this->dimZ; ++k)
  {
    unsigned int i = 0;
    cond[(i * this->dimY + j) * this->dimZ + k] += this->x_n[(i * this->dimY + j) * this->dimZ + k] * w_x;
    i = this->dimX - 1;
    cond[(i * this->dimY + j) * this->dimZ + k] += this->x_n[(i * this->dimY + j) * this->dimZ + k] * w_x;
  };

  for (unsigned int i = 0; i < this->dimX; ++i)
  for (unsigned int k = 0; k < this->dimZ; ++k)
  {
    unsigned int j = 0;
    cond[(i * this->dimY + j) * this->dimZ + k] += this->x_n[(i * this->dimY + j) * this->dimZ + k] * w_y;
    j = this->dimY - 1;
    cond[(i * this->dimY + j) * this->dimZ + k] += this->x_n[(i * this->dimY + j) * this->dimZ + k] * w_y;
  };

  for (unsigned int i = 0; i < this->dimX; ++i)
  for (unsigned int j = 0; j < this->dimY; ++j)
  {
    unsigned int k = 0;
    cond[(i * this->dimY + j) * this->dimZ + k] += this->valueOnTheEdge * w_z;
    k = this->dimZ - 1;
    cond[(i * this->dimY + j) * this->dimZ + k] += this->x_n[(i * this->dimY + j) * this->dimZ + k] * w_z;
  };


  /* CONJUGUATE GRADIENT */
  array<double> b = this->x_n + (this->dt * this->f) + cond; // b array in Ax = b
  this->x_np1.init(1.);

  array<double> r = b - this->productMatrixArray(this->x_np1);
  array<double> p = b - this->productMatrixArray(this->x_np1);

  for (unsigned int cmpt = 0; cmpt < this->x_np1.size(); ++cmpt)
  {
    double alpha = (r|r) / (p|this->productMatrixArray(p));
    this->x_np1 = this->x_np1 + (alpha * p);

    array<double> rPlus1 = r - (alpha * this->productMatrixArray(p));

    if (sqrt(rPlus1|rPlus1) < pow(10, -15))
    {
      break;
    };

    double beta = (rPlus1|rPlus1) / (r|r);

    p = rPlus1 + (beta * p);
    r = rPlus1;
  };

  this->error = sqrt((this->x_np1 - this->x_n)|(this->x_np1 - this->x_n)); // ||x_(n+1) - x_n||_2
  return this->x_np1;
}

array<double> heatUnsteady_3D::iteration()
{
  this->x_n = this->conjuguateGradient();
  return this->x_np1;
}

int heatUnsteady_3D::resolve(unsigned int maxIt, double eps)
{
  assert(system("mkdir -p files") != -1);
  unsigned int iterNumber = 0;
  cout << "maximum number of iterations = " << maxIt << endl;
  cout << "epsilon                      = " << eps << endl;

  while ((iterNumber < maxIt) && (this->error > eps))
  {
    this->iteration();
    cout << "iterNumber = " << iterNumber + 1 << endl;
    cout << "error      = " << this->error << endl;

    string nameFile = "files//coalStoveUnsteady_3D_" + to_string(iterNumber) + ".vtk";
    if (iterNumber % 10 == 0)
    {
      saveSolVTK(nameFile);
    }
    iterNumber++;

    this->x_n = this->x_np1;
  }
  return iterNumber;
}
