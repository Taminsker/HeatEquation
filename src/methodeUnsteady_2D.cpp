#include "heatClassUnsteady_2D.h"

heatUnsteady_2D::heatUnsteady_2D( double _deltaX,
                                  double _deltaY,
                                  unsigned int _dimX,
                                  unsigned int _dimY,
                                  double _dt):
                                  deltaX(_deltaX),
        /* constructor */         deltaY(_deltaY),
                                  dt(_dt),
                                  dimX(_dimX),
                                  dimY(_dimY),
                                  x_n(_dimX * _dimY),
                                  x_np1(_dimX * _dimY),
                                  f(_dimX * _dimY)
                                  {};


void heatUnsteady_2D::init(double value) // put 'valueOnTheEdge' on the edge, and initialize solution array and call #functionContribution
{
  this->valueOnTheEdge = value;
  this->error = 100;
  this->x_n.init(0.);
  this->x_np1.init(0.);
  this->functionContribution();
};

heatUnsteady_2D::~heatUnsteady_2D(){}; // destructor

void heatUnsteady_2D::functionContribution() // put one in the center
{
  for(unsigned int i = 0; i < this->dimX; i++)
  for(unsigned int j = 0; j < this->dimX; j++)
  {
    double distance = sqrt(pow(i * this->deltaX - centerCoalStove_coorX, 2) + pow(j * this->deltaY - centerCoalStove_coorY, 2));

    if (distance < 0.1)
    {
      // this->f[i * this->dimY + j] = sin(i * this->deltaX + j * this->deltaY);
      // this->f[i * this->dimY + j] = 1.;
    };
  };
};

void heatUnsteady_2D::saveSolVTK(string nameFile) // save the solution array in a vtk file
{
  ofstream file;
  file.open(nameFile.c_str());

  file << "# vtk DataFile Version 2.0" << '\n';
  file << nameFile << '\n';
  file << "ASCII" << '\n';
  file << "DATASET UNSTRUCTURED_GRID" << '\n';

  file << "POINTS " << (this->dimX * this->dimY) << " double" << '\n';
  for (unsigned int i = 0; i < this->dimX; ++i)
  for (unsigned int j = 0; j < this->dimY; ++j)
  {
    array<double> point(3);
    point[0] = i * this->deltaX;
    point[1] = j * this->deltaY;
    point[2] = 0.;
    file << function_1(point); // transformFunction
  };

  file << "CELLS " << ((this->dimX-1) * (this->dimY-1)) << " " << (4+1) * ((this->dimX-1) * (this->dimY-1)) << '\n';
  for (unsigned int i = 0; i < this->dimX-1; ++i)
  for (unsigned int j = 0; j < this->dimY-1; ++j)
  {
    file << "4 " << i*(this->dimY)+j << " "  /* node 1 */
                 << (i+1)*(this->dimY)+j << " "  /* node 2 */
                 << (i+1)*(this->dimY)+j+1 << " "  /* node 3 */
                 << i*(this->dimY)+j+1 << '\n';  /* node 4 */
  };

  file << "CELL_TYPES " << ((this->dimX-1) * (this->dimY-1)) << '\n';
  for (unsigned int i = 0; i < this->dimX-1; ++i)
  for (unsigned int j = 0; j < this->dimY-1; ++j)
  {
    file << "9" << '\n';
  };

  file << "POINT_DATA " << (this->dimX * this->dimY) << '\n';
  file << "SCALARS cell double" << '\n';
  file << "LOOKUP_TABLE default" << '\n' << flush;
  for(unsigned int i = 0; i < this->dimX; i++)
  {
    for(unsigned int j = 0; j < this->dimY; j++)
    {
      file << this->x_n[i * this->dimY + j] << " ";
    };
    file << '\n';
  };

  file << '\n' << flush;
  file.close();
}


array<double> heatUnsteady_2D::getSolution()
{
  return this->x_n;
}

array<double> heatUnsteady_2D::productMatrixArray(array<double> vec)
{
  array<double> res(this->dimX * this->dimY, 0.);

  double w_x = this->dt / (this->deltaX * this->deltaX);
  double w_y = this->dt / (this->deltaY * this->deltaY);

  double diagonalFactor = 1. + 2. * w_x + 2. * w_y;

  unsigned int nodeNumber = 0;
  unsigned int sizeTo = this->dimX * this->dimY;

  for (unsigned int i = 0; i < this->dimX; ++i)
  for (unsigned int j = 0; j < this->dimY; ++j)
  {
    unsigned int itsTheNode = i * this->dimY + j;

    res[itsTheNode] = diagonalFactor * vec[itsTheNode];

    nodeNumber = (i-1) * this->dimY + j; // node before on x axis : i-1
    if (nodeNumber < sizeTo && i != 0) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_x * vec[nodeNumber];
    };
    nodeNumber = (i+1) * this->dimY + j; // node after on x axis : i+1
    if (nodeNumber < sizeTo && i != (this->dimX - 1)) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_x * vec[nodeNumber];
    };

    nodeNumber = i * this->dimY + (j-1); // node before on y axis : j-1
    if (nodeNumber < sizeTo && j != 0) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_y * vec[nodeNumber];
    };
    nodeNumber = i * this->dimY + (j+1); // node after on y axis : j+1
    if (nodeNumber < sizeTo && j != (this->dimY - 1)) // not just next to the edge :: carefull
    {
      res[itsTheNode] -= w_y * vec[nodeNumber];
    };
  };

  return res;
}
array<double> heatUnsteady_2D::conjuguateGradient()
{
  // ----- build cond-vector for the edge contribution ------- //

  double w_x = this->dt / (this->deltaX * this->deltaX);
  double w_y = this->dt / (this->deltaY * this->deltaY);

  array<double> cond(this->dimX * this->dimY, 0.);

  /*Impose value on the whole edge */

  for (unsigned int j = 0; j < this->dimY; ++j)
  {
    unsigned int i = 0;
    cond[i * this->dimY + j] += this->valueOnTheEdge * w_x;
    i = this->dimX - 1;
    cond[i * this->dimY + j] += this->valueOnTheEdge * w_x;
  };

  for (unsigned int i = 0; i < this->dimX; ++i)
  {
    unsigned int j = 0;
    cond[i * this->dimY + j] += this->valueOnTheEdge * w_y;
    j = this->dimY - 1;
    cond[i * this->dimY + j] += this->valueOnTheEdge * w_y;
  };

  /* Impose value on the a part of the edge
     Impose the other part with the last compute value x_{n-1} at (i, j) */

  // for (unsigned int j = 0; j < this->dimY; ++j)
  // {
  //   unsigned int i = 0;
  //   cond[i * this->dimY + j] += this->x_n[i * this->dimY + j] * w_x;
  //   i = this->dimX - 1;
  //   cond[i * this->dimY + j] += this->valueOnTheEdge * w_x;
  // };
  //
  // for (unsigned int i = 0; i < this->dimX; ++i)
  // {
  //   unsigned int j = 0;
  //   cond[i * this->dimY + j] += this->x_n[i * this->dimY + j] * w_y;
  //   j = this->dimY - 1;
  //   cond[i * this->dimY + j] += this->x_n[i * this->dimY + j] * w_y;
  // };

  /* CONJUGUATE GRADIENT */
  array<double> b = this->x_n + (this->dt * this->f) + cond;

  this->x_np1.init(1.);

  array<double> r = b - this->productMatrixArray(this->x_np1);
  array<double> p = b - this->productMatrixArray(this->x_np1);

  unsigned int k = 0;

  for (k = 0; k < this->x_np1.size(); ++k)
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
  this->error = sqrt((this->x_np1 - this->x_n)|(this->x_np1 - this->x_n));
  return this->x_np1;
}

array<double> heatUnsteady_2D::iteration()
{
  this->x_n = this->conjuguateGradient();
  return this->x_np1;
}

int heatUnsteady_2D::resolve(unsigned int maxIt, double eps)
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

    string nameFile = "files//coalStoveUnsteady_2D_" + to_string(iterNumber) + ".vtk";
    if (iterNumber % 5 == 0)
    {
      saveSolVTK(nameFile);
    }
    iterNumber++;

    this->x_n = this->x_np1;
  }
  return iterNumber;
}
