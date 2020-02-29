#include "heatClassSteady.h"

heatSteady::heatSteady(double _deltaX,
  double _deltaY,
  unsigned int _dimX,
  unsigned int _dimY,
  double _dt):
  deltaX(_deltaX),
  deltaY(_deltaY),
  dt(_dt),
  dimX(_dimX),
  dimY(_dimY),
  x_n(_dimX * _dimY),
  x_np1(_dimX * _dimY),
  f(_dimX * _dimY)
  {};

  void heatSteady::init(double valueOnTheEdge)
  {
    for(unsigned int i = 0; i < this->dimX; ++i)
    {
      this->x_n[i * this->dimY + 0]              = valueOnTheEdge;
      this->x_n[i * this->dimY + this->dimY - 1] = valueOnTheEdge;
    }
    for(unsigned int j = 0; j < this->dimY; ++j)
    {
      this->x_n[0 * this->dimY + j]                = valueOnTheEdge;
      this->x_n[(this->dimX - 1) * this->dimY + j] = valueOnTheEdge;
    }
  }

  heatSteady::~heatSteady(){};

  void heatSteady::saveSolVTK(string nomFichier)
  {
    ofstream file;
    file.open(nomFichier.c_str());

    file << "# vtk DataFile Version 3.0" << '\n';
    file << "cell" << '\n';
    file << "ASCII" << '\n';
    file << "DATASET STRUCTURED_POINTS" << '\n';
    file << "DIMENSIONS " << this->dimX << " " << this->dimY << " " << 1 << '\n';
    file << "ORIGIN " << 0 << " " << 0 << " " << 0 << '\n';
    file << "SPACING " << 1.0 << " " << 1.0 << " " << 1.0 << '\n';
    file << "POINT_DATA " << this->dimX * this->dimY << '\n';
    file << "SCALARS cell double" << '\n';
    file << "LOOKUP_TABLE default" << '\n' << flush;

    for(unsigned int i = 0; i < this->dimX; i++)
    {
      for(unsigned int j = 0; j < this->dimY; j++)
      {
        file << this->x_n[i * this->dimY + j] << " ";
      }
      file << '\n';
    }
    file << '\n' << flush;
    file.close();
  }

  double heatSteady::iteration()
  {
    double error = 0.;
    double w_x = this->dt / (this->deltaX * this->deltaX);
    double w_y = this->dt / (this->deltaY * this->deltaY);
    double diagonalFactor = 1. - 2. * w_x - 2. * w_y;

    for (unsigned int i = 1; i < this->dimX - 1; ++i)
    for (unsigned int j = 1; j < this->dimY - 1; ++j)
    {
      this->x_np1[i * this->dimY + j] = diagonalFactor * this->x_n[i * this->dimY + j]
      + w_x * (this->x_n[(i - 1) * this->dimY + j] + this->x_n[(i + 1) * this->dimY + j])
      + w_y * (this->x_n[i * this->dimY + j - 1]   + this->x_n[i * this->dimY + j + 1]);

      error += pow(this->x_np1[i * this->dimY + j] - this->x_n[i * this->dimY + j], 2);
    }
    return sqrt(error);
  }

  int heatSteady::resolve(unsigned int maxIt, double eps)
  {
    assert(system("mkdir -p files") != -1);
    unsigned int iterNumber = 0;

    double error = 100;

    cout << "maximum number of iterations = " << maxIt << endl;
    cout << "epsilon                      = " << eps  << endl;

    this->x_np1 = this->x_n;


    while ((iterNumber < maxIt) && (error > eps))
    {
      error = iteration();
      string nameFile = "files//heatSteady_" + to_string(iterNumber) + ".vtk";
      if (iterNumber % 5 == 0)
      {
        saveSolVTK(nameFile);
      }
      cout << " iteration number = " << iterNumber << " error = " << error << endl;
      iterNumber++;

      this->x_n = this->x_np1;
    }
    return iterNumber;
  }

  array<double> heatSteady::getSolution()
  {
    return this->x_n;
  }

  array<double> heatSteady::productMatrixArray(array<double> matrix)
  {

    double w_x = this->dt / (this->deltaX * this->deltaX);
    double w_y = this->dt / (this->deltaY * this->deltaY);
    double diagonalFactor = 1. + 2. * w_x + 2. * w_y;

    array<double> res(this->dimX * this->dimY, 0.);

    for (unsigned int i = 1; i < this->dimX - 1; ++i)
    for (unsigned int j = 1; j < this->dimY - 1; ++j)
    {
      res[i * this->dimY + j] = diagonalFactor * matrix[i * this->dimY + j]
      - w_x * (matrix[(i-1) * this->dimY + j]+ matrix[(i+1) * this->dimY + j])
      - w_y * (matrix[i * this->dimY + j - 1] + matrix[i * this->dimY + j + 1]);
    }

    for(unsigned int i = 0; i < this->dimX; ++i)
    {
      res[i * this->dimY + 0] = 1.;
      res[i * this->dimY + this->dimY -1 ] =  1.;
    }

    for(unsigned int j = 0; j < this->dimY; ++j)
    {
      res[0 * this->dimY + j] =  1.;
      res[(this->dimX - 1) * this->dimY + j] =  1.;
    }
    return res;
  }
