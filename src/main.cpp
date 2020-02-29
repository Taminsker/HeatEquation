#include "heatClassSteady.h"
#include "heatClassUnsteady_2D.h"
#include "heatClassUnsteady_3D.h"

#include "transformFunction.h"

double centerCoalStove_coorX, centerCoalStove_coorY, centerCoalStove_coorZ;

int main(int argc, char * argv[])
{
  if (argc < 4)
  {
    cerr << "usage  : " << argv[0] << " Nx Ny maxIt prec" << endl <<  endl;
    cout << "Nx     : number of nodes in the X dimension," << endl;
    cout << "Ny     : number of nodes in the Y dimension," << endl;
    cout << "maxIt  : number of iteration maximal," << endl;
    cout << "prec   : epsilon to stop." << endl;
    return -1;
  }
  if (argc == 5)
  {
    unsigned int Nx = atoi(argv[1]);
    unsigned int Ny = atoi(argv[2]);

    // array<double> a(3);
    // a[0] = Nx;
    // a[1] = Ny;
    // a = function_1(a);
    // cout << "Nx = " << Nx << " Ny = " << Ny << endl;
    //
    // Nx = ((a[0] != a[0]) ? Nx : abs(a[0]));
    // Ny = ((a[1] != a[1]) ? Ny : abs(a[1]));
    // cout << "Nx = " << Nx << " Ny = " << Ny << endl;



    unsigned int maxIt = atoi(argv[3]);
    double eps = atof(argv[4]);

    double hx = 1. / Nx;
    double hy = 1. / Ny;
    double dt = min (hx * hx / 16., hy * hy / 16.);

    cout << "------- Steady problem -----------------" << endl;
    heatSteady coalStoveSteady(hx, hy, Nx + 2, Ny + 2, dt);

    coalStoveSteady.init(1.);
    coalStoveSteady.resolve(maxIt, eps);

    cout << "------- Unsteady problem 2D---------------" << endl;
    heatUnsteady_2D coalStoveUnsteady_2D(hx, hy, Nx, Ny, dt);
    centerCoalStove_coorX = (1-hx) / 2;
    centerCoalStove_coorY = (1-hy) / 2;

    coalStoveUnsteady_2D.init(1.);
    coalStoveUnsteady_2D.resolve(maxIt, eps);
  } else {
    unsigned int Nx = atoi(argv[1]);
    unsigned int Ny = atoi(argv[2]);
    unsigned int Nz = atoi(argv[3]);

    // array<double> a(3);
    // a[0] = Nx;
    // a[1] = Ny;
    // a[2] = Nz;
    // a = function_1(a);
    // cout << "Nx = " << Nx << " Ny = " << Ny << " Nz = " << Nz << endl;
    //
    // Nx = ((a[0] != a[0]) ? Nx : abs(a[0]));
    // Ny = ((a[1] != a[1]) ? Ny : abs(a[1]));
    // Nx = ((a[2] != a[2]) ? Nx : abs(a[2]));
    // cout << "Nx = " << Nx << " Ny = " << Ny << " Nz = " << Nz << endl;

    unsigned int maxIt = atoi(argv[4]);
    double eps = atof(argv[5]);

    double hx = 1. / Nx;
    double hy = 1. / Ny;
    double hz = 1. / Nz;

    double dt = min(hx * hx / 4., min(hy * hy / 4., hz * hz / 4.));

      cout << "------- Unsteady problem 3D---------------" << endl;
    heatUnsteady_3D coalStoveUnsteady_3D(hx, hy, hz, Nx, Ny, Nz, dt);

    centerCoalStove_coorX = (1-hx) / 2;
    centerCoalStove_coorY = (1-hy) / 2;
    centerCoalStove_coorZ = (1-hz) / 2;


    coalStoveUnsteady_3D.init(1.);
    coalStoveUnsteady_3D.resolve(maxIt, eps);
  }


  return 0;
};
