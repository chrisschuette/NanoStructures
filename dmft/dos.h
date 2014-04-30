#ifndef DOS_H
#define DOS_H
namespace dmft {
namespace dos {
double M(double eps);

double N(double eps);

double rho3(double omega);
double rho2(double eps);
double rhoXX(double eps);
double rhoXY(double eps);
}
}
#endif // DOS_H
