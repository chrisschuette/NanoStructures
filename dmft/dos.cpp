#include "dos.h"
#include <cmath>

namespace dmft {
namespace dos {
double M(double eps) {
    return 0.0746756967363403 + (8675.0 * eps * eps) / 6283008;
}

double N(double eps) {
    return -0.10744121287128713 - (51.0 * eps) / 25856.0 - (81.0 * eps * eps) / 51712.0;
}

double rho3(double omega) {
    double dos = M(omega) * sqrt(36 - omega * omega);
    if (omega < -2.0)
        dos += N(-omega) * sqrt(4.0 - (omega + 4.0) * (omega + 4.0));
    if (omega > 2.0)
        dos += N(omega) * sqrt(4.0 - (omega - 4.0) * (omega - 4.0));

    return dos / M_PI;
}

double rho2(double eps) {
    return 0.14282921950711333 + (-0.051299342112440104 + 0.0003538996863908572 * eps * eps) * log(fabs(eps));
}

double rhoXX(double eps) {
    if(std::abs(eps) > 4.0)
        return 0.0;
    if(std::abs(eps) < 1e-20)
        return 0.4053206358913109 - 0.04161446460559167*pow(eps,2);
    return 0.4053206358913109 - 0.04161446460559167*pow(eps,2) + 2.0 * 0.00625*pow(eps,2)*log(std::abs(eps));
}

double rhoXY(double eps) {
    if(std::abs(eps) > 4.0)
        return 0.0;
    if(std::abs(eps) < 1e-20)
        return 0.0;

    return -0.36051079779815587*eps + 0.010232094641393531*pow(eps,3) + 0.00005767847611739453*pow(eps,5) - 6.854604420339911e-7*pow(eps,7) + 0.2*eps*log(std::abs(eps)) - 2.0 * 0.0020833333333333333*pow(eps,3)*log(std::abs(eps));
}

}
}
