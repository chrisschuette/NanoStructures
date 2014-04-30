#include "physics.h"

#include <cmath>

namespace math {
namespace physics {
/* Andrew

    function fermi(b,x)
      real(r8) , intent(in)  :: b,x
      real(r8)               :: fermi
      if (x.ge.zero) then
         fermi=exp(-b*x)/(one+exp(-b*x))
      else
         fermi=one/(one+exp(b*x))
      end if
      return
      */
double fermi(double beta, double x) {
    if(x >= 0.0)
        return std::exp(-beta * x) / (1.0 + std::exp(-beta * x));
    else
        return 1.0 / ( 1.0 + std::exp( beta * x ) );
}
double dfermi(double beta, double x) {
    return -beta / ( exp(-beta * x) + 2.0 + exp(beta * x) );
}
}
}
