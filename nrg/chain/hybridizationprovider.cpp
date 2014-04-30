#include "hybridizationprovider.h"
#include "../../math/cfunction.h"
#include "../../math/exceptions.h"

#include <fstream>

namespace nrg {
namespace chain {
HybridizationProvider::HybridizationProvider() : ChainProvider()
{
}

void HybridizationProvider::buildChain() {
    m_up = buildChain(m_hybridizationUp);
    m_down = m_up;
}

mpfr::mpreal HybridizationProvider::sign(mpfr::mpreal a, mpfr::mpreal b)
{
    mpfr::mpreal pref;
    if(b >= 0.0)
        return a;
    else
        return -a;
}

tWilsonChain HybridizationProvider::buildChain(math::PFunction& h) {
    mpfr::mpreal::set_default_prec(2048); // super super high precision

    bool negativeValues = false;

    //check that we don't have negative hybridization values
    for(int i = 0; i < h.getSize(); i++) {
        if(h.getValue(i) < 0) {
            std::cout << h.getArgument(i) << " " << h.getValue(i) << std::endl;
            negativeValues = true;
        }
    }

    if(negativeValues) {
        std::ofstream hProblem("h_prob.dat");
        for(int n = 0; n < h.getSize(); n++)
            hProblem << h.getArgument(n) << " " << h.getValue(n) << std::endl;
        hProblem.close();
        std::cerr << "WARNING: negative values in hybridization function" << std::endl;
        std::cerr << "integrals probably ill convergent!" << std::endl;
        //        throw HybridizationException("Hybridization function has neg. values.");
    }

    // avoid problems
    for(int n = 0; n < h.getSize(); n++)
        if(h.getValue(n) < 1e-4)
            h.setValue(n, 1e-4);

    //calculate correction factor
    double A = (m_lambda + 1.0)/(2.0 * m_lambda - 2.0) * std::log(m_lambda);

    int nmax = (m_length+1);

    mpfr::mpreal* x_minus = new mpfr::mpreal[2*nmax];
    mpfr::mpreal* x_plus = new mpfr::mpreal[2*nmax];

    for(int n = 0; n < 2*nmax; n++)
    {
        x_plus[n] = mpfr::pow(m_lambda,-n);
        x_minus[n] = -x_plus[n];
        if(mpfr::_isnan(x_plus[n])) {
            std::cerr << "x_plus is nan." << std::endl;
            throw error::math::Nan("x_plus is nan.");
        }
        if(mpfr::_isnan(x_minus[n])) {
            std::cerr << "x_minus is nan." << std::endl;
            throw error::math::Nan("x_minus is nan.");
        }
    }

    mpfr::mpreal xi_0 = h.integrate(-1,1);
//    std::cout << "xi_0: " << xi_0 << std::endl;
    mpfr::mpreal* gamma_plus = new mpfr::mpreal[2*nmax];
    mpfr::mpreal* gamma_minus = new mpfr::mpreal[2*nmax];

    mpfr::mpreal* xi_plus = new mpfr::mpreal[2*nmax];
    mpfr::mpreal* xi_minus = new mpfr::mpreal[2*nmax];

    mpfr::mpreal* rho_plus = new mpfr::mpreal[2*nmax];
    mpfr::mpreal* rho_minus = new mpfr::mpreal[2*nmax];

    for(int n = 0; n < ( 2*nmax - 1 ); n++)
    {
        rho_plus[n] = h.integrate(x_plus[n+1],x_plus[n]);
        gamma_plus[n] = mpfr::sqrt(rho_plus[n]);
        if(mpfr::_isnan(gamma_plus[n])) {
            std::cerr << "gamma_plus is nan." << std::endl;
           // throw error::math::Nan("gamma_plus is nan.");
        }
        xi_plus[n] = h.integrate_arg(x_plus[n+1],x_plus[n]) / rho_plus[n];
        if(mpfr::_isnan(xi_plus[n])) {
            std::cerr << "xi_plus is nan." << std::endl;
           // throw error::math::Nan("xi_plus is nan.");
        }
        rho_minus[n] = h.integrate(x_minus[n],x_minus[n+1]);
        gamma_minus[n] = mpfr::sqrt(rho_minus[n]);
        if(mpfr::_isnan(gamma_minus[n])) {
            std::cerr << "gamma_minus is nan." << std::endl;
         //   throw error::math::Nan("gamma_minus is nan.");
        }
        xi_minus[n] = h.integrate_arg(x_minus[n],x_minus[n+1]) / rho_minus[n];
        if(mpfr::_isnan(xi_minus[n])) {
            std::cerr << "xi_minus is nan." << std::endl;
          //  throw error::math::Nan("xi_minus is nan.");
        }
        if(rho_plus[n] == 0)
        {
            std::cerr << "rho_plus[n] == 0" << std::endl;
            gamma_plus[n] = 0.0;
            xi_plus[n] = 0.0;
        }
        if(rho_minus[n] == 0)
        {
            std::cerr << "rho_minus[n] == 0" << std::endl;
            gamma_minus[n] = 0.0;
            xi_minus[n] = 0.0;
        }
    }

    mpfr::mpreal* um1 = new mpfr::mpreal[2*nmax];
    mpfr::mpreal* vm1 = new mpfr::mpreal[2*nmax];

    mpfr::mpreal* um2 = new mpfr::mpreal[2*nmax];
    mpfr::mpreal* vm2 = new mpfr::mpreal[2*nmax];

    mpfr::mpreal* normu = new mpfr::mpreal[2*nmax];
    mpfr::mpreal* normv = new mpfr::mpreal[2*nmax];

    mpfr::mpreal* u = new mpfr::mpreal[2*nmax];
    mpfr::mpreal* v = new mpfr::mpreal[2*nmax];

    for(int m = 0; m < 2*nmax; m++)
    {
        normu[m] = 0.0;
        normv[m] = 0.0;
    }

    mpfr::mpreal epsn, eps_loc, sum,epsnm1,epsn_prev;
    mpfr::mpreal epssum = 0.0;

    tWilsonChain wilsonChain;


    for(int n = -1; n < nmax; n++) //do n=-1,Nmax
    {

        if(n==-1)
        {
            epsn = mpfr::sqrt(xi_0/M_PI);
            if(mpfr::_isnan(epsn))
                throw error::math::Nan("epsn nan.");

            eps_loc = 0.0;
        }


        if(n==0) {
            epsn = 0.0;
            eps_loc = 0.0;
            sum = 0.0;
            for(int m = 0; m < (2*nmax -1); m++) { //do m=1,2*Nmax
                um2[m] = gamma_plus[m]/mpfr::sqrt(xi_0);
                vm2[m] = gamma_minus[m]/mpfr::sqrt(xi_0);
                if(mpfr::_isnan(um2[m]))
                    throw error::math::Nan("um2[m] nan.");
                if(mpfr::_isnan(vm2[m]))
                    throw error::math::Nan("vm2[m] nan.");
                normu[m] = normu[m] + mpfr::pow(um2[m],2);
                normv[m] = normv[m] + mpfr::pow(vm2[m],2);
                sum = sum + mpfr::pow(um2[m],2) + mpfr::pow(vm2[m],2);
                epsn= epsn + mpfr::pow((xi_plus[m] * um2[m]),2) + mpfr::pow((xi_minus[m] * vm2[m]),2);
                eps_loc= eps_loc + xi_plus[m] * mpfr::pow((um2[m]),2) + xi_minus[m] * mpfr::pow((vm2[m]),2);
            }
            epsn = mpfr::sqrt(epsn-mpfr::pow((eps_loc),2));
            if(mpfr::_isnan(epsn))
                throw error::math::Nan("epsn is nan.");
            epssum = epssum + mpfr::pow((epsn),2);
            sum = 0.0;
            for(int m = 0; m < 2*nmax-1;m++) { //do m=1,2*Nmax
                um1[m] = um2[m] * (xi_plus[m] - eps_loc)/ epsn;
                vm1[m] = vm2[m] * (xi_minus[m] - eps_loc) / epsn;
                normu[m] = normu[m] + mpfr::pow((um1[m]),2);
                normv[m] = normv[m] + mpfr::pow((vm1[m]),2);
                sum = sum + mpfr::pow((um1[m]),2) + mpfr::pow((vm1[m]),2);
            }
            epsnm1 = epsn;
        }

        if(n>0) {
            sum = 0.0;
            eps_loc = 0.0;
            for(int m = 0; m < (2*nmax -1); m++) { //do m=1,2*Nmax
                sum = sum + mpfr::pow((xi_plus[m] * um1[m]),2) + mpfr::pow((xi_minus[m] * vm1[m]),2);
                eps_loc = eps_loc + xi_plus[m] * mpfr::pow((um1[m]),2) + xi_minus[m] * mpfr::pow(( vm1[m]),2);
            }
            if(mpfr::_isnan(eps_loc)) {
                eps_loc = 0.0;
            }
            epsn = mpfr::sqrt(sum - mpfr::pow((epsnm1),2) - mpfr::pow((eps_loc),2)); // wird negativ
            if(mpfr::_isnan(epsn))
            {
                std::cout << "epsn nan." << std::endl;
                epsn = 1.0;
            }
            epssum = epssum + mpfr::pow((epsn),2);
            sum = 0.0;
            for(int m = 0; m < (2*nmax - 1); m++) { //do m=1,2*Nmax
                u[m] = um1[m] * (xi_plus[m] - eps_loc)/ epsn - epsnm1* um2[m]/epsn;
                v[m] = vm1[m] * (xi_minus[m] - eps_loc)/ epsn - epsnm1* vm2[m]/epsn;
                sum = sum + mpfr::pow((u[m]),2) + mpfr::pow((v[m]),2);
            }
            // Korrektur der Summenregel
            for(int m = 0; m < (2*nmax - 1); m++) { //do m=1,2*Nmax
                u[m] = u[m]/mpfr::sqrt(sum);
                v[m] = v[m]/mpfr::sqrt(sum);
                if (normu[m] > 1.0) normu[m] = 1.0;
                if (normv[m] > 1.0) normv[m] = 1.0;
                if (mpfr::pow(u[m],2) > (1.0-normu[m])) u[m] = sign(mpfr::sqrt(1.0-normu[m]),u[m]);
                if (mpfr::pow(v[m],2) > (1.0-normv[m])) v[m] = sign(mpfr::sqrt(1.0-normv[m]),v[m]);
            }

            sum = 0.0;
            for(int m = 0; m < (2*nmax - 1); m++) { //do m=1,2*Nmax
                sum = sum + mpfr::pow((u[m]),2) + mpfr::pow((v[m]),2);
            }
            for(int m = 0; m < (2*nmax - 1); m++) { //do m=1,2*Nmax
                u[m] = u[m]/mpfr::sqrt(sum);
                v[m] = v[m]/mpfr::sqrt(sum);
                normu[m] = normu[m] + mpfr::pow((u[m]),2);
                normv[m] = normv[m] + mpfr::pow((v[m]),2);
            }

            epsnm1 = epsn;
            for(int m = 0; m < (2*nmax - 1); m++) { //do m=1,2*Nmax
                um2[m] = um1[m];
                um1[m] = u[m];
                vm2[m] = vm1[m];
                vm1[m] = v[m];
            }
        }



        if (n != -1) {
            double oolambda = 1.0/m_lambda;
            double pref = std::sqrt(oolambda) * (2.0 / (1.0 + oolambda));
            //double pref = 1.0;

            if(n == 0)
                epsn_prev *= std::sqrt(A);

            if(mpfr::_isnan(epsn_prev))
                throw error::math::Nan("epsn_prev is nan.");
            if(mpfr::_isnan(eps_loc))
                throw error::math::Nan("eps_loc is nan.");

            wilsonChain.push_back(tOrbital(
                                      pref*epsn_prev*std::pow(m_lambda,((double) n)/2.0),
                                      pref*eps_loc*std::pow(m_lambda,((double) n)/2.0)
                                      ));
        }
        epsn_prev = epsn;

    }
    if(m_symmetry_PH == true) {
        for(unsigned int i = 0; i < wilsonChain.size(); i++)
            wilsonChain.at(i).second = 0.0;
    }


    //free
    delete [] um1;
    delete [] vm1;

    delete [] um2;
    delete [] vm2;

    delete [] normu;
    delete [] normv;

    delete [] u;
    delete [] v;

    delete []  gamma_plus;
    delete []  gamma_minus;

    delete []  xi_plus;
    delete []  xi_minus;

    delete []  rho_plus;
    delete []  rho_minus;


    return wilsonChain;
}


void HybridizationProvider::setHybridization(const math::CFunction& hUp, const math::CFunction& hDown) {
    m_hybridizationUp = hUp;
    m_hybridizationDown = hDown;
}

void HybridizationProvider::setHybridization(const math::CFunction& hUp) {
    m_hybridizationUp = hUp;
    m_hybridizationDown = math::PFunction();
}

HybridizationProvider::~HybridizationProvider() {

}


}
}
