#ifndef FFTBROADENER_H
#define FFTBROADENER_H

#include "../broadener.h"

namespace nrg {
namespace broadening {
class FFTBroadener : public Broadener
{
public:
    FFTBroadener();
    FFTBroadener(const FFTBroadener& orig);
    FFTBroadener& operator=(const FFTBroadener& orig);
    virtual ~FFTBroadener();

    virtual math::CFunction broaden();
    virtual Broadener* clone() { return new FFTBroadener(*this); }
    virtual void showInfo();
protected:
    int getPadding(double eps);
    double logGaussian(int j);
    double getCrossOver(double frequency);
    double gaussian(double frequency);

    double contribution(double omega, double peakPosition, double peakWeight);
    std::vector<double> m_frequencies;
    double m_sqrtPi;
    double m_epsPadding;
};
}
}
#endif // FFTBROADENER_H
