#ifndef LAYERPARAMETERS_H
#define LAYERPARAMETERS_H

#include <string>

namespace nano {
class LayerParameters
{
public:
    LayerParameters();
    double U;
    double mu;
    double rhoB;
    double initialV;
    std::string selfEnergyFile;
    std::string occupancyFile;
};
}
#endif // LAYERPARAMETERS_H
