#ifndef HILBERTSUBSPACE_H
#define HILBERTSUBSPACE_H

#include "../math/matrix.h"
#include "../math/vector.h"

namespace nrg {
class HilbertSubSpace
{
public:
    HilbertSubSpace(int Q, int Sz);

    //getters and setters
    inline int getQ() { return m_Q; }
    inline int getSz() { return m_Sz; }
    inline int getKeptStates() { return m_nKeptStates; }
    void setKeptStates(int keptStates) { m_nKeptStates = keptStates; }
    inline math::Matrix& getHamiltonian() { return m_H; }
    void setHamiltonian(math::Matrix& H) { m_H = H; }
    inline math::Vector& getEnergies() { return m_E; }
    void setEnergies(math::Vector& E) { m_E = E; }
    void setR(int i, int r) { m_r[i] = r; }
    int getR(int i) { return m_r[i]; }
    inline math::Matrix& getChainElementsUp() { return m_chainOperatorElementsUp; }
    inline math::Matrix& getChainElementsDown() { return m_chainOperatorElementsDown; }
    inline math::Matrix& getDensityMatrixEigenBasis() { return m_densityMatrixEigenBasis; }
    inline math::Matrix& getDensityMatrix() { return m_densityMatrix; }
    inline math::Matrix& getLocalMatrixElementUp() { return m_fUp; }
    inline math::Matrix& getLocalMatrixElementDown() { return m_fDown; }
    inline math::Matrix& getLocalMatrixElementUp2() { return m_fffUp; }
    inline math::Matrix& getLocalMatrixElementDown2() { return m_fffDown; }
protected:
    int m_Q;
    int m_Sz;
    int m_r[4];
    int m_nKeptStates;

    math::Matrix m_H;
    math::Vector m_E;
    math::Matrix m_chainOperatorElementsUp;
    math::Matrix m_chainOperatorElementsDown;
    math::Matrix m_densityMatrix;
    math::Matrix m_densityMatrixEigenBasis;
    math::Matrix m_fUp;
    math::Matrix m_fDown;
    math::Matrix m_fffUp;
    math::Matrix m_fffDown;
};
}

#endif // HILBERTSUBSPACE_H
