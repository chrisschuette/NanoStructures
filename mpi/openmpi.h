#ifndef MPI_H
#define MPI_H

#include <memory>

namespace math {
    class CFunction;
}

namespace mpi {
class OpenMPI
{
public:
    static OpenMPI& getInstance();
    friend class std::auto_ptr<OpenMPI>;

    int getRank() const { return m_rank; }
    int getSize() const { return m_size; }

    //network operations
    void send(int target, int source, math::CFunction& f);
    void send(int target, int source, double &v);
    void combine(math::CFunction& f);
    void sync(math::CFunction& f, int master);
    void sync(double& v, int master);
    void sync();


protected:
    OpenMPI();
    virtual ~OpenMPI();
    static std::auto_ptr< OpenMPI > m_ptr;
    int m_rank;
    int m_size;
};
#define HEAD if(!mpi::OpenMPI::getInstance().getRank())
}
#endif // MPI_H
