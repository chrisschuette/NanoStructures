#include "openmpi.h"
#include <mpi.h>
#include <iostream>
#include "../math/cfunction.h"

namespace mpi {

std::auto_ptr<OpenMPI> OpenMPI::m_ptr = std::auto_ptr<OpenMPI>(0);

OpenMPI::OpenMPI()
{
    MPI::Init();
    m_rank = MPI::COMM_WORLD.Get_rank();
    m_size = MPI::COMM_WORLD.Get_size();
}

OpenMPI& OpenMPI::getInstance() {
    if(m_ptr.get() == 0) {
        m_ptr = std::auto_ptr<OpenMPI>(new OpenMPI());
    }
    return *m_ptr.get();
}

void OpenMPI::send(int target, int source, math::CFunction& f) {
    if(target == source)
        return;
    if(m_rank == source) {
        MPI::COMM_WORLD.Send(&f.m_length, 1, MPI_INT, target, 0);
        MPI::COMM_WORLD.Send(f.m_arguments, f.m_length, MPI_DOUBLE, target, 1);
        MPI::COMM_WORLD.Send(f.m_values, f.m_length, MPI_DOUBLE_COMPLEX, target, 2);
    }
    if(m_rank == target ) {
        int length;
        MPI::Status status;
        MPI::COMM_WORLD.Recv( &length, 1, MPI::INT, source, 0,
        status );
        f.resize(length);
        MPI::COMM_WORLD.Recv( f.m_arguments, length, MPI_DOUBLE, source, 1,
        status );
        MPI::COMM_WORLD.Recv( f.m_values, length, MPI_DOUBLE_COMPLEX, source, 2,
        status );
    }
}

void OpenMPI::sync() {
    MPI::COMM_WORLD.Barrier();
}

void OpenMPI::send(int target, int source, double& v) {
    if(target == source)
        return;
    if(m_rank == source)
        MPI::COMM_WORLD.Send(&v, 1, MPI_DOUBLE, target, 0);
    if(m_rank == target ) {
        MPI::Status status;
        MPI::COMM_WORLD.Recv( &v, 1, MPI::DOUBLE, source, 0,
        status );
    }
}

void OpenMPI::combine(math::CFunction& f) {
    if(m_rank == 0) {
        //gather
        int size = 0;
        math::CFunction* parts = new math::CFunction[m_size];
        parts[0] = f;
        size += f.getSize();
        for(int r = 1; r < m_size; r++) {
            // get the parts
            send(0, r, parts[r]);
            size += parts[r].getSize();
        }
        //allocate
        f.resize(size);
        //combine
        int counter = 0;
        for(int i = 0; i < parts[0].getSize(); i++) {
            for(int r = 0; r < m_size; r++) {
                if(i < parts[r].getSize()) {
                    f.set(counter, parts[r].getArgument(i), parts[r].getValue(i));
                    counter++;
                }
            }
        }
	delete [] parts;
        // distribute
        for(int r = 1; r < m_size; r++)
            send(r, 0, f);
    } else {
        // send our part to head node
        send(0, m_rank, f);
        // receive completed function
        send(m_rank, 0, f);
    }
}

void OpenMPI::sync(double &v, int master) {
    if(m_rank == master) {
        for(int r = 0; r < m_size; r++)
            if(r != master)
                send(r, master, v);
    } else {
        send(m_rank, master, v);
    }
}

void OpenMPI::sync(math::CFunction& f, int master) {
    if(m_rank == master) {
        for(int r = 0; r < m_size; r++)
            if(r != master)
                send(r, master, f);
    } else {
        send(m_rank, master, f);
    }
}

OpenMPI::~OpenMPI() {
    MPI::Finalize();
}

}
