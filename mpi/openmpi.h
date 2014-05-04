#ifndef MPI_H
#define MPI_H

/**
 * @file
 *
 * @ingroup mpi
 *
 * @brief Singleton in charge of data exchange via OpenMPI 
 */

/*
 * Copyright (c) 2014 Christoph Schuette.
 *
 * The license and distribution terms for this file may be
 * found in the file LICENSE in this distribution
 */

#include <memory>

namespace math {
    class CFunction;
}

namespace mpi {
    /**
     * @brief The MPI singleton is responsible for data exchange between
     *        compute nodes via OpenMPI.
     * 
     * First and foremost, the MPI class provides convenience functions for the exchange of math::CFunction
     * instances. Secondly, combine allows to reassemble a function which has been
     * jointly calculated on the participating compute nodes. Third, it provides barrier functions
     * to synchronize the nodes.
     */
    class OpenMPI {
    public:
        /**
         * @brief returns the only instance of the OpenMPI class.
         * 
         * If OpenMPI has not been accessed previously an instance is created.
         * The auto_ptr ensures that the object is properly destructed when the 
         * the application closes. OpenMPI is initialized when the object is
         * first constructed.
         * 
         * @return the only instance of the OpenMPI class
         */
        static OpenMPI& getInstance();
        friend class std::auto_ptr<OpenMPI>;

        /**
         * @brief returns the rank of the node.
         * @return rank of the node
         */
        int getRank() const {
            return m_rank;
        }
        /**
         * @brief returns the total number of nodes.
         * @return total number of nodes.
         */
        int getSize() const {
            return m_size;
        }

        //network operations
        /**
         * @brief sends a copy of the math::CFunction from node <i>source</i> to
         * node <i>target</i>.
         * 
         * If the return value of getRank() is neither equal
         * to <i>source</i> nor to <i>target</i> the function
         * does nothing. For getRank()==<i>source</i> a copy of the math::CFunction
         * instance referenced by <i>f</i> is send to node <i>target</i>. The contents
         * of the math::CFunction instance referenced by <i>f</i> on <i>target</i>
         * is replaced with the received data. Function blocks until operation is complete.
         * 
         * @param[in] target target node
         * @param[in] source source node
         * @param[in,out] f data to send (source) / storage location (target)
         */
        void send(int target, int source, math::CFunction& f);
        
        /**
         * @brief sends a copy of double value <i>v</i> from node <i>source</i> to
         * node <i>target</i>.
         * 
         * If the return value of getRank() is neither equal
         * to <i>source</i> nor to <i>target</i> the function
         * does nothing. For getRank()==<i>source</i> the double value <i>v</i>]
         * is send to node <i>target</i>. The contents of the double value referenced
         * by <i>f</i> on <i>target</i> is replaced with the received data.
         * Function blocks until operation is complete.
         * 
         * @param[in] target target node
         * @param[in] source source node
         * @param[in,out] v data to send (source) / storage location (target)
         */
        void send(int target, int source, double &v);
        
        /**
         * @brief assembles the partial information of the individual compute nodes
         *        to a full copy on each one.
         * 
         * Typically (in this program), the calculation of discretized functions 
         * \f$f(x_i)\f$ is parallelized by assigning the individual calculations
         * for different arguments \f$x_i\f$ in a round-robin fashion among the
         * compute nodes. Therefore compute node \f$p\f$ holds the data points
         * \f$ f(x_p), f(x_{p+N}),f(x_{p+2N}), \dots  \f$ with \f$N\f$ the total
         * number of compute nodes. This function reassembles the complete function
         * by collecting all partial information on the master node, rebuilding the 
         * full function and distributing it to all nodes.
         * Function blocks until operation is complete.
         * 
         * @param[in,out] f data contribution (in) / storage for assembled data (out)
         */
        void combine(math::CFunction& f);
        
        /**
         * @brief distributes the math::CFunction instance <i>f</i> on node getRank()==<i>master</i>
         *        to all other nodes.
         * 
         * Function blocks until operation is complete.
         * 
         * @param[in,out] f data to send (master) / storage location (others)
         * @param[in] master node with master copy 
         */
        void sync(math::CFunction& f, int master);
        
        /**
         * @brief distributes the double value <i>v</i> on node getRank()==<i>master</i>
         *        to all other nodes.
         * 
         * Function blocks until operation is complete.
         * 
         * @param[in,out] f data to send (master) / storage location (others)
         * @param[in] master node with master copy 
         */
        void sync(double& v, int master);
        
        /**
         * @brief blocks until all nodes have reached this function call.
         */
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
