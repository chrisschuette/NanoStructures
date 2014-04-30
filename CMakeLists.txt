cmake_minimum_required (VERSION 2.6)
project (dmft)

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-reorder -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wno-reorder -O3 -DNDEBUG")

SITE_NAME(mySite)
MESSAGE( ${MPI_INCLUDE_PATH} )
MESSAGE("mySite='${mySite}'")

IF("${mySite}" STREQUAL "iChris.local")
   set(LIBS ${LIBS} config++ gsl mpfr gmp fftw3 /opt/local/lib/libboost_program_options-mt.a )
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework vecLib")
   MESSAGE( "MacBook Pro" )
   include_directories(/opt/local/include )
   link_directories(/opt/local/lib )
ENDIF("${mySite}" STREQUAL "iChris.local")

IF("${mySite}" STREQUAL "cheops0")
   set(LIBS ${LIBS} config++ mkl_intel_lp64 mkl_sequential mkl_core boost_program_options fftw3 mpfr gmp gsl )
   MESSAGE( "Cheops Cluster" )
   include_directories(/home/cschuet0/include )
   link_directories(/home/cschuet0/lib )
ENDIF("${mySite}" STREQUAL "cheops0")

add_executable (${PROJECT_NAME} 
./cmd/cmd.h
./config/configuration.h
./debug/timer.h
./dmft/cfcoefficients.h
./dmft/dmft.h
./dmft/dodmft.h
./dmft/dos.h
./err/exception.h
./math/cfunction.h
./math/exceptions.h
./math/fourier.h
./math/function.h
./math/hilberttransformer.h
./math/integral.h
./math/lapack.h
./math/matrix.h
./math/mkl_types.h
./math/pfunction.h
./math/physics.h
./math/vector.h
./mpi/openmpi.h
./nano/donano.h
./nano/layerparameters.h
./nano/nanostructure.h
./nrg/broadener.h
./nrg/broadening/fftbroadener.h
./nrg/chain/flatband.h
./nrg/chain/hybridizationprovider.h
./nrg/chainprovider.h
./nrg/donrg.h
./nrg/hilbertspacetable.h
./nrg/hilbertsubspace.h
./nrg/nrg.h
./utils/complex.h
./utils/mpreal.h
./utils/utils.h
./utils/utils_io.h
./cmd/cmd.cpp
./config/configuration.cpp
./dmft/dmft.cpp
./dmft/dodmft.cpp
./dmft/dos.cpp
./err/exception.cpp
./main.cpp
./math/cfunction.cpp
./math/fourier.cpp
./math/function.cpp
./math/hilberttransformer.cpp
./math/integral.cpp
./math/matrix.cpp
./math/pfunction.cpp
./math/physics.cpp
./math/vector.cpp
./mpi/openmpi.cpp
./nano/donano.cpp
./nano/layerparameters.cpp
./nano/nanostructure.cpp
./nrg/broadener.cpp
./nrg/broadening/fftbroadener.cpp
./nrg/chain/flatband.cpp
./nrg/chain/hybridizationprovider.cpp
./nrg/chainprovider.cpp
./nrg/donrg.cpp
./nrg/hilbertspacetable.cpp
./nrg/hilbertsubspace.cpp
./nrg/nrg.cpp
./utils/mpreal.cpp
./utils/utils.cpp
./utils/utils_io.cpp)

target_link_libraries (${PROJECT_NAME} ${LIBS} ${MPI_LIBRARIES})
