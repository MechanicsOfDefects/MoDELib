/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <MPI/ModelMPIbase.h>
#include <iostream>
#include <mpi.h>
#include <assert.h>
#include <TerminalColors.h>

namespace model
{

    ModelMPIbase::ModelMPIbase(int argc, char *argv[])
    { /*! Constructor
       */
        initMPI(argc, argv);
    }

    ModelMPIbase::ModelMPIbase()
    {
    }

    ModelMPIbase::~ModelMPIbase()
    { /*! Destructor
       */
        MPI_Finalize();
    }

    auto ModelMPIbase::initMPI(int argc, char *argv[]) -> void
    {
        MPI_Initialized(&_mpiInitialized);
        if (!_mpiInitialized)
        {
            MPI_Init(&argc, &argv);

            // Verify
            MPI_Initialized(&_mpiInitialized);
            assert(_mpiInitialized && "MPI not initialized.");

            MPI_Comm_rank(MPI_COMM_WORLD, &_mpiRank);
            MPI_Comm_size(MPI_COMM_WORLD, &_mpiProcs);

            int tempLen;
            char tempName[MPI_MAX_PROCESSOR_NAME];
            MPI_Get_processor_name(tempName, &tempLen);
            _processorName = tempName;

            std::cout << greenBoldColor << "MPI process " << _mpiRank << " of "
                      << _mpiProcs
                      /*     */
                      << " (on node " << _processorName << "), mpiInitialized=" << _mpiInitialized << defaultColor
                      << std::endl;

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    auto ModelMPIbase::mpiRank() -> const int &
    {
        return _mpiRank;
    }

    auto ModelMPIbase::mpiProcs() -> const int &
    {
        return _mpiProcs;
    }

    auto ModelMPIbase::mpiInitialized() -> const int &
    {
        return _mpiInitialized;
    }

    auto ModelMPIbase::processorName() -> const std::string &
    {
        return _processorName;
    }

    // static data
    int ModelMPIbase::_mpiInitialized = 0;
    int ModelMPIbase::_mpiRank = 0;
    int ModelMPIbase::_mpiProcs = 1;
    std::string ModelMPIbase::_processorName = "NO NAME ASSIGNED YET";

} // namespace model
