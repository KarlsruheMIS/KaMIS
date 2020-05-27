/******************************************************************************
 * mpi_tools.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#include <mpi.h>
#include <unistd.h>

#include "mpi_tools.h"

mpi_tools::mpi_tools() {
                
}

mpi_tools::~mpi_tools() {
                

}

void mpi_tools::non_active_wait_for_root() {
        int rank   = MPI::COMM_WORLD.Get_rank();
        int size   = MPI::COMM_WORLD.Get_size();
        int MASTER = 0;

        if(rank == MASTER) {
                //wake up call
                bool wakeup = true;
                for( int to = 1; to < size; to++) {
                        MPI::COMM_WORLD.Send(&wakeup, 1, MPI::BOOL, to, 0);
                }
        } else {
                //non-busy waiting:
                bool stop = false;
                do {
                        usleep(5000);
                        stop = MPI::COMM_WORLD.Iprobe(MASTER,0);
                } while(!stop);

                bool wakeup = true;
                MPI::COMM_WORLD.Recv(&wakeup, 1, MPI::BOOL, MASTER, 0);
        }
}

