/*
 * Parallel Version of the Floyd algorithm
 * using a Row Striped Block decomposition
 */

#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv){

    int rank;
    int processes;

    int **matrix, *matrix_storage, i;
    int m,n;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&processes);



    MPI_Finalize();
    return 0;
}