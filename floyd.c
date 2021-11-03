/*
 * Parallel Version of the Floyd algorithm
 * using a Row Striped Block decomposition
 */

#include <mpi.h>
#include <stdio.h>
#include <malloc.h>


#define MIN(a,b)           ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))
#define PTR_SIZE           (sizeof(void*))
#define CEILING(i,j)       (((i)+(j)-1)/(j))

int main(int argc, char **argv){

    int rank;
    int processes;

    int **matrix, *matrix_storage, i;
    int n;

    FILE* matrix_file;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&processes);


    if(!rank) {
        matrix_file = fopen("../matrix","r");
        fscanf(matrix_file,"%d",n);

        matrix_storage = (int *) malloc(n * sizeof(int));
        matrix = (int **) malloc(n * n * sizeof(int));
        for (i = 0; i < n; i++) matrix[i] = &matrix_storage[i * n];



    }


    MPI_Finalize();
    return 0;
}