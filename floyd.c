/*
 * Parallel Version of the Floyd algorithm
 * using a Row Striped Block decomposition
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


#define MIN(a,b)           ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))
#define PTR_SIZE           (sizeof(void*))
#define CEILING(i,j)       (((i)+(j)-1)/(j))

int main(int argc, char **argv){

    int rank;
    int processes;

    int **matrix, *matrix_storage;
    int n;
    int size;
    char buffer[255];

    MPI_Status  status;

    FILE* matrix_file;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&processes);

    // Reading matrix size
    if(rank==(processes-1)) {
        matrix_file = fopen("./matrix","r");
        fscanf(matrix_file,"%d",&n);
        //n = strtol(buffer,NULL,10);
    }
    MPI_Bcast(&n,1,MPI_INT,processes-1,MPI_COMM_WORLD);

    //Preparing array
    size = BLOCK_SIZE(rank,processes,n);
    matrix = (int **) malloc(size*sizeof(int));
    if(matrix == NULL){
        MPI_Abort(MPI_COMM_WORLD,-1);
    }

    matrix_storage = (int *) malloc(size*n*sizeof(int));
    if(matrix_storage == NULL){
        MPI_Abort(MPI_COMM_WORLD,-1);
    }

    for(int i=0; i<size; i++){
        matrix[i] = &matrix_storage[i*n];
    }

    //Loading array
    if(rank == processes-1){
        for(int i = 0; i < processes-1; i++){
            int local_size = BLOCK_SIZE(i,processes,n);
          //  fread(matrix_storage,sizeof(int),local_size*n,matrix_file);
            for(int i = 0; i < local_size*n; i++)
                fscanf(matrix_file,"%d",&matrix_storage[i]);
            MPI_Send(matrix_storage,local_size*n,MPI_INT,i,0,MPI_COMM_WORLD);
        }
        for(int i = 0; i < size*n; i++)
            fscanf(matrix_file,"%d",&matrix_storage[i]);
    } else{
        MPI_Recv(matrix_storage,size*n,MPI_INT,processes-1, 0,MPI_COMM_WORLD, &status);
    }

#ifdef _DEBUG
    fprintf(stdout,"Process: %d\n",rank);
    for(int i =0; i <size; i ++){
        for(int k =0; k<n; k++){
            fprintf(stdout,"%d ",matrix[i][k]);
        }
        fprintf(stdout,"\n");
    }
#endif

    // Free matrix and exit
    free(matrix),matrix=NULL;
    free(matrix_storage),matrix_storage=NULL;

    MPI_Finalize();
    return 0;
}