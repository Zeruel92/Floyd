/*
 * Parallel Version of the Floyd algorithm
 * using a Row Striped Block decomposition
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#define MIN(a, b)           ((a)<(b)?(a):(b))
#define BLOCK_LOW(id, p, n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j, p, n) (((p)*((j)+1)-1)/(n))

#define MAX_ITERATIONS 1

int main(int argc, char **argv) {

    int rank;
    int processes;

    int **matrix, *matrix_storage, *k_row;
    int n;
    int size;
    int low_value;

    double elapsed_time;

    MPI_Status status;
    MPI_Request request;

    FILE *matrix_file;

#ifdef _DEBUG
    int debug = 1;
#endif

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    // Reading matrix size
    if (rank == (processes - 1)) {
        matrix_file = fopen("./matrix", "r");
        fscanf(matrix_file, "%d", &n);
        if (processes > n) {
            fprintf(stderr, "The number of processors exceed matrix dimensions\n");
            MPI_Abort(MPI_COMM_WORLD, -2);
        }
    }
    MPI_Bcast(&n,1,MPI_INT,processes-1,MPI_COMM_WORLD);

    //Preparing array
    size = BLOCK_SIZE(rank, processes, n);
    low_value = BLOCK_LOW(rank, processes, n);
    matrix = (int **) malloc(size * sizeof(int));
    if(matrix == NULL){
        fprintf(stderr, "Failed allocating matrix\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    matrix_storage = (int *) malloc(size*n*sizeof(int));
    if(matrix_storage == NULL){
        fprintf(stderr, "Failed allocating matrix_storage\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    for(int i=0; i<size; i++){
        matrix[i] = &matrix_storage[i*n];
    }

    //Loading array
    if(rank == processes-1){
        for(int i = 0; i < processes-1; i++) {
            int local_size = BLOCK_SIZE(i, processes, n);
            for (int j = 0; j < local_size * n; j++)
                fscanf(matrix_file, "%d", &matrix_storage[j]);
            if (processes > 1) MPI_Irsend(matrix_storage, local_size * n, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
        }
        for (int i = 0; i < size * n; i++)
            fscanf(matrix_file, "%d", &matrix_storage[i]);
    } else {
        MPI_Recv(matrix_storage, size * n, MPI_INT, processes - 1, 0, MPI_COMM_WORLD, &status);
    }

#ifdef _DEBUG

    if (!rank) {
        fprintf(stdout, "INPUT MATRIX:\n");
        fflush(stdout);
    } else {
        MPI_Recv(&debug, 1, MPI_INT, (rank - 1) % processes, 1, MPI_COMM_WORLD, &status);
    }
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < n; k++) {
            fprintf(stdout, "%d ", matrix[i][k]);
            fflush(stdout);
        }
        fprintf(stdout, "\n");
        fflush(stdout);
    }
    if (processes > 1) MPI_Send(&debug, 1, MPI_INT, (rank + 1) % processes, 1, MPI_COMM_WORLD);
#endif

    k_row = (int *) malloc(sizeof(int) * n);
    if (k_row == NULL) {
        fprintf(stderr, "Failed allocating k_row\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //Starting computation
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    for (int iterations = 0; iterations < MAX_ITERATIONS; iterations++) {

        for (int k = 0; k < n; k++) {
            int owner = BLOCK_OWNER(k, processes, n);
            if (rank == owner) {
#ifdef _DEBUG
                fprintf(stdout, "K Value %d process %d k local %d size %d\n", k, rank, (k - low_value), size);
#endif
                for (int j = 0; j < n; j++)
                    k_row[j] = matrix[k - low_value][j];
            }
            MPI_Bcast(k_row, n, MPI_INT, owner, MPI_COMM_WORLD);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < n; j++) {
                    matrix[i][j] = MIN(matrix[i][j], matrix[i][k] + k_row[j]);
                }
            }
#ifdef _DEBUG
            sleep(1);
            if (!rank) {
                fprintf(stdout, "K ROW\n");
                fflush(stdout);
                for (int j = 0; j < n; j++) {
                    fprintf(stdout, "%d ", k_row[j]);
                    fflush(stdout);
                }
                fprintf(stdout, "\nIntermediate MATRIX:\n");
                fflush(stdout);
            } else {
                MPI_Recv(&debug, 1, MPI_INT, (rank - 1) % processes, 1, MPI_COMM_WORLD, &status);
            }
            for (int i = 0; i < size; i++) {
                for (int k = 0; k < n; k++) {
                    fprintf(stdout, "%d ", matrix[i][k]);
                    fflush(stdout);
                }
                fprintf(stdout, "\n");
                fflush(stdout);
            }
            if (processes > 1) MPI_Send(&debug, 1, MPI_INT, (rank + 1) % processes, 1, MPI_COMM_WORLD);
#endif
        }
    }

#ifdef _DEBUG
    sleep(1);
    if (!rank) {
        fprintf(stdout, "OUTPUT MATRIX:\n");
        fflush(stdout);
    } else {
        MPI_Recv(&debug, 1, MPI_INT, (rank - 1) % processes, 1, MPI_COMM_WORLD, &status);
    }
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < n; k++) {
            fprintf(stdout, "%d ", matrix[i][k]);
            fflush(stdout);
        }
        fprintf(stdout, "\n");
        fflush(stdout);
    }
    if (processes > 1) MPI_Send(&debug, 1, MPI_INT, (rank + 1) % processes, 1, MPI_COMM_WORLD);
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    elapsed_time /= MAX_ITERATIONS;

    if (!rank) fprintf(stdout, "Processors: %d Iterations: %d, Time: %f", processes, MAX_ITERATIONS, elapsed_time);

    // Free matrix and exit
    free(matrix), matrix = NULL;
    free(matrix_storage), matrix_storage = NULL;
    free(k_row), k_row = NULL;

    MPI_Finalize();
    return 0;
}