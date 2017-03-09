/* C-only benchmark for the kdtree library */

#include <stdlib.h>
#include <stdio.h>
#include "kdtree.h"
#include "mpi.h"

#define SIZE 4000000

int main(int argc, char *argv[]) {

    printf("argv %d\n", argv[1]);
    int prov;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &prov);

    printf("argv %d\n", argv[1]);
    int num_threads = 12;
    //MPI_Init(NULL,NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(!rank) {
        if(prov == MPI_THREAD_SINGLE)
            printf("Providing SINGLE\n");
        else if (prov == MPI_THREAD_FUNNELED)
            printf("Providing FUNNELED\n");
        else if (prov == MPI_THREAD_SERIALIZED)
            printf("Providing SERIALIZED\n");
        else if (prov == MPI_THREAD_MULTIPLE)
            printf("Providing MULTIPLE\n");
    }

    int i;
    int n = SIZE;

    FLOAT *x = (FLOAT *)malloc(n*sizeof(FLOAT));
    FLOAT *y = (FLOAT *)malloc(n*sizeof(FLOAT));
    FLOAT *z = (FLOAT *)malloc(n*sizeof(FLOAT));

    for(i=0; i<n; i++) {
        x[i] = ((FLOAT) rand())/RAND_MAX;
        y[i] = ((FLOAT) rand())/RAND_MAX;
        z[i] = ((FLOAT) rand())/RAND_MAX;
    }
    
    printf("Generated random data...\n");
    kdtree_t tree = tree_construct(n, x, y, z);
    printf("Constructed k-d tree...\n");

    FLOAT radius = 0.05;
    long long output = two_point_correlation(tree, x, y, z, n, radius, num_threads,
                                                            MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
