#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "kdtree.h"

#define SIZE 2000000

int main(int argc, char *argv[]) {

    int prov;
    //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &prov);
    MPI_Init(NULL,NULL);

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

    double *x = (double *)malloc(n*sizeof(double));
    double *y = (double *)malloc(n*sizeof(double));
    double *z = (double *)malloc(n*sizeof(double));

    for(i=0; i<n; i++) {
        x[i] = ((double) rand())/RAND_MAX;
        y[i] = ((double) rand())/RAND_MAX;
        z[i] = ((double) rand())/RAND_MAX;
    }
    
    kdtree_t tree = tree_construct(n, x, y, z);

    double radius = 0.05;
    long long output = two_point_correlation(tree, x, y, z, n, radius,
                                                            MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
