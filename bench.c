/* C-only benchmark for the kdtree library */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "kdtest.h"

#ifndef KDTREE_H
#define KDTREE_H
#include "kdtree.h"
#endif

#define SIZE 400000

int main(int argc, char *argv[]) {

    int num_threads = 4;

    int i;
    int n = SIZE;

    FLOAT *x = (FLOAT *)malloc(n*sizeof(FLOAT));
    FLOAT *y = (FLOAT *)malloc(n*sizeof(FLOAT));
    FLOAT *z = (FLOAT *)malloc(n*sizeof(FLOAT));

    srand(2);

    for(i=0; i<n; i++) {
        x[i] = ((FLOAT) rand())/RAND_MAX;
        y[i] = ((FLOAT) rand())/RAND_MAX;
        z[i] = ((FLOAT) rand())/RAND_MAX;
    }
    
    printf("Generated random data...\n");
    kdtree_t data_tree = tree_construct(n, x, y, z);
    printf("Constructed k-d tree...\n");

    FLOAT radius = 0.05;
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
    long long output = two_point_correlation(data_tree, x, y, z, n, radius, num_threads);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    printf("Sum: %lld\n", output);
    printf("Completed query in %f sec\n", (finish.tv_sec-start.tv_sec)
                                + (finish.tv_nsec-start.tv_nsec)/1e9);
    printf("A node is %d bytes.\n", nodesize());
    verify_main(data_tree,0);
    printf("Done verifying\n");
    printf("The tree has %d nodes.\n", count_main(data_tree));
    return 0;
}
