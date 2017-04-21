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

    int num_threads = 32;

    int i;
    int n = SIZE;

    double *x = (double *)malloc(n*sizeof(double));
    double *y = (double *)malloc(n*sizeof(double));
    double *z = (double *)malloc(n*sizeof(double));

    srand(2);

    for(i=0; i<n; i++) {
        x[i] = ((double) rand())/RAND_MAX;
        y[i] = ((double) rand())/RAND_MAX;
        z[i] = ((double) rand())/RAND_MAX;
    }

    array3d_t data;
    data.x = x;
    data.y = y;
    data.z = z;
    data.size = n;
    
    printf("Generated random data...\n");
    kdtree_t data_tree = tree_construct(data, NULL);
    printf("Constructed k-d tree...\n");

    double radius = 0.05;
    struct timespec start, finish;

    clock_gettime(CLOCK_MONOTONIC, &start);
    long long output = pair_count_jackknife(data_tree, data, radius, -1, num_threads);
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
