/* C-only benchmark for the kdtree library */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "kdtest.h"

#ifndef KDTREE_H
    #include "kdtree.h"
#endif

#define SIZE 400000

int main(int argc, char *argv[]) {

    int num_threads = 16;

    int i;
    int n = SIZE;

    datum_t * data = (datum_t *)malloc(n*sizeof(datum_t));
    int *f = (int *)malloc(n*sizeof(int));

    srand(2);

    for(i=0; i<n; i++) {
        data[i].value[0] = ((double) rand())/RAND_MAX;
        data[i].value[1] = ((double) rand())/RAND_MAX;
        data[i].value[2] = ((double) rand())/RAND_MAX;
    }

    int num_fields = 4;

    for(i=0; i<n; i++) {
        f[i] = 1 + (long long) rand() * num_fields / ((long long)RAND_MAX+1);
    }

    printf("Generated random data...\n");
    kdtree_t data_tree = tree_construct((double *) data, f, n, num_fields);
    printf("Constructed k-d tree...\n");

    verify_tree(data_tree);
    printf("Done verifying\n");
    printf("The tree has %d nodes.\n", count_tree(data_tree));

    double radius = 0.05;
    struct timespec start, finish;

    clock_gettime(CLOCK_MONOTONIC, &start);
    //long long * output = pair_count_jackknife(data_tree, data, radius, num_threads);
    long long * output = pair_count_jackknife(data_tree, (double *) data, f, SIZE, num_fields, radius, num_threads);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    for(i=0; i<num_fields+1; i++) {
        printf("%lld ", output[i]);
    }
    printf("\n");
    //printf("Sum: %lld\n", output);
    printf("Completed query in %f sec\n", (finish.tv_sec-start.tv_sec)
                                + (finish.tv_nsec-start.tv_nsec)/1e9);
    printf("A node is %d bytes.\n", nodesize());
    return 0;
}
