/* C-only benchmark for the kdtree library */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "kdtest.h"
#include "kdtree.h"

#define SIZE 400000

int main(int argc, char *argv[]) {

    int num_threads = 4;

    int i;

    double * data = (double *) malloc(3*SIZE*sizeof(double));
    int *f = (int *)malloc(SIZE*sizeof(int));

    double x,y,z;

    srand(2);

    /*
    for(i=0; i<n; i++) {
        data[3*i] = (float) (((double) rand())/RAND_MAX);
        data[3*i+1] = (float) (((double) rand())/RAND_MAX);
        data[3*i+2] = (float) (((double) rand())/RAND_MAX);
    }
    */

    /*
    for(i=0; i<SIZE/2; i++) {
        x =  (((double) rand())/RAND_MAX);
        y =  (((double) rand())/RAND_MAX);
        z =  (((double) rand())/RAND_MAX);
        data[3*i] = x;
        data[3*(i+SIZE/2)] = x;
        data[3*i+1] = y;
        data[3*(i+SIZE/2)+1] = y;
        data[3*i+2] = z;
        data[3*(i+SIZE/2)+2] = z;
    }
    */

    for(i=0; i<SIZE; i++) {
        x =  (((double) rand())/RAND_MAX);
        y =  (((double) rand())/RAND_MAX);
        z =  (((double) rand())/RAND_MAX);
        data[3*i] = x;
        data[3*i+1] = y;
        data[3*i+2] = z;
    }

    int num_fields = 4;

    for(i=0; i<SIZE; i++) {
        f[i] = rand() % num_fields + 1;
    }

    printf("Generated random data...\n");
    kdtree_t data_tree = tree_construct(data, f, SIZE, num_fields);
    printf("Constructed k-d tree...\n");

    verify_tree(data_tree);
    printf("Done verifying\n");
    printf("The tree has %d nodes.\n", count_tree(data_tree));

    //exit(1);

    double radius = 0.05;
    struct timespec start, finish;

    clock_gettime(CLOCK_MONOTONIC, &start);
    //long long * output = pair_count_jackknife(data_tree, data, radius, num_threads);
    long long * output = pair_count_jackknife(data_tree, data, f, SIZE, num_fields, radius, num_threads);
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
