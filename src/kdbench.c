/* 
 * Copyright (C) 2016-2017 Andrew Pellegrino
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* C-only benchmark for the kdtree library */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "kdtest.h"
#include "kdtree.h"

#define SIZE 400000
#define NUM_THREADS 4
#define NUM_FIELDS 4

int main(void) {

    int i;

    double * data = (double *) malloc(3*SIZE*sizeof(double));
    int *f = (int *)malloc(SIZE*sizeof(int));

    double x,y,z;

    srand(2);

    /* Create duplicates to check for handling */

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

    for(i=0; i<SIZE; i++) {
        f[i] = rand() % NUM_FIELDS;
    }
    printf("Generated random data...\n");

    kdtree_t data_tree = tree_construct(data, f, SIZE, NUM_FIELDS);
    printf("Constructed k-d tree...\n");

    verify_tree(data_tree);
    printf("Done verifying\n");
    printf("The tree has %d nodes.\n", count_tree(data_tree));

    double radius = 0.05;
    struct timespec start, finish;

    clock_gettime(CLOCK_MONOTONIC, &start);
    long long * output = pair_count(data_tree, data, f, SIZE, NUM_FIELDS, radius, NUM_THREADS);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    long long sum = 0;
    for(i=0; i<NUM_FIELDS*NUM_FIELDS; i++) {
        printf("%lld ", output[i]);
        sum += output[i];
    }
    printf("Total = %lld\n", sum);
    printf("\n");

    printf("Completed query in %f sec\n", (finish.tv_sec-start.tv_sec)
                                + (finish.tv_nsec-start.tv_nsec)/1e9);

    free(output);

    clock_gettime(CLOCK_MONOTONIC, &start);
    output = pair_count(data_tree, data, NULL, SIZE, 1, radius, NUM_THREADS);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    sum = 0;
    for(i=0; i<NUM_FIELDS; i++) {
        printf("%lld ", output[i]);
        sum += output[i];
    }
    printf("Total = %lld\n", sum);
    printf("\n");

    printf("Completed query in %f sec\n", (finish.tv_sec-start.tv_sec)
                                + (finish.tv_nsec-start.tv_nsec)/1e9);
    printf("A node is %lu bytes.\n", sizeof(node_t));

    free(data);
    free(f);

    free(data_tree.node_data);
    free(output);

    return 0;
}
