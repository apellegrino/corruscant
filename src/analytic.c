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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "kdtree.h"

#define NUM_THREADS 4
#define NTRIALS 20

int main(void) {

    FILE * fp = fopen("analytic.csv", "w");

    long long size;
    for(size=4000L; size<=4096000L; size *= 2) {

        fprintf(fp, "%lld,", size);    

        int trial;
        for(trial=0; trial<NTRIALS; trial++) {
            srand(trial+2);

            double radius = 0.05;
            int i, j;

            long long * output;
            double *data, *newdata;


            //struct timespec start, finish;

            double x,y,z, analytic;

            data = (double *) malloc(3*size*sizeof(double));

            //make subset of querying data
            newdata = (double *) malloc(3*size*sizeof(double));


            for(i=0; i<size; i++) {
                x =  (((double) rand())/RAND_MAX);
                y =  (((double) rand())/RAND_MAX);
                z =  (((double) rand())/RAND_MAX);
                data[3*i] = x;
                data[3*i+1] = y;
                data[3*i+2] = z;
            }

            printf("Generated random data...\n");

            kdtree_t data_tree;
            data_tree = tree_construct(data, NULL, size, 1);
            printf("Constructed k-d tree...\n");

            j = 0;
            for(i=0; i<size; i++) {
                if(data[3*i] > radius && data[3*i] < 1.0 - radius &&
                    data[3*i+1] > radius && data[3*i+1] < 1.0 - radius &&
                    data[3*i+2] > radius && data[3*i+2] < 1.0 - radius) {
                    // inside safe zone

                    newdata[3*j] = data[3*i];
                    newdata[3*j+1] = data[3*i+1];
                    newdata[3*j+2] = data[3*i+2];
                    j++;
                }
            }

            printf("%d points in safe zone\n", j+1);

            //clock_gettime(CLOCK_MONOTONIC, &start);
            output = pair_count(data_tree, newdata, NULL, j+1, 1, radius, NUM_THREADS);
            //clock_gettime(CLOCK_MONOTONIC, &finish);

            analytic = pow(1.0 - 2*radius, 3) *
                       (4.0 / 3.0 * M_PI * pow(size,2) * pow(radius,3) + size);

            printf("\n");
            printf("Total           = %lld\n", *output);
            printf("Analytic result = %f\n", analytic);
            printf("\n");

            fprintf(fp, "%lld,", *output);

            //printf("Completed query in %f sec\n", (finish.tv_sec-start.tv_sec)
            //                            + (finish.tv_nsec-start.tv_nsec)/1e9);

            free(data_tree.node_data);
            free(output);

            free(data);
            free(newdata);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 0;
}
