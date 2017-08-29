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
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include "kdtree.h"

static double * _data_query;
static int *_field_query;
static node_t * _tree_data;
static double rad;

static inline void count_node(double n, double rsq, field_counter_t * counter,
                       int this_field, int query_field, int this_weight)
{
    if (n < rsq) {
        counter->array[query_field * counter->num_fields + this_field] += this_weight;
    }
    return;
}

static inline double norm2(datum_t * a, datum_t * b)
{
    double sum = 0.0;
    int i;
    for(i=0; i<NDIM; i++) {
        sum += (a->value[i] - b->value[i])*(a->value[i] - b->value[i]);
    }
    return sum;
}

static inline int left_child(int p)
{
    return 2*p;
}

static inline int right_child(int p)
{
    return 2*p+1;
}

/*
 * Query how many points in the tree with head p lie within radius r of point q
 */
static void radius(int pi, int qi, field_counter_t * counter)
{
    node_t p = *(_tree_data+pi);
    int this_field = p.data.field;

    int query_field = 0;
    if(_field_query != NULL) {
        query_field = _field_query[qi];
    }

    double val;

    double rsq;
    double pos_upper = 0.0, pos_lower = 0.0, point = 0.0;

    rsq = rad*rad;

    datum_t datum;
    int i;
    for(i=0;i<NDIM;i++) {
        datum.value[i] = _data_query[NDIM*qi+i];
    }
    double n = norm2(&p.data, &datum);

    if( !p.has_lchild ) { /* no children */
        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
        return;

    } else if ( !p.has_rchild ) { /* one child */
        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
        radius(left_child(pi),qi,counter);
        return;

    } else { /* two children */
        val = _data_query[NDIM*qi+p.dim];
        pos_upper = val + rad;
        pos_lower = val - rad;
        point = p.data.value[p.dim];

        if (pos_upper < point) {
            radius(left_child(pi),qi,counter);
        } else if (pos_lower > point) {
            radius(right_child(pi),qi,counter);
        } else {
            count_node(n, rsq, counter, this_field, query_field, p.data.weight);
            radius(left_child(pi),qi,counter);
            radius(right_child(pi),qi,counter);
        }
        return;
    }
    
}

typedef struct thread_args {
    int rank;
    int start;
    int stop; // actually one after the last index to stop at
    int * finished;
    field_counter_t * counter;
} thread_args_t;

static inline int min(int a, int b)
{
    return a < b ? a : b;
}

/*
void assign_idx(int rank, int size, int n_array, int * start, int * stop)
{
    int n_array_local = n_array / size + ( n_array % size > rank );

    *start = rank * ( n_array / size ) + min( rank, n_array % size );
    *stop = *start + n_array_local;
}
*/

field_counter_t * init_field_counter(int qsize, int tsize)
{
    field_counter_t * c = malloc(sizeof(field_counter_t));
    c->array = calloc(qsize*tsize, sizeof(long long));
    c->num_fields = qsize;
    return c;
}

void free_field_counter(field_counter_t * c)
{
    free(c->array);
    free(c);
}

void * twopoint_wrap(void *voidargs)
{

    int i;
    thread_args_t *targs;
    targs = (thread_args_t *) voidargs;

    //int rank = targs->rank;
    int start = targs->start;
    int stop = targs->stop;

    //printf("%d %d %d\n", targs->rank, start, stop);
   
    for(i=start; i<stop; i++) {
        radius(1, i, targs->counter);
    }

    *(targs->finished) = 1;

    return NULL;
}

/* 
 * Compute the sum of the number of data points in the tree which lie within
 * radius r of each of the points in the array
 */
long long * pair_count(kdtree_t tree, double * data,
                                int * fields, int length, int num_fields,
                                double r, int num_threads)
{
    int i, j;
    int next_assign = 0;
    _tree_data = tree.node_data;
    _data_query = data;
    _field_query = fields;
    rad = r;

    if(num_threads > length / 20) {
        printf("More than %d threads!\n", length / 20);
        exit(EXIT_FAILURE);
    }

    int step = length / num_threads / 20;
    
    pthread_t threads[num_threads];
    thread_args_t targs[num_threads];
    int * started = calloc(num_threads, sizeof(int));
    int * finished = calloc(num_threads, sizeof(int));
    for(i=0; i<num_threads; i++) {
        targs[i].counter = init_field_counter(num_fields, tree.num_fields);
        // give each thread a unique index, like an MPI worker
        targs[i].rank = i;
        targs[i].start = next_assign;
        targs[i].stop = next_assign + step;
        targs[i].finished = finished+i;
        next_assign += step;
    }

    // create threads
    for(i=0; i<num_threads; i++) {
        pthread_create(threads+i, NULL, twopoint_wrap, targs+i);
        started[i] = 1;
    }

    //while (next_assign < length) {
    while(1) {

        for(i=0; i<num_threads; i++) {
            if(finished[i]) {
                started[i] = 0;
                finished[i] = 0;
                pthread_join(threads[i], NULL);

                if (next_assign >= length) break;

                targs[i].start = next_assign;
                targs[i].stop = min(length, next_assign+step);
                next_assign = min(length, next_assign+step);

                started[i] = 1;
                pthread_create(threads+i, NULL, twopoint_wrap, targs+i);
            }
        }

        if (next_assign >= length) break;
        usleep(5000);
    }

    long long * results = calloc(num_fields * tree.num_fields, sizeof(long long));
    // join threads, sum the array for each thread into one array
    for(i=0; i<num_threads; i++) {
        if(started[i]) {
            pthread_join(threads[i], NULL);
        }

        for(j=0; j<num_fields*tree.num_fields; j++) {
            results[j] += ((targs[i].counter)->array)[j];
        }
        free_field_counter(targs[i].counter);
    }
    free(finished);
    return results;
}
