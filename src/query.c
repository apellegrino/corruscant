#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#ifndef KDTREE_H
#define KDTREE_H
#include "kdtree.h"
#endif

static double *_x_query, *_y_query, *_z_query;
static int *_field_query;
static node_t * _tree_data;
static int _include;
static int _use_err;

static void count_node(FLOAT n, FLOAT rsq, field_counter_t * counter, int this_field, int query_field)
{
    if (n < rsq) {
        // index 0 represents whole dataset
        counter->array[0]++;

        // no error
        if(!_use_err) return;

        // field-to-field error -- increment counter for this_field if it
        // coincides with query_field
        if (_include) {
            if (this_field == query_field) (counter->array[this_field])++;
        }
        else { //jackknife error -- increment all but this_field and query_field
            int i;
            // needs vectorization?
            for(i=1; i<counter->size; i++) { (counter->array[i])++; }
            /* 
             * un-count this pair for the querying point field or that of the
             * point on the tree, but don't double count if they have the same
             * field
             */
            (counter->array[this_field])--;
            if (this_field != query_field) (counter->array[query_field])--;
        }
    }
    return;
}
static inline int min(int a, int b)
{
    return a < b ? a : b;
}

static inline FLOAT norm2(FLOAT a, FLOAT b, FLOAT c)
{
    return a*a+b*b+c*c;
}

inline int get_field(node_t p)
{
    return (p.flags & ID_MASK) >> 2;
}

inline int has_lchild(node_t p)
{
    return p.flags & HAS_LCHILD;
}

inline int has_rchild(node_t p)
{
    return p.flags & HAS_RCHILD;
}

/*
 * Query how many points in the tree with head p lie within radius r of point
 * (x, y, z). Recursive.
 */
static void radius(int pi, enum dim d, int qi, FLOAT r, field_counter_t * counter)
{
    //int i;
    FLOAT x, y, z;
    node_t p = *(_tree_data+pi);

    int this_field = get_field(p);
    
    x = (FLOAT) _x_query[qi]; y = (FLOAT) _y_query[qi]; z = (FLOAT) _z_query[qi];
    int query_field = _field_query[qi];

    FLOAT rsq, dx, dy, dz;
    FLOAT pos_upper = 0.0, pos_lower = 0.0, point = 0.0;

    rsq = r*r;
    dx = p.x - x; dy = p.y - y; dz = p.z - z;

    FLOAT n = norm2(dx,dy,dz);
    if( !has_lchild(p) ) { /* no children */

        count_node(n, rsq, counter, this_field, query_field);
        return;

    } else if ( !has_rchild(p) ) { /* one child */

        count_node(n, rsq, counter, this_field, query_field);
        radius(left_child(pi),next_dim(d),qi,r,counter);
        return;

    } else { /* two children */
        switch(d) {
        case X:
            pos_upper = x+r;
            pos_lower = x-r;
            point = p.x;
            break;
        case Y:
            pos_upper = y+r;
            pos_lower = y-r;
            point = p.y;
            break;
        case Z:
            pos_upper = z+r;
            pos_lower = z-r;
            point = p.z;
            break;
        }

        if (pos_upper < point) {
            radius(left_child(pi),next_dim(d),qi,r,counter);
        } else if (pos_lower > point) {
            radius(right_child(pi),next_dim(d),qi,r,counter);
        } else {
            count_node(n, rsq, counter, this_field, query_field);
            radius(left_child(pi),next_dim(d),qi,r,counter);
            radius(right_child(pi),next_dim(d),qi,r,counter);
        }
        return;
    }
    
}

typedef struct shared_args {
    int n;
    int node_start;
    int node_stop; // actually one after the last index to stop at
    int num_fields;
    int num_threads;
    field_counter_t ** counters;
    FLOAT r;
} shared_args_t;

typedef struct thread_args {
    int thread_rank;
    struct shared_args * shared_args_p;
} thread_args_t;


void assign_idx(int rank, int size, int n_array, int * start, int * stop)
{
    int n_array_local = n_array / size + ( n_array % size > rank );

    *start = rank * ( n_array / size ) + min( rank, n_array % size );
    *stop = *start + n_array_local;
}

void * twopoint_wrap(void *voidargs)
{

    int start, stop, i;
    thread_args_t *targs;
    targs = (thread_args_t *) voidargs;
    shared_args_t * args;
    args = targs->shared_args_p;
    
    int rank = targs->thread_rank;
    assign_idx(rank, args->num_threads, args->n, &start, &stop);

    field_counter_t * new_counter = malloc(sizeof(field_counter_t));
    new_counter->array = calloc(args->num_fields+1, sizeof(long long));
    new_counter->size = args->num_fields+1;
    args->counters[rank] = new_counter;
   
    for(i=start; i<stop; i++) {
        radius(1, 0, i, args->r, args->counters[rank]);
    }

    return NULL;
}

/* 
 * Compute the sum of the number of data points in the tree which lie within
 * radius r of each of the (x,y,z) points in the array. Result may easily
 * exceed the size of a 32-bit int, so we return a long long.
 */
static long long * pair_count(kdtree_t tree, array3d_t data, double r, int num_threads)
{
    int i, j;
    _tree_data = tree.node_data;
    _x_query = data.x; _y_query = data.y; _z_query = data.z;
    _field_query = data.fields;
    int num_fields = data.num_fields;

    if(num_threads > 1<<16) {
        printf("More than %d threads!\n", 1<<16);
        exit(1);
    }
    
    shared_args_t ss;
    ss.counters = (field_counter_t **) malloc(num_threads * sizeof(field_counter_t *));
    ss.n = data.size;
    ss.r = r;
    ss.num_fields = num_fields;
    ss.num_threads = num_threads;
    
    thread_args_t targs[num_threads];

    // give each thread a unique index, like an MPI worker
    for(i=0; i<num_threads; i++) {
        targs[i].shared_args_p = &ss;
        targs[i].thread_rank = i;
    }

    pthread_t threads[num_threads];

    long long * results = calloc(num_fields+1, sizeof(long long));

    // create threads
    for(i=0; i<num_threads; i++) 
        pthread_create(threads+i, NULL, twopoint_wrap, targs+i);

    // join threads, sum the array for each thread into one array
    for(i=0; i<num_threads; i++) {
        pthread_join(threads[i], NULL);
        for(j=0; j<num_fields+1; j++) {
            results[j] += ((ss.counters[i])->array)[j];
        }
        free(ss.counters[i]);
    }
    return results;
}

long long * pair_count_jackknife(kdtree_t tree, array3d_t data,
                               double r, int num_threads)
{
    // exclude only the field of the querying point
    _include = 0;
    _use_err = 1;
    return pair_count(tree, data, r, num_threads);
}

long long * pair_count_ftf(kdtree_t tree, array3d_t data,
                         double r, int num_threads)
{
    // include only the field of the querying point
    _include = 1;
    _use_err = 1;
    return pair_count(tree, data, r, num_threads);
}

long long * pair_count_noerr(kdtree_t tree, array3d_t data,
                            double r, int num_threads)
{
    _use_err = 0;
    return pair_count(tree, data, r, num_threads);
}
