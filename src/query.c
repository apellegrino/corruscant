#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#ifndef KDTREE_H
#define KDTREE_H
#include "kdtree.h"
#endif

static double *_x_query, *_y_query, *_z_query;
static node_t * tree_data;

static inline int min(int a, int b)
{
    return a < b ? a : b;
}

static inline FLOAT norm2(FLOAT a, FLOAT b, FLOAT c)
{
    return a*a+b*b+c*c;
}

/*
 * Query how many points in the tree with head p lie within radius r of point
 * (x, y, z). Recursive.
 */
static int radius(int pi, enum dim d, int qi, FLOAT r)
{
    int i;
    double x, y, z;
    node_t p = *(tree_data+pi);

    x = _x_query[qi]; y = _y_query[qi]; z = _z_query[qi];

    FLOAT rsq, dx, dy, dz;
    FLOAT pos_upper = 0.0, pos_lower = 0.0, point = 0.0;

    rsq = r*r;
    dx = p.x - x; dy = p.y - y; dz = p.z - z;

    FLOAT n = norm2(dx,dy,dz);
    if( !(p.flags & HAS_LCHILD) ) { /* no children */
        return (n < rsq);
    } else if ( !(p.flags & HAS_RCHILD) ) { /* one child */
        return (n < rsq) + radius(left_child(pi),next_dim(d),qi,r);
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
            i=radius(left_child(pi),next_dim(d),qi,r);
        } else if (pos_lower > point) {
            i=radius(right_child(pi),next_dim(d),qi,r);
        } else {
            i = (n < rsq) +
                radius(left_child(pi),next_dim(d),qi,r) +
                radius(right_child(pi),next_dim(d),qi,r);
        }
        return i;
    }
    
}

typedef struct shared_args {
    int n;
    int node_start;
    int node_stop;
    int num_threads;
    long long *sum;
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

    args->sum[rank] = 0;

    for(i=start; i<stop; i++) {
        args->sum[rank] += radius(1, 0, i, args->r);
    }

    return NULL;
}

/* 
 * Compute the sum of the number of data points in the tree which lie within
 * radius r of each of the (x,y,z) points in the array. Result may easily
 * exceed the size of a 32-bit int, so we return a long long.
 */
long long pair_count(kdtree_t tree, FLOAT x[], FLOAT y[], FLOAT z[],
                                int n, FLOAT r, int num_threads)
{
    //printf("In pair_count\n");
    tree_data = tree.node_data;

    _x_query = x; _y_query = y; _z_query = z;

    shared_args_t ss;
    ss.n = n;
    ss.r = r;
    ss.num_threads = num_threads;
    ss.sum = (long long *)calloc(num_threads, sizeof(long long));
    
    thread_args_t targs[num_threads];

    int i;
    for(i=0; i<num_threads; i++) {
        targs[i].shared_args_p = &ss;
        targs[i].thread_rank = i;
    }

    pthread_t threads[num_threads];
    threads[0] = pthread_self();

    long long result = 0;

    for(i=0; i<num_threads; i++) 
        pthread_create(threads+i, NULL, twopoint_wrap, targs+i);

    for(i=0; i<num_threads; i++)
        pthread_join(threads[i], NULL);

    for (i=0; i<num_threads; i++)
        result += ss.sum[i];

    //printf("Leaving pair_count\n");
    return result;
}
