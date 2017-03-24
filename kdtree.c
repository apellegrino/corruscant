#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#ifndef KDTREE_H
#include "kdtree.h"
#endif


inline void swapCustomFloat(FLOAT * a, FLOAT * b)
{
    FLOAT temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

inline void swapInt(int * a, int * b)
{
    int temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

inline int min(int a, int b)
{
    return a < b ? a : b;
}

inline FLOAT norm2(FLOAT a, FLOAT b, FLOAT c)
{
    return a*a+b*b+c*c;
}

/* 
 *  perform a quicksort on the array a[left] to a[right] which performs
 *  identical operations on another array b[left] to b[right]
 */
void quick_argsort(FLOAT *a, int *b, int left, int right)
{
    int i, j;
    FLOAT pivot;

    if( right <= left ) return;

    pivot = a[left];
    i = left;
    j = right+1;

    while(1) {
        do ++i; while( a[i] <= pivot && i <= right );
        do --j; while( a[j] > pivot );
        if( i >= j ) break;
        swapCustomFloat(a+i,a+j);
        swapInt(b+i,b+j);
    }
    swapCustomFloat(a+left,a+j);
    swapInt(b+left,b+j);
    quick_argsort(a,b,left,j-1);
    quick_argsort(a,b,j+1,right);
}

int * argsort(FLOAT *a, int size)
{

    /* copy a to keep a unchanged */
    FLOAT *acpy = (FLOAT *) malloc(sizeof(FLOAT)*size);
    memcpy(acpy, a, sizeof(FLOAT)*size);

    int *ind = (int*) malloc(sizeof(int)*size);
    int i;

    /* initialize array of indices to a */
    for(i=0; i<size; i++)
        ind[i] = i;

    /*
     * sort the copy of "a" while performing duplicate operations on "ind".
     * "ind" becomes an array of indices to "a" which sort "a"
     */
    quick_argsort(acpy,ind, 0, size-1);
    free(acpy);
    return ind;
}

void partition(FLOAT * vals, int * args, FLOAT key, int left, int right, 
                int median, int med_index)
{
    int i, j, k, arg, size, split;
    split = med_index - left;
    j = 0;
    k = split+1;
    size = right-left+1;
    int *args_new = (int*) malloc(sizeof(int)*size);

    args_new[split] = -1;
    for(i=left; i<=right; i++) {
        arg = args[i];
        if (arg == median) continue;
        if (vals[arg] < key) {
            args_new[j] = arg;
            j++;
        } else if (vals[arg] > key) {
            args_new[k] = arg;
            k++;
        }
    }

    memcpy(args+left,args_new,size*sizeof(int));
}

node_t * build(FLOAT *x, FLOAT *y, FLOAT *z,int *x_arg, int *y_arg, int *z_arg,
                                            int left, int right, enum dim d)
{
    d=d%3;

    int med, med_arg;
    node_t *parent = (node_t *) malloc( sizeof(node_t) );
    /* Median index of the sub-array. Rounds down for even sized lists */
    med = (left+right)/2;

    /* Find index of the median in appropriate position list */
    switch(d) {
        case X:
            med_arg = x_arg[med];
            break;
        case Y:
            med_arg = y_arg[med];
            break;
        case Z:
            med_arg = z_arg[med];
            break;
    }

    /* this node is the median */
    parent->x = x[med_arg];
    parent->y = y[med_arg];
    parent->z = z[med_arg];

    parent->lchild = NULL;
    parent->rchild = NULL;

    /* 
     * Base cases: subtree of size 1, 2 or 3
     */

    switch(right-left) {
    case 0: /* array length 1, no children */
        return parent;
    case 1: /* array length 2, one child */
        parent->rchild = build(x,y,z,x_arg,y_arg,z_arg,right,right,d);
        return parent;
    case 2: /* array length 3, two children */
        parent->lchild = build(x,y,z,x_arg,y_arg,z_arg,left,left,d);
        parent->rchild = build(x,y,z,x_arg,y_arg,z_arg,right,right,d);
        return parent;
    }

    /*
     * Recursive cases: size > 2
     */

    /*  partition index array of other dims w.r.t. current dim */
    switch(d) {
    case X:
        partition(x, y_arg, x[med_arg], left, right, med_arg, med);
        partition(x, z_arg, x[med_arg], left, right, med_arg, med);
        break;
    case Y:
        partition(y, z_arg, y[med_arg], left, right, med_arg, med);
        partition(y, x_arg, y[med_arg], left, right, med_arg, med);
        break;
    case Z:
        partition(z, x_arg, z[med_arg], left, right, med_arg, med);
        partition(z, y_arg, z[med_arg], left, right, med_arg, med);
        break;
    }

    parent->lchild = build(x,y,z,x_arg,y_arg,z_arg,left,med-1,d+1);
    parent->rchild = build(x,y,z,x_arg,y_arg,z_arg,med+1,right,d+1);

    return parent;
}


void destroy(node_t *p)
{
    if(p->lchild != NULL) destroy(p->lchild);
    if(p->rchild != NULL) destroy(p->rchild);
    free(p);
}

/*
 * Query how many points in the tree with head p lie within radius r of point
 * (x, y, z). Recursive.
 */
int radius(node_t *p, enum dim d, FLOAT x, FLOAT y, FLOAT z, FLOAT r)
{

    d=d%3;

    int i;
    FLOAT rsq, dx, dy, dz;
    FLOAT pos_upper, pos_lower, point;

    rsq = r*r;
    dx = p->x - x;
    dy = p->y - y;
    dz = p->z - z;

    FLOAT n = norm2(dx,dy,dz);

    if(p->lchild == NULL && p->rchild == NULL) { /* no children */
        return (n < rsq);
    } else if (p->lchild == NULL) { /* one child */
        return (n < rsq) + radius(p->rchild,d+1,x,y,z,r);
    } else { /* two children */

        switch(d) {
        case X:
            pos_upper = x+r;
            pos_lower = x-r;
            point = p->x;
            break;
        case Y:
            pos_upper = y+r;
            pos_lower = y-r;
            point = p->y;
            break;
        case Z:
            pos_upper = z+r;
            pos_lower = z-r;
            point = p->z;
            break;
        }
        if (pos_upper < point) {
            i=radius(p->lchild,d+1,x,y,z,r);
        } else if (pos_lower > point) {
            i=radius(p->rchild,d+1,x,y,z,r);
        } else {
            i = (n < rsq) +
                radius(p->lchild,d+1,x,y,z,r) +
                radius(p->rchild,d+1,x,y,z,r);
        }
        return i;
    }
    
}

kdtree_t tree_construct(int size, FLOAT x[], FLOAT y[], FLOAT z[])
{

    kdtree_t tree;
    tree.size = size;
    tree.x = x; tree.y = y; tree.z = z;

    /* Argsort the inputs */
    int *x_arg, *y_arg, *z_arg;

    x_arg = argsort(x, size); y_arg = argsort(y, size); z_arg = argsort(z, size);

    tree.root = build(x,y,z,x_arg,y_arg,z_arg,0,size-1,X);

    free(x_arg); free(y_arg); free(z_arg);

    return tree;
}

typedef struct shared_args {
    kdtree_t tree;
    FLOAT *x, *y, *z;
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


void assign_idx(int this_thread, int num_threads, int n_array, int * start,
                                                                int * stop)
{

    int mpi_rank, mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
    int size = mpi_size * num_threads;
    int rank = mpi_rank * num_threads + this_thread;

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
        args->sum[rank] += radius(args->tree.root, 0, args->x[i],
                                            args->y[i], args->z[i], args->r);
    }

    return NULL;
}

/* 
 * Compute the sum of the number of data points in the tree which lie within
 * radius r of each of the (x,y,z) points in the array. Result may easily
 * exceed the size of a 32-bit int, so we return a long long.
 */
long long two_point_correlation(kdtree_t tree, FLOAT x[], FLOAT y[],
                FLOAT z[], int n, FLOAT r, int num_threads, MPI_Comm comm)
{

    int i, rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    shared_args_t ss;
    ss.tree = tree;
    ss.x = x;
    ss.y = y;
    ss.z = z;
    ss.n = n;
    ss.r = r;
    ss.num_threads = num_threads;
    ss.sum = (long long *)malloc(num_threads * sizeof(long long));
    
    thread_args_t targs[num_threads];
    for(i=0; i<num_threads; i++) {
        targs[i].shared_args_p = &ss;
        targs[i].thread_rank = i;
    }

    pthread_t threads[num_threads];
    threads[0] = pthread_self();

    long long result = 0;

    double t1, t2;
    t1 = MPI_Wtime();

    for(i=1; i<num_threads; i++) 
        pthread_create(threads+i, NULL, twopoint_wrap, targs+i);

    twopoint_wrap(targs);

    for(i=1; i<num_threads; i++)
        pthread_join(threads[i], NULL);

    for (i=0; i<num_threads; i++)
        result += ss.sum[i];

    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

    t2 = MPI_Wtime();

    if(!rank) {
        printf("Time on rank 0: %f sec\n", t2 - t1);
        printf("Sum: %lld\n", result);
    }

    //return gresult;
    return result;
}
