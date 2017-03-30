#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#ifndef KDTREE_H
#define KDTREE_H
#include "kdtree.h"
#endif

static int *x_arg, *y_arg, *z_arg;

static inline void swapCustomFloat(FLOAT * a, FLOAT * b)
{
    FLOAT temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

static inline void swapInt(int * a, int * b)
{
    int temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

static inline int min(int a, int b)
{
    return a < b ? a : b;
}

/* 
 *  perform a quicksort on the array a[left] to a[right] which performs
 *  identical operations on another array b[left] to b[right]
 */
static void quick_argsort(FLOAT *a, int *b, int left, int right)
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

static int * argsort(FLOAT *a, int size)
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
    quick_argsort(acpy, ind, 0, size-1);
    free(acpy);
    return ind;
}

static void partition(FLOAT * vals, int * args, FLOAT key, int left, int right, 
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

inline int left_child(int p)
{
    return 2*p;
}

inline int right_child(int p)
{
    return 2*p+1;
}

enum dim inline next_dim(enum dim d)
{
    return (d+1)%3;
}

void build(kdtree_t tree, int ind, int left, int right, enum dim d)
{
    double *x, *y, *z;
    int med, med_arg;

    node_t * parent = tree.node_data + ind;
    parent->flags = 0;

    x = tree.x_data; y = tree.y_data; z = tree.z_data;

    /* Median index of the sub-array. Rounds up for even sized lists */
    med = (left+right+1)/2;

    med_arg = 0;
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
    parent->x = x[med_arg]; parent->y = y[med_arg]; parent->z = z[med_arg];

    /* 
     * Base cases: subtree of size 1, 2 or 3
     */

    switch(right-left) {
    case 0: /* array length 1, no children */
        return;
    case 1: /* array length 2, one child */
        parent->flags |= HAS_LCHILD;
        build(tree, left_child(ind), left,left,d);
        return;
    case 2: /* array length 3, two children */
        parent->flags |= HAS_LCHILD;
        parent->flags |= HAS_RCHILD;
        build(tree, left_child(ind), left,left,d);
        build(tree, right_child(ind), right,right,d);
        return;
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

    parent->flags |= HAS_LCHILD;
    parent->flags |= HAS_RCHILD;
    build(tree, left_child(ind), left,med-1,next_dim(d));
    build(tree, right_child(ind), med+1,right,next_dim(d));
    return;
}

static int pow2ceil(int x)
{
    int i = 1;
    while(x > 1) {
        x /= 2;
        i++;
    }

    while(i > 0) {
        x *= 2;
        i--;
    }

    return x;
}

kdtree_t tree_construct(int size, FLOAT x[], FLOAT y[], FLOAT z[])
{
    kdtree_t tree;
    tree.memsize = pow2ceil(size);
    tree.node_data = (node_t *) calloc( tree.memsize, sizeof(node_t) );
    tree.x_data = x; tree.y_data = y; tree.z_data = z;

    /* Argsort the inputs */
    x_arg = argsort(x, size); y_arg = argsort(y, size); z_arg = argsort(z, size);

    build(tree, 1, 0, size-1, X );

    free(x_arg); free(y_arg); free(z_arg);

    return tree;
}

