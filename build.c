#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#ifndef KDTREE_H
#define KDTREE_H
#include "kdtree.h"
#endif

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
 *  Return median index of sorted array given edges. With choice of rounding to
 *  the left or right for even lengths, we choose right so that nodes with one
 *  child have a left child.
 */
static inline int median(int left, int right)
{
    return (left+right+1)/2;
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

/*
static enum dim choose_dim()
{
}
*/

FLOAT * get_data_array(kdtree_t tree, enum dim d)
{
    switch(d) {
    case X:
        return tree.x_data;
    case Y:
        return tree.y_data;
    case Z:
        return tree.z_data;
    default:
        return NULL;
    }

}

int * get_arg_array(kdtree_t tree, enum dim d)
{
    switch(d) {
    case X:
        return tree.x_arg;
    case Y:
        return tree.y_arg;
    case Z:
        return tree.z_arg;
    default:
        return NULL;
    }

}

/*
 *  Partition a section of the array "args" from indicies "left" to "right" as
 *  being either indicies to values in "vals" smaller or larger than the key.
 *  Leave a null value in the middle in place of the index that locates "key"
 *  inside "vals".
 */
static void partition(kdtree_t tree, int left, int right, enum dim d_key,
                      enum dim d_part)
{
    FLOAT * vals;
    int * args;
    vals = get_data_array(tree, d_key);
    args = get_arg_array(tree, d_part);

    int med_index = median(left,right);
    int med_arg = *(get_arg_array(tree,d_key) + med_index);
    FLOAT key = vals[med_arg];

    int split = med_index - left;

    int i, j, k, size;

    j = 0;
    k = split+1;
    size = right-left+1;
    int *args_new = (int*) malloc(sizeof(int)*size);

    args_new[split] = -1;

    int arg;
    for(i=left; i<=right; i++) {
        arg = args[i];
        if (arg == med_arg) continue;
        if (vals[arg] < key) {
            args_new[j] = arg;
            j++;
        } else if (vals[arg] > key) {
            args_new[k] = arg;
            k++;
        }
    }

    memcpy(args+left,args_new,size*sizeof(int));
    free(args_new);
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

void destroy(kdtree_t * t)
{
    free(t->node_data);
}

void build(kdtree_t tree, int ind, int left, int right, enum dim d)
{
    double *x, *y, *z;
    int med, med_arg;

    node_t * parent = tree.node_data + ind;
    parent->flags = 0;

    x = tree.x_data; y = tree.y_data; z = tree.z_data;

    /* Median index of the sub-array. Rounds up for even sized lists */
    med = median(left,right);

    /* Find index of the median in appropriate position list */
    med_arg = *( get_arg_array(tree, d) + med );

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
        partition(tree, left, right, X, Y);
        partition(tree, left, right, X, Z);
        break;
    case Y:
        partition(tree, left, right, Y, X);
        partition(tree, left, right, Y, Z);
        break;
    case Z:
        partition(tree, left, right, Z, X);
        partition(tree, left, right, Z, Y);
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
    tree.size = size;
    tree.memsize = pow2ceil(size);
    tree.node_data = (node_t *) calloc( tree.memsize, sizeof(node_t) );
    tree.x_data = x; tree.y_data = y; tree.z_data = z;

    /* Argsort the inputs */
    tree.x_arg = argsort(x, size);
    tree.y_arg = argsort(y, size);
    tree.z_arg = argsort(z, size);

    build(tree, 1, 0, size-1, X );

    /* argsort key arrays are only necessary for building the tree */
    free(tree.x_arg); free(tree.y_arg); free(tree.z_arg);

    return tree;
}
