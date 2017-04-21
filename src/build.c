/*  Andrew Pellegrino, 2017
 *
 *  K-d tree building algorithm adapted from Russell A. Brown, Building a
 *  Balanced k-d Tree in O(kn log n) Time, Journal of Computer Graphics
 *  Techniques (JCGT), vol. 4, no. 1, 50-68, 2015
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#ifndef KDTREE_H
#define KDTREE_H
#include "kdtree.h"
#endif

static inline void swapDouble(double * a, double * b)
{
    double temp;
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
static void quick_argsort(double *a, int *b, int left, int right)
{
    int i, j;
    double pivot;

    if( right <= left ) return;

    pivot = a[left];
    i = left;
    j = right+1;

    while(1) {
        do ++i; while( a[i] <= pivot && i <= right );
        do --j; while( a[j] > pivot );
        if( i >= j ) break;
        swapDouble(a+i,a+j);
        swapInt(b+i,b+j);
    }
    swapDouble(a+left,a+j);
    swapInt(b+left,b+j);
    quick_argsort(a,b,left,j-1);
    quick_argsort(a,b,j+1,right);
}

static int * argsort(double *a, int size)
{

    /* copy a to keep a unchanged */
    double *acpy = (double *) malloc(sizeof(double)*size);
    memcpy(acpy, a, sizeof(double)*size);

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

double * get_data_array(kdtree_t tree, enum dim d)
{
    switch(d) {
    case X:
        return tree.data.x;
    case Y:
        return tree.data.y;
    case Z:
        return tree.data.z;
    default:
        return NULL;
    }

}

int * get_arg_array(kdtree_t tree, enum dim d)
{
    switch(d) {
    case X:
        return tree.arg_data.x;
    case Y:
        return tree.arg_data.y;
    case Z:
        return tree.arg_data.z;
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
    double * vals;
    int * args;
    vals = get_data_array(tree, d_key);
    args = get_arg_array(tree, d_part);

    int med_index = median(left,right);
    int med_arg = *(get_arg_array(tree,d_key) + med_index);
    double key = vals[med_arg];

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

void destroy(kdtree_t t)
{
    free(t.node_data);

    // Have these already been freed?
    //free(t.data.x);
    //free(t.data.y);
    //free(t.data.z);
}

inline void set_id(node_t * node, int id)
{
    node->flags |= (ID_MASK & id << 2);
}

inline void set_lchild(node_t * node)
{
    node->flags |= HAS_LCHILD;
}

inline void set_rchild(node_t * node)
{
    node->flags |= HAS_RCHILD;
}

void build(kdtree_t tree, int ind, int left, int right, enum dim d)
{
    double *x, *y, *z;
    int *ids;
    int med, med_arg, this_id;

    node_t * parent = tree.node_data + ind;
    parent->flags = 0;

    x = tree.data.x; y = tree.data.y; z = tree.data.z;

    /* Median index of the sub-array. Rounds up for even sized lists */
    med = median(left,right);

    /* Find index of the median in appropriate position list */
    med_arg = *( get_arg_array(tree, d) + med );

    /* this node is the median */
    parent->x = x[med_arg]; parent->y = y[med_arg]; parent->z = z[med_arg];

    ids = tree.field_data;
    if (ids != NULL) {
        this_id = ids[med_arg];
        set_id(parent, this_id);
    }

    /* 
     * Base cases: subtree of size 1, 2 or 3
     * Skip partitioning step
     */

    switch(right-left) {
    case 0: /* array length 1, no children */
        return;
    case 1: /* array length 2, one child */
        set_lchild(parent);
        build(tree, left_child(ind), left,left,d);
        return;
    case 2: /* array length 3, two children */
        set_lchild(parent);
        set_rchild(parent);
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

    set_lchild(parent);
    set_rchild(parent);
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

array3d_t form_array(double *x, double *y, double *z, int size)
{
    array3d_t data;
    data.x = x;
    data.y = y;
    data.z = z;
    data.size = size;
    return data;
}

argarray3d_t array3d_argsort(array3d_t data)
{
    argarray3d_t arg;
    arg.x = argsort(data.x, data.size);
    arg.y = argsort(data.y, data.size);
    arg.z = argsort(data.z, data.size);
    arg.size = data.size;
    return arg;
}

kdtree_t tree_construct(array3d_t data, int * fid)
{
    kdtree_t tree;
    tree.size = data.size;
    tree.memsize = pow2ceil(data.size);
    tree.node_data = (node_t *) calloc( tree.memsize, sizeof(node_t) );

    tree.data = data;
    tree.field_data = fid;

    /* Argsort the inputs */
    tree.arg_data = array3d_argsort(tree.data);

    build(tree, 1, 0, data.size-1, X );

    /* argsort key arrays are only necessary for building the tree */
    free(tree.arg_data.x);
    free(tree.arg_data.y);
    free(tree.arg_data.z);

    return tree;
}
