/*  Andrew Pellegrino, 2017
 *
 *  K-d tree building algorithm adapted from Russell A. Brown, Building a
 *  Balanced k-d Tree in O(kn log n) Time, Journal of Computer Graphics
 *  Techniques (JCGT), vol. 4, no. 1, 50-68, 2015
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef KDTREE_H
#define KDTREE_H
#include "kdtree.h"
#endif

#define ROOT 1

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
 *  perform a merge sort on the array a[left] to a[right] on doubles,
 *  which performs identical operations on another array b[left] to b[right]
 *  on ints
 */
static void merge_argsort(double *a, int *b, int left, int right)
{
    if (right == left) return;
    if (right - left == 1) {
        if(a[right] < a[left]) {
            swapDouble(a+left,a+right);
            swapInt(b+left,b+right);
            return;
        }
    }

    int m = median(left,right);

    merge_argsort(a, b, left, m-1);
    merge_argsort(a, b, m, right);

    double * new_a = (double *) malloc( (right-left+1) * sizeof(double) );
    int * new_b = (int *) malloc( (right-left+1) * sizeof(int) );

    // cursors to current elements in each array
    int lc = left;
    int rc = m;
    int nc = 0;

    while(lc < m && rc <= right) {
        if( a[lc] < a[rc] ) {
            new_a[nc] = a[lc];
            new_b[nc] = b[lc];
            lc++; nc++;
        } else if ( a[lc] > a[rc] ) {
            new_a[nc] = a[rc];
            new_b[nc] = b[rc];
            rc++; nc++;
        } else {
            new_a[nc] = a[lc];
            new_b[nc] = b[lc];
            lc++; nc++;
            new_a[nc] = a[rc];
            new_b[nc] = b[rc];
            rc++; nc++;
        }
    }

    while(lc < m) {
        new_a[nc] = a[lc];
        new_b[nc] = b[lc];
        lc++; nc++;
    }
    while(rc <= right) {
        new_a[nc] = a[rc];
        new_b[nc] = b[rc];
        rc++; nc++;
    }

    int i;
    for(i=0; i<right-left+1; i++) {
        a[left+i] = new_a[i];
        b[left+i] = new_b[i];
    }

    free(new_a);
    free(new_b);

    return;
}

/*
 *  Return a pointer to an array of ints which index the input array in
 *  ascending order of values
 */
static int * argsort(double *a, int size)
{

    /* copy "a" to keep "a" unchanged */
    double *acpy = (double *) malloc(sizeof(double)*size);
    memcpy(acpy, a, sizeof(double)*size);

    int * ind = (int*) malloc(sizeof(int)*size);
    int i;

    /* initialize array of indices */
    for(i=0; i<size; i++)
        ind[i] = i;

    /*
     * sort the copy of "a" while performing duplicate operations on "ind"
     */
    merge_argsort(acpy, ind, 0, size-1);
    free(acpy);
    return ind;
}

//static inline double * get_data_array(kdtree_t tree, enum dim d)
double * get_data_array(kdtree_t tree, enum dim d)
{
    switch(d) {
    case X:
        return tree.x;
    case Y:
        return tree.y;
    case Z:
        return tree.z;
    default:
        return NULL;
    }

}

//static inline int * get_arg_array(kdtree_t tree, enum dim d)
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
static void partition(kdtree_t * tree, int left, int right, enum dim d_key,
                      enum dim d_part)
{
    double * vals;
    int * args;
    vals = get_data_array(*tree, d_key);
    args = get_arg_array(*tree, d_part);

    int med_index = median(left,right);
    int med_arg = (get_arg_array(*tree,d_key))[med_index];
    double key = vals[med_arg];

    int i, size;
    int split = med_index - left;
    size = right-left+1;

    int *args_new = (int*) malloc(sizeof(int)*size);

    // this value should not be used again
    args_new[split] = -1;

    // partition `args` on the range `left` to `right` compared to the key
    int arg;
    int j = 0;
    int k = split+1;

    for(i=left; i<=right; i++) {

        arg = args[i];

        if (arg == med_arg) continue;
        if (vals[arg] < key) {
            args_new[j] = arg;
            j++;
        } else if (vals[arg] > key) {
            args_new[k] = arg;
            k++;
        } else {
            // edge case of float equality. Partition cannot handle
            // duplicate values
            fprintf(stderr, "Two equal floating point numbers found!\n");
            exit(1);
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
}

static inline void set_id(node_t * node, int id)
{
    node->flags |= ID_MASK & (id << 2);
}

/* Flag node as having a left child */
static inline void set_lchild(node_t * node)
{
    node->flags |= HAS_LCHILD;
}

/* Flag node as having a right child */
static inline void set_rchild(node_t * node)
{
    node->flags |= HAS_RCHILD;
}

void build(kdtree_t * tree, int ind, int left, int right, enum dim d)
{
    double *x, *y, *z;
    int *ids;
    int med, med_arg, this_id;

    node_t * parent = tree->node_data + ind;
    parent->flags = 0;

    x = tree->x; y = tree->y; z = tree->z;

    /* Median index of the sub-array. Rounds up for even sized lists */
    med = median(left,right);

    /* Find index of the median in appropriate position list */
    med_arg = (get_arg_array(*tree, d))[med];

    /* the median point in dim d has index med_arg */
    parent->x = x[med_arg]; parent->y = y[med_arg]; parent->z = z[med_arg];

    ids = tree->fields;
    if (ids != NULL) {
        this_id = ids[med_arg];
        set_id(parent, this_id);
    }

    /* 
     *  Base cases: subtree of size 1, 2 or 3
     *  Skip partitioning step
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
     *  Recursive cases: size > 2
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

/*
argarray3d_t array3d_argsort(array3d_t data)
{
    argarray3d_t arg;
    arg.x = argsort(data.x, data.size);
    arg.y = argsort(data.y, data.size);
    arg.z = argsort(data.z, data.size);
    arg.size = data.size;
    return arg;
}
*/

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

/*
{
    kdtree_t tree;
    tree.size = data.size;
    tree.memsize = pow2ceil(data.size);
    tree.node_data = (node_t *) calloc( tree.memsize, sizeof(node_t) );

    tree.data = data;

    tree.arg_data = array3d_argsort(tree.data);

    build(&tree, ROOT, 0, data.size-1, X );

    // argsort key arrays are only necessary for building the tree
    free(tree.arg_data.x);
    free(tree.arg_data.y);
    free(tree.arg_data.z);

    return tree;
}
*/

kdtree_t tree_construct(double * x, double * y, double * z, int * fields, int length, int num_fields)
{
    kdtree_t tree;
    tree.size = length;
    tree.memsize = pow2ceil(length);
    tree.node_data = (node_t *) calloc( tree.memsize, sizeof(node_t) );

    tree.x = x;
    tree.y = y;
    tree.z = z;
    tree.fields = fields;
    tree.num_fields = num_fields;

    /* Argsort the inputs */
    tree.x_arg = argsort(x, length);
    tree.y_arg = argsort(y, length);
    tree.z_arg = argsort(z, length);

    build( &tree, ROOT, 0, length-1, X );

    /* argsort key arrays are only necessary for building the tree */
    free(tree.x_arg);
    free(tree.y_arg);
    free(tree.z_arg);

    return tree;
}
