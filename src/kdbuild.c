/*  Andrew Pellegrino, 2017
 *
 *  K-d tree building algorithm adapted from Russell A. Brown, Building a
 *  Balanced k-d Tree in O(kn log n) Time, Journal of Computer Graphics
 *  Techniques (JCGT), vol. 4, no. 1, 50-68, 2015
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "kdtree.h"

#define ROOT 1

static inline void swapDatum(datum_t * a, datum_t * b)
{
    datum_t temp;
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

/*
 *  Return median index of sorted array given edges. With choice of rounding to
 *  the left or right for even lengths, we choose right so that nodes with one
 *  child have a left child.
 */
static inline int median(int left, int right)
{
    return (left+right+1)/2;
}

// Return 1 if superkey of a < that of b, return 0 if superkey of a > that of b
// i.e. whether the elements are in order
int superkey_compare(datum_t * data, int msd, int a, int b)
{
    int i, dim;
    double aval, bval;

    for(i=msd; i<msd+NDIM; i++) {
        dim = i < NDIM ? i : i - NDIM;
        aval = (data[a]).value[dim];
        bval = (data[b]).value[dim];
        if(aval == bval) continue;

        return aval < bval;
    }

    printf("ERROR: Duplicate points found %d %d\n", a, b);
    for(dim=0; dim<NDIM; dim++) {
        aval = (data[a]).value[dim];
        bval = (data[b]).value[dim];
        printf("DIM %d : %f %f\n", dim, aval, bval);
    }
    exit(EXIT_FAILURE);
} 

int superkey_dupesafe(datum_t * data, int a, int b)
{
    int i;
    double aval, bval;

    for(i=0; i<NDIM; i++) {
        aval = (data[a]).value[i];
        bval = (data[b]).value[i];
        if(aval == bval) continue;

        return aval < bval ? 1 : -1;
    }
    
    return 0;
} 

static void mergesort_x(datum_t * data, int left, int right)
{
    if(right == left) return;
    if(right - left == 1) {
        if(superkey_dupesafe(data, left, right) == -1) {
            swapDatum(data+left, data+right);
        }
        return;
    }

    int m = median(left, right);

    mergesort_x(data, left, m-1);
    mergesort_x(data, m, right);

    datum_t * new_data = (datum_t *) malloc((right-left+1) * sizeof(datum_t));

    int lc = left;
    int rc = m;
    int nc = 0;

    while(lc < m && rc <= right) {
        if(superkey_dupesafe(data, lc, rc) == 1) { // if data[lc] < data[rc]
            new_data[nc] = data[lc];
            lc++; nc++;
        } else if(superkey_dupesafe(data, lc, rc) == -1) {
            new_data[nc] = data[rc];
            rc++; nc++;
        } else { // equal points, i.e. duplicates
            new_data[nc] = data[lc];
            lc++; nc++;
            new_data[nc] = data[rc];
            rc++; nc++;
            
        }
    }

    while(lc < m) {
        new_data[nc] = data[lc];
        lc++; nc++;
    }
    while(rc <= right) {
        new_data[nc] = data[rc];
        rc++; nc++;
    }

    int i;
    for(i=0; i<right-left+1; i++) {
        data[left+i] = new_data[i];
    }

    free(new_data);

    return;
}
/* 
 *  perform a merge sort on the array data[left] to data[right],
 *  which performs identical operations on another array ind[left] to ind[right]
 *  on integer indices to the data array
 */
static void merge_supersort(datum_t * data, int * ind, int left, int right, int dim)
{
    if (right == left) return;
    if (right - left == 1) {
        if(!superkey_compare(data, dim, left, right)) {
            swapDatum(data+left,data+right);
            swapInt(ind+left,ind+right);
        }
        return;
    }

    int m = median(left,right);

    merge_supersort(data, ind, left, m-1, dim);
    merge_supersort(data, ind, m, right, dim);

    datum_t * new_data = (datum_t *) malloc( (right-left+1) * sizeof(datum_t) );
    int * new_ind = (int *) malloc( (right-left+1) * sizeof(int) );

    // cursors to current elements in each array
    int lc = left;
    int rc = m;
    int nc = 0;

    while(lc < m && rc <= right) {
        if(superkey_compare(data, dim, lc, rc)) { // if data[lc] < data[rc]
            new_data[nc] = data[lc];
            new_ind[nc] = ind[lc];
            lc++; nc++;
        } else { // comparision should never be equal
            new_data[nc] = data[rc];
            new_ind[nc] = ind[rc];
            rc++; nc++;
        }
    }

    while(lc < m) {
        new_data[nc] = data[lc];
        new_ind[nc] = ind[lc];
        lc++; nc++;
    }
    while(rc <= right) {
        new_data[nc] = data[rc];
        new_ind[nc] = ind[rc];
        rc++; nc++;
    }

    int i;
    for(i=0; i<right-left+1; i++) {
        data[left+i] = new_data[i];
        ind[left+i] = new_ind[i];
    }

    free(new_data);
    free(new_ind);

    return;
}

/*
 *  Return a pointer to an array of ints which index the input array in
 *  ascending order of values
 */
static int * supersort(datum_t * data, int size, int dim)
{
    datum_t * cpy = (datum_t *) malloc(sizeof(datum_t) * size);
    memcpy(cpy, data, sizeof(datum_t) * size);

    int * indices = (int *) malloc(sizeof(int) * size);
    int i;

    /* Initialize array of indices */
    for(i=0; i<size; i++) indices[i] = i;

    merge_supersort(cpy, indices, 0, size-1, dim);

    for(i=0; i<size; i++) {
        if(indices[i] < 0 || indices[i] >= size) {
            printf("ERROR: arg array not initialized properly\n");
            exit(EXIT_FAILURE);
        }
    }

    free(cpy);
    return indices;
}

/*
 *  Partition a section of the array `args` from indices `left` to `right` as
 *  being either indicies to values in `vals` smaller or larger than the key.
 *  Leave a null value in the middle in place of the index that locates `key`
 *  inside `vals`.
 */
static void partition(kdtree_t * tree, int left, int right, int d_key)
{
    int med_index = median(left,right);
    int med_arg = (tree->args[d_key])[med_index];

    int i, j, k, dim, size, arg;
    size = right-left+1;

    int * args;
    int * args_new = (int *) malloc(sizeof(int) * size);

    for(dim=0; dim<NDIM; dim++) {
        if(dim == d_key) continue;

        // relative index of left
        j = 0;
        // relative index of median + 1
        k = med_index-left+1;

        args = tree->args[dim];

        for(i=left; i<=right; i++) {
            arg = args[i];

            if(arg == med_arg) continue;

            // if data[arg] < data[med_arg]
            if( superkey_compare(tree->data, d_key, arg, med_arg) ) {
                args_new[j] = arg;
                j++;
            } else {
                args_new[k] = arg;
                k++;
            }

        }

        memcpy(args+left, args_new, size*sizeof(int));
    }
    free(args_new);
}

static inline int left_child(int p)
{
    return 2*p;
}

static inline int right_child(int p)
{
    return 2*p+1;
}

void destroy(kdtree_t t)
{
    free(t.node_data);
}

int max_variance_dim(kdtree_t * tree, datum_t * data, int begin, int end)
{
    double avg[NDIM]  = {0};
    double var[NDIM]  = {0};
    double temp[NDIM] = {0};

    int i, j, arg;
    int * args;

    // Maybe reverse loop nesting later
    for(j=0; j<NDIM; j++) {
        args = tree->args[j];

        for(i=begin; i<=end; i++) {
            arg = args[i];
            avg[j] += (data[arg]).value[j];
        }
    }

    for(j=0; j<NDIM; j++) {
        avg[j] /= (end-begin+1);
    }

    for(j=0; j<NDIM; j++) {
        args = tree->args[j];
        for(i=begin; i<=end; i++) {
            arg = args[i];
            temp[j] = (data[arg].value[j] - avg[j]);
            var[j] += temp[j]*temp[j];
        }
    }
    
    int dim = 0;
    double maxvar = var[0];

    for(j=1; j<NDIM; j++) {
        if(var[j] > maxvar) {
            maxvar = var[j];
            dim = j;
        }
    }

    return dim;
}

double arg_variance(int * args, double * data, int begin, int end)
{
    double avg = 0;
    double var = 0;
    int i;
    for(i=begin; i<=end; i++) {
        avg += data[args[i]];
    }

    avg /= (end-begin+1);

    for(i=begin; i<=end; i++) {
        var += (data[args[i]] - avg)*(data[args[i]] - avg);
    }

    // we do not need a coefficient in front, because only data with same
    // N will be compared
    return var;
}

void build(kdtree_t * tree, int ind, int left, int right)
{
    int med, med_arg;

    node_t * parent = tree->node_data + ind;

    /* Median index of the sub-array. Rounds up for even sized lists */
    med = median(left,right);

    /* Choose dim with max. variance */
    int tdim = 0;
    if(right-left > 2) {
        tdim = max_variance_dim(tree, tree->data, left, right);
    }

    parent->dim = tdim;

    /* Find index of the median in appropriate position list */
    med_arg = (tree->args[tdim])[med];

    /* the median point in dim tdim has index med_arg */
    parent->data = tree->data[med_arg];

    /* 
     *  Base cases: subtree of size 1, 2 or 3
     *  Skip partitioning step
     *  TODO: make true base cases
     */

    switch(right-left) {
    case 0: /* array length 1, no children */
        parent->has_lchild = 0;
        parent->has_rchild = 0;
        return;
    case 1: /* array length 2, one child */
        parent->has_lchild = 1;
        parent->has_rchild = 0;
        build(tree, left_child(ind), left, left);
        return;
    case 2: /* array length 3, two children */
        parent->has_lchild = 1;
        parent->has_rchild = 1;
        build(tree, left_child(ind), left, left);
        build(tree, right_child(ind), right, right);
        return;
    }

    /*
     *  Recursive cases: size > 2
     */

    /*  partition index array of other dims w.r.t. current dim */
    partition(tree, left, right, tdim);

    parent->has_lchild = 1;
    parent->has_rchild = 1;
    build(tree, left_child(ind), left, med-1);
    build(tree, right_child(ind), med+1, right);

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

kdtree_t tree_construct(double * xyz, int * fields, int length, int num_fields)
{
    // replace once dupe checking is done
    //kdtree_t tree;
    //tree.size = length;
    //tree.memsize = pow2ceil(length);
    //tree.node_data = (node_t *) calloc(tree.memsize, sizeof(node_t));

    int i;

    // copy unchecked data
    datum_t * uniquedata = (datum_t *) malloc(length * sizeof(datum_t));

    for(i=0; i<length; i++) {
        uniquedata[i].value[0] = xyz[3*i];
        uniquedata[i].value[1] = xyz[3*i+1];
        uniquedata[i].value[2] = xyz[3*i+2];
        uniquedata[i].field = 0;
        uniquedata[i].weight = 1;
    }

    if(fields != NULL) {
        for(i=0; i<length; i++) {
            uniquedata[i].field = fields[i];
        }

    }        

    // mergesort data, 0th dim as most significant
    mergesort_x(uniquedata, 0, length-1);

    /* remove duplicates */

    int old = 0;
    int drop = 1;
    int new;

    for(new=1; new<length; new++) {
        if(superkey_dupesafe(uniquedata, old, new) == 0) {
            uniquedata[old].weight++;
            continue;
        }

        uniquedata[drop] = uniquedata[new];
        old++; drop++;
    }

    int newlength = old+1;
    uniquedata = realloc(uniquedata, newlength * sizeof(datum_t));

    int sumweights = 0;
    for(i=0; i<newlength; i++) {
        sumweights += uniquedata[i].weight;
    }

    //printf("%d duplicates consolidated\n", length-newlength);

    if(sumweights != length) {
        printf("Invalid weights: %d out of %d\n", sumweights, length);
        exit(EXIT_FAILURE);
    }

    kdtree_t tree;
    tree.size = newlength;
    tree.memsize = pow2ceil(newlength);
    tree.node_data = (node_t *) calloc(tree.memsize, sizeof(node_t));

    tree.data = uniquedata;
    //tree.fields = fields;
    tree.num_fields = num_fields;

    /* Argsort the inputs */
    int dim;
    for(dim=0; dim<NDIM; dim++) {
        tree.args[dim] = supersort(tree.data, newlength, dim);
    }

    build( &tree, ROOT, 0, newlength-1 );

    /* arg arrays are only necessary for building the tree */
    for(dim=0; dim<NDIM; dim++) {
        free(tree.args[dim]);
    }
    free(uniquedata);

    return tree;
}
