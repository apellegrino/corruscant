//#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include "mpi.h"
#include "kdtree.h"

static inline void swapDouble(double * a, double * b) {
	double temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

static inline void swapInt(int * a, int * b) {
	int temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

static inline double max(double a, double b) {
	return a > b ? a : b;
}

static inline double min(double a, double b) {
	return a < b ? a : b;
}

static inline double max3(double a, double b, double c) {
	if(a < b) {
		if(a < c) { return a; }
		else { return c; }
	} else {
		if(b < c) { return b; }
		else { return c; }
	}
}

static inline double min3(double a, double b, double c) {
	if(a > b) {
		if(a > c) { return a; }
		else { return c; }
	} else {
		if(b > c) { return b; }
		else { return c; }
	}
}

static inline double norm2(double a, double b, double c) {
	return a*a+b*b+c*c;
}

void quicksort(double *a, int *b, int left, int right) {
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
	quicksort(a,b,left,j-1);
	quicksort(a,b,j+1,right);
}

int * argsort(double *a, int size) {

	//copy a to keep a unchanged
	double *acpy = (double*) malloc(sizeof(double)*size);
	memcpy(acpy, a, sizeof(double)*size);

	int *ind = (int*) malloc(sizeof(int)*size);

	int i;
	//initialize array of indices to a
	for(i=0; i<size; i++) {
		ind[i] = i;
	}

	//sort a while performing duplicate operation on ind
	//ind becomes an array of indices to a which sorts a
	quicksort(acpy,ind, 0, size-1);
	free(acpy);
	return ind;
}

enum dim { X=0, Y=1, Z=2 };

void partition(double * vals, int * args, double key, int left, int right, 
				int median, int med_index) {
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

node_t * build(double *x, double *y, double *z,
					int *x_arg, int *y_arg, int *z_arg,
					int left, int right, enum dim d) {


	d=d%3;

	int med, med_arg;
	node_t *parent = (node_t *) malloc( sizeof(node_t) );
	//Median index of the sub-array. Rounds down for even sized lists
	med = (left+right)/2;

	//Find index of the median in appropriate position list
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

	//this node is the median
	parent->x = x[med_arg];
	parent->y = y[med_arg];
	parent->z = z[med_arg];

	parent->lchild = NULL;
	parent->rchild = NULL;

	/* Base cases: subtree of size 0, 1 or 2 */

	/* notes: separate max and min functions do more comparisons than necessary */
	if(right == left) {

		return parent;

	} else if(right-left == 1) { //length 2, one child

		parent->rchild = build(x,y,z,x_arg,y_arg,z_arg,right,right,d);
		return parent;

	} else if (right-left == 2) { //length 3, two children

		parent->lchild = build(x,y,z,x_arg,y_arg,z_arg,left,left,d);
		parent->rchild = build(x,y,z,x_arg,y_arg,z_arg,right,right,d);
		return parent;

	}


	//partition index array of other dims w.r.t. current dim
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


void destroy(node_t *p) {
	if(p->lchild != NULL) destroy(p->lchild);
	if(p->rchild != NULL) destroy(p->rchild);
	free(p);
}

/*  Query how many points in the tree with head p lie within radius r of point
    (x, y, z). Recursive. */
int radius(node_t *p, enum dim d, double x, double y, double z, double r) {
	d=d%3;
	int i;
	double rsq, dx, dy, dz;
	double dxmin, dxmax, dymin, dymax, dzmin, dzmax;
	double pos_upper, pos_lower, point;
	rsq = r*r;
	dx = p->x - x;
	dy = p->y - y;
	dz = p->z - z;

	double n = norm2(dx,dy,dz);

	if(p->lchild == NULL && p->rchild == NULL) {
		return (n < rsq);
	} else if (p->lchild == NULL) {
		return (n < rsq) + radius(p->rchild,d+1,x,y,z,r);
	} else {
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

kdtree_t tree_construct(int size, double x[], double y[], double z[]) {

	kdtree_t tree;
	tree.size = size;
	tree.x = x;
	tree.y = y;
	tree.z = z;

	// Argsort the inputs
	int *x_arg, *y_arg, *z_arg;
	x_arg = argsort(x, size);
	y_arg = argsort(y, size);
	z_arg = argsort(z, size);

	tree.root = build(x,y,z,x_arg,y_arg,z_arg,0,size-1,X);

	return tree;
}

long long two_point_correlation(kdtree_t tree, double x[], double y[],
    						double z[], int n, double r, MPI_Comm comm) {


    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);


	long long result = 0;

    int nlocal = n / mpi_size + (n % mpi_size > mpi_rank);

    int idx_min = 0;
    if (mpi_rank > 0)
        idx_min = mpi_rank * ( n / mpi_size ) + min( mpi_rank, n % mpi_size );

    int idx_max_plus_one = idx_min + nlocal;


    double t1, t2;
    int i;
    
    t1 = MPI_Wtime();
    for(i=idx_min; i<idx_max_plus_one; i++) {
		result += radius(tree.root, 0, x[i], y[i], z[i], r);
	}
    t2 = MPI_Wtime();

    printf("Time: %f sec\n", t2 - t1);

    //printf("before: %lld\n", result);
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INTEGER, MPI_SUM, comm);
    //printf("after: %lld\n", result);

	return result;
}
