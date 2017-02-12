//#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include "kdstruct.h"

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

static double max3(double a, double b, double c) {
	if(a < b) {
		if(a < c) { return a; }
		else { return c; }
	} else {
		if(b < c) { return b; }
		else { return c; }
	}
}

static double min3(double a, double b, double c) {
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

/*
struct node_sph {
	double theta, phi;
	struct node * lchild;
	struct node * rchild;
}

int data_data(struct kdtree t) {
	int i;
	for(i=0; i<size; i++) {
		
	}	
}
*/
void partition(double * vals, int * args, double key, int left, int right, 
				int median, int med_index) {
	int i, j, k, arg, size, split;
	split = med_index - left;
	j = 0;
	k = split+1;
	size = right-left+1;
	int *args_new = (int*) malloc(sizeof(int)*size);

	//for(i=0; i<size; i++) args_new[i] = -1;
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
	//free(args_new);
}

struct node * build(double *x, double *y, double *z,
					int *x_arg, int *y_arg, int *z_arg,
					int left, int right, enum dim d) {
	d=d%3;

	int med, med_arg;
	struct node *parent = (struct node *) malloc( sizeof(struct node) );
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
		parent->xmin = parent->x;
		parent->xmax = parent->x;
		parent->ymin = parent->y;
		parent->ymax = parent->y;
		parent->zmin = parent->z;
		parent->zmax = parent->z;
		parent->size = 1;
		return parent;
	} else if(right-left == 1) { //length 2, one child
		parent->rchild = build(x,y,z,x_arg,y_arg,z_arg,right,right,d);

		parent->xmin = min(parent->x,parent->rchild->x);
		parent->xmax = max(parent->x,parent->rchild->x);
		parent->ymin = min(parent->y,parent->rchild->y);
		parent->ymax = max(parent->y,parent->rchild->y);
		parent->zmin = min(parent->z,parent->rchild->z);
		parent->zmax = max(parent->z,parent->rchild->z);

		parent->size = 1 + parent->rchild->size;
	/*
		printf("%i %i (%f %f %f) -> (%f %f %f)\n", med_arg, d, parent->x,parent->y,parent->z,
							parent->rchild->x,parent->rchild->y,parent->rchild->z);
	*/
		return parent;

	} else if (right-left == 2) { //length 3, two children
		parent->lchild = build(x,y,z,x_arg,y_arg,z_arg,left,left,d);
		parent->rchild = build(x,y,z,x_arg,y_arg,z_arg,right,right,d);

        // compute minimum coords of parent and all children
		parent->xmin = min3(parent->x,parent->lchild->xmin,parent->rchild->xmin);
		parent->xmax = max3(parent->x,parent->lchild->xmax,parent->rchild->xmax);
		parent->ymin = min3(parent->y,parent->lchild->ymin,parent->rchild->ymin);
		parent->ymax = max3(parent->y,parent->lchild->ymax,parent->rchild->ymax);
		parent->zmin = min3(parent->z,parent->lchild->zmin,parent->rchild->zmin);
		parent->zmax = max3(parent->z,parent->lchild->zmax,parent->rchild->zmax);

		parent->size = 1 + parent->lchild->size + parent->rchild->size;
		/*
		printf("%i %i (%f %f %f) -> (%f %f %f) (%f %f %f)\n", med_arg, d, parent->x,parent->y,parent->z,
			parent->lchild->x,parent->lchild->y,parent->lchild->z,
			parent->rchild->x,parent->rchild->y,parent->rchild->z);
		*/
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

	parent->size = 1 + parent->lchild->size + parent->rchild->size;
	/*
	printf("%i %i (%f %f %f) -> (%f %f %f) (%f %f %f)\n", med_arg, d, parent->x,parent->y,parent->z,
			parent->lchild->x,parent->lchild->y,parent->lchild->z,
			parent->rchild->x,parent->rchild->y,parent->rchild->z);
	*/
			
	return parent;
}

void verify(struct node *root, enum dim d) {
	d=d%3;
	struct node *p = root;
    char *errmsg = "node %p should not be %s child of %p\n";

	if (root->lchild != NULL) {
		switch(d) {
			case X:
				if(root->x <= root->lchild->x) {
					fprintf(stderr,errmsg,
                    (void*)root->lchild,"left",(void*)root);
				}
				break;
			case Y:
				if(root->y <= root->lchild->y) {
					fprintf(stderr,errmsg,
                    (void*)root->lchild,"left",(void*)root);
				}
				break;
			case Z:
				if(root->z <= root->lchild->z) {
					fprintf(stderr,errmsg,
                    (void*)root->lchild,"left",(void*)root);
				}
				break;
		}
		verify(root->lchild,d+1);
	}
	if (root->rchild != NULL) {
		switch(d) {
			case X:
				if(root->x >= root->rchild->x) {
					fprintf(stderr,errmsg,
                    (void*)root->lchild,"right",(void*)root);
				}
				break;
			case Y:
				if(root->y >= root->rchild->y) {
					fprintf(stderr,errmsg,
                    (void*)root->lchild,"right",(void*)root);
				}
				break;
			case Z:
				if(root->z >= root->rchild->z) {
					fprintf(stderr,errmsg,
                    (void*)root->lchild,"right",(void*)root);
				}
				break;
		}
		verify(root->rchild,d+1);
	}
}

void destroy(struct node *p) {
	if(p->lchild != NULL) destroy(p->lchild);
	if(p->rchild != NULL) destroy(p->rchild);
	free(p);
}

int count(struct node *p) {
    if(p->rchild == NULL) return (p->lchild == NULL);
    return (count(p->lchild) + count(p->rchild));
}

int radius(struct node *p, enum dim d, double x, double y, double z, double r) {
	d=d%3;
	int i;
	double rsq, dx, dy, dz;
	double dxmin, dxmax, dymin, dymax, dzmin, dzmax;
	double pos_upper, pos_lower, point;
	rsq = r*r;
	dx = p->x - x;
	dy = p->y - y;
	dz = p->z - z;

	dxmin = p->xmin - x;
	dxmax = p->xmax - x;
	dymin = p->ymin - y;
	dymax = p->ymax - y;
	dzmin = p->zmin - z;
	dzmax = p->zmax - z;
	
	double n = norm2(dx,dy,dz);

	if(p->lchild == NULL && p->rchild == NULL) {
		return (n < rsq);
	} else if (p->lchild == NULL) {
		return (n < rsq) + radius(p->rchild,d+1,x,y,z,r);
	} else {
		int contained = ( norm2(dxmin,dymin,dzmin) < rsq &&
						norm2(dxmin,dymin,dzmax) < rsq &&
						norm2(dxmin,dymax,dzmin) < rsq &&
						norm2(dxmin,dymax,dzmax) < rsq &&
						norm2(dxmax,dymin,dzmin) < rsq &&
						norm2(dxmax,dymin,dzmax) < rsq &&
						norm2(dxmax,dymax,dzmin) < rsq &&
						norm2(dxmax,dymax,dzmax) < rsq );

		//printf("%d %d\n", p->size, contained);
		if(contained) {
			return p->size;
		}

		/*
		switch(d) {
			case X:
				if (x+r < p->x) {
					i = radius(p->lchild,d+1,x,y,z,r);
				} else if (x-r > p->x) {
					i = radius(p->rchild,d+1,x,y,z,r);
				} else {
					i= (int) (n < rsq) +
						radius(p->lchild,d+1,x,y,z,r) +
						radius(p->rchild,d+1,x,y,z,r);
				}
				break;
			case Y:
				if (y+r < p->y) {
					i=radius(p->lchild,d+1,x,y,z,r);
				} else if (y-r > p->y) {
					i=radius(p->rchild,d+1,x,y,z,r);
				} else {
					i= (int) (n < rsq) +
						radius(p->lchild,d+1,x,y,z,r) +
						radius(p->rchild,d+1,x,y,z,r);
				}
				break;
			case Z:
				if (z+r < p->z) {
					i=radius(p->lchild,d+1,x,y,z,r);
				} else if (z-r > p->z) {
					i=radius(p->rchild,d+1,x,y,z,r);
				} else {
					i= (int) (n < rsq) +
						radius(p->lchild,d+1,x,y,z,r) +
						radius(p->rchild,d+1,x,y,z,r);
				}
				break;
		}
		*/
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
