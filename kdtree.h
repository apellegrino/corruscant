//#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include "kdstruct.h"

static inline void swapDouble(double * a, double * b);
static inline void swapInt(int * a, int * b);
void quicksort(double *a, int *b, int left, int right);
int * argsort(double *a, int size);
enum dim { X=0, Y=1, Z=2 };
void partition(double * vals, int * args, double key, int left, int right, 
				int median, int med_index);
struct node * build(double *x, double *y, double *z,
					int *x_arg, int *y_arg, int *z_arg,
					int left, int right, enum dim d);
void verify(struct node *root, enum dim d);
void destroy(struct node *p);
int radius(struct node *p, enum dim d, double x, double y, double z, double r);
