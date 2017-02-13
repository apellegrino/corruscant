//#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>

enum dim;

typedef struct kdtree kdtree_t;

typedef struct node node_t;

static inline double max(double, double);
static inline double min(double, double);
static inline double max3(double, double, double);
static inline double min3(double, double, double);
static inline double norm2(double, double, double);

static inline void swapDouble(double *, double *);
static inline void swapInt(int *, int *);
void quicksort(double *, int *, int, int);
int * argsort(double *, int);

void partition(double *, int *, double, int, int, int, int);
node_t * build(double *, double *, double *,
					int *, int *, int *,
					int, int, enum dim);
void verify(node_t *, enum dim);
void destroy(node_t *);
int radius(node_t *, enum dim, double, double, double, double);
kdtree_t tree_construct(int, double [], double [], double []);
long long two_point_correlation(kdtree_t tree, double [], double [],
									double [], int, double, MPI_Comm); 
