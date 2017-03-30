#define FLOAT double

#define KDTREE_H

#ifndef MPI_H
#define MPI_H
#include "mpi.h"
#endif

#define HAS_LCHILD 1
#define HAS_RCHILD 2

enum dim { X=0, Y=1, Z=2 };


typedef struct node {
    double x, y, z;
    unsigned short flags;
} node_t;

typedef struct kdtree {
    node_t * node_data;
    int size;
    int memsize;
    double *x_data, *y_data, *z_data;
} kdtree_t;

int left_child(int);
int right_child(int);
enum dim next_dim(enum dim);

kdtree_t tree_construct(int, FLOAT [], FLOAT [], FLOAT []);

long long two_point_correlation(kdtree_t, FLOAT [], FLOAT [],
                                    FLOAT [], int, FLOAT, int, MPI_Comm);
