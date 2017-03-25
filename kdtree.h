#define FLOAT double

#define KDTREE_H

#ifndef MPI_H
#define MPI_H
#include "mpi.h"
#endif

#define HAS_LCHILD 1
#define HAS_RCHILD 2

enum dim { X=0, Y=1, Z=2 };

typedef struct kdtree {
    struct node* root;
    int size;
} kdtree_t;

typedef struct node {
    double x, y, z;
    unsigned short flags;
} node_t;

node_t * tree_construct(int, FLOAT [], FLOAT [], FLOAT []);

long long two_point_correlation(node_t *, FLOAT [], FLOAT [],
                                    FLOAT [], int, FLOAT, int, MPI_Comm);
