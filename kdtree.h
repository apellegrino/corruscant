#define FLOAT double

#define KDTREE_H

#ifndef MPI_H
#define MPI_H
#include "mpi.h"
#endif

enum dim { X=0, Y=1, Z=2 };

typedef struct kdtree {
    struct node* root;
    int size;
} kdtree_t;

typedef struct node {
    double x, y, z;
    struct node * lchild;
    struct node * rchild;
} node_t;

kdtree_t tree_construct(int, FLOAT [], FLOAT [], FLOAT []);

long long two_point_correlation(kdtree_t tree, FLOAT [], FLOAT [],
                                    FLOAT [], int, FLOAT, int, MPI_Comm);
