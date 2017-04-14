/*  Andrew Pellegrino, 2017
 *
 *  K-d tree building algorithm adapted from Russell A. Brown, Building a
 *  Balanced k-d Tree in O(kn log n) Time, Journal of Computer Graphics
 *  Techniques (JCGT), vol. 4, no. 1, 50-68, 2015
 */

#define FLOAT double

#define KDTREE_H

#define HAS_LCHILD 1<<0
#define HAS_RCHILD 1<<1

enum dim { X=0, Y=1, Z=2 };


typedef struct node {
    FLOAT x, y, z;
    unsigned short flags;
} node_t;

typedef struct kdtree {
    node_t * node_data;
    int size;
    int memsize;
    double *x_data, *y_data, *z_data;
    int *x_arg, *y_arg, *z_arg;
} kdtree_t;

int left_child(int);
int right_child(int);
enum dim next_dim(enum dim);

kdtree_t tree_construct(int, double [], double [], double []);

long long pair_count(kdtree_t, double [], double [],
                     double [], int, double, int);
