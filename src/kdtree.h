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
#define ID_MASK 0b1111111100

enum dim { X=0, Y=1, Z=2 };


typedef struct node {
    FLOAT x, y, z;
    unsigned short flags;
} node_t;

typedef struct array3d {
    double *x; 
    double *y; 
    double *z;
    int size;
} array3d_t;

typedef struct argarray3d {
    int *x;
    int *y;
    int *z;
    int size;
} argarray3d_t;

typedef struct kdtree {
    node_t * node_data;
    int size;
    int memsize;
    array3d_t data;
    argarray3d_t arg_data;
    int * field_data;
} kdtree_t;

int left_child(int);
int right_child(int);
enum dim next_dim(enum dim);

array3d_t form_array(double *, double *, double *, int);
kdtree_t tree_construct(array3d_t, int *);

long long pair_count_jackknife(kdtree_t, array3d_t, double, int, int);
long long pair_count_noerr(kdtree_t, array3d_t, double, int);
