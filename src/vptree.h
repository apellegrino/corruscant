/*  Andrew Pellegrino, 2017
 */

#define FLOAT double

#define VPTREE_H

#define HAS_LCHILD (1<<0)
#define HAS_RCHILD (1<<1)

#define ID_MASK_MAXINT ((1<<9) - 1)
// bits 2 through 10 (indexing from 0)
#define ID_MASK (ID_MASK_MAXINT << 2)

enum dim { X=0, Y=1, Z=2 };

typedef struct vpoint {
    FLOAT lat, lon;
} vpoint_t;

typedef struct node {
    vpoint_t point;
    double r;
    unsigned short flags;
} node_t;

typedef struct field_counter {
    long long * array;
    int size;
} field_counter_t;

typedef struct vptree {
    node_t * node_data;
    vpoint_t * data;
    int size;
    int memsize;
} vptree_t;

int left_child(int);
int right_child(int);

vptree_t make_vp_tree(vpoint_t *, int);

long long pair_count(vptree_t, double * lat, double * lon, int, double, int);

/*
long long * pair_count_jackknife(kdtree_t, array3d_t, double, int);
long long * pair_count_ftf(kdtree_t tree, array3d_t data, double, int);
long long * pair_count_noerr(kdtree_t, array3d_t, double, int);
*/
