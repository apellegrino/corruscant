/*  Andrew Pellegrino, 2017
 */

#define VPTREE_H

#define HAS_LCHILD (1<<0)
#define HAS_RCHILD (1<<1)

#define ID_MASK_MAXINT ((1<<9) - 1)
// bits 2 through 10 (indexing from 0)
#define ID_MASK (ID_MASK_MAXINT << 2)

typedef struct vpoint {
    double lat, lon;
    int field;
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
    int size;
    int memsize;
} vptree_t;

typedef struct vpdata {
    double * lat;
    double * lon;
    int * fields;
    int length;
    int num_fields;
} vpdata_t;

int left_child(int);
int right_child(int);

vptree_t make_vp_tree(double *, double *, int *, int);

long long * pair_count_noerr(vptree_t, double *, double *, int, double, int);
long long * pair_count_jackknife(vptree_t, double *, double *, int *, int, int, double, int);
long long * pair_count_ftf(vptree_t, double *, double *, int *, int, int, double, int);

/*
long long * pair_count_jackknife(kdtree_t, array3d_t, double, int);
long long * pair_count_ftf(kdtree_t tree, array3d_t data, double, int);
long long * pair_count_noerr(kdtree_t, array3d_t, double, int);
*/
