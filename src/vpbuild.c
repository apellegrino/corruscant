/*
 * Andrew Pellegrino, 2017
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vptree.h"

static inline void swapPoint(vpoint_t * a, vpoint_t * b)
{
    vpoint_t temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

inline int left_child(int p)
{
    return 2*p;
}

inline int right_child(int p)
{
    return 2*p+1;
}

void destroy(vptree_t t)
{
    free(t.node_data);
}

static inline void set_id(node_t * node, int id)
{
    node->flags |= ID_MASK & (id << 2);
}

/* Flag node as having a left child */
static inline void set_lchild(node_t * node)
{
    node->flags |= HAS_LCHILD;
}

/* Flag node as having a right child */
static inline void set_rchild(node_t * node)
{
    node->flags |= HAS_RCHILD;
}

static inline int has_lchild(node_t p)
{
    return p.flags & HAS_LCHILD;
}

static inline int has_rchild(node_t p)
{
    return p.flags & HAS_RCHILD;
}

double vincenty(vpoint_t a, vpoint_t b)
{
    double c, d, e;
    double dlon = a.lon - b.lon;
    c = cos(b.lat) * sin(dlon);
    d = cos(a.lat) * sin(b.lat) - sin(a.lat) * cos(b.lat) * cos(dlon);
    e = sin(a.lat) * sin(b.lat) + cos(a.lat) * cos(b.lat) * cos(dlon);
    return atan2( sqrt( c*c + d*d ), e );
}

double great_circle(vpoint_t a, vpoint_t b)
{
    return acos( sin(a.lat) * sin(b.lat)
                + cos(a.lat) * cos(b.lat) * cos(a.lon - b.lon) );
}

double distance(vpoint_t a, vpoint_t b)
{
    double c, d, e;
    double dlon = a.lon - b.lon;
    c = cos(b.lat) * sin(dlon);
    d = cos(a.lat) * sin(b.lat) - sin(a.lat) * cos(b.lat) * cos(dlon);
    e = sin(a.lat) * sin(b.lat) + cos(a.lat) * cos(b.lat) * cos(dlon);
    return atan2( sqrt( c*c + d*d ), e );

    /*
    double dx = cos(b.lat) * cos(b.lon) - cos(a.lat) * cos(a.lon);
    double dy = cos(b.lat) * sin(b.lon) - cos(a.lat) * sin(a.lon);
    double dz = sin(b.lat) - sin(a.lat);
    return 2 * asin( sqrt(dx*dx+dy*dy+dz*dz) / 2.0 );
    */
}


/*
 * Shuffle `list` so that the nth element is the same as if `list` were sorted.
 * Sorting is based on distance to `vantage`. Uses quickselect, which acts like
 * quicksort except function only recurses to the side that n is in
 */
void nth_element(vpoint_t * list, int left, int right, vpoint_t vantage, int n)
{
    if ( right <= left ) return;

    vpoint_t pivot = list[left];
    int i = left+1;
    int j = right;

    while(1) {
        while( distance(vantage,list[i]) <= distance(vantage,pivot) && i <= right ) i++;
        while( distance(vantage,list[j]) >  distance(vantage,pivot) ) j--;

        if ( i >= j ) break;
        swapPoint(list+i, list+j);
    }
    swapPoint(&list[left], &list[j]);

    if (j > n) {
        nth_element(list, left, j-1, vantage, n);
    } else if (j < n) {
        nth_element(list, j+1, right, vantage, n);
    }
    return;
}

static int pow2ceil(int x)
{
    int i = 1;
    while(x > 1) { x /= 2; i++; }
    while(i > 0) { x *= 2; i--; }
    return x;
}

void build(vptree_t * tree, int node_id, vpoint_t * data, int left, int right)
{
    node_t * node = &((tree->node_data)[node_id]);
    // select root -- 0th for now
    node->point = data[left];
    if (left == right) {
        node->r = 0.0;
        return;
    }

    if (right - left == 1) {
        set_rchild(node);
        node_t * rc = &((tree->node_data)[right_child(node_id)]);
        rc->point = data[right];
        rc->r = 0.0;
        node->r = distance(node->point, rc->point);
        return;
    }

    int median = ((left+1)+right+1)/2;
    nth_element(data, left+1, right, node->point, median);
    node->r = distance(node->point, data[median]);

    set_lchild(node);
    set_rchild(node);
    build(tree, left_child(node_id), data, left + 1, median - 1);
    build(tree, right_child(node_id), data, median, right);

    return;
}

/*
void verify(node_t * data, int index)
{
    
    node_t *parent, *lc, *rc;
    //char *errmsg = "node %p should not be %s child of %p\n";

    parent = data + index;

	if (parent->flags & HAS_LCHILD) {
        lc = data + left_child(index);
		verify(data,left_child(index));
	}
	if (parent->flags & HAS_RCHILD) {
        rc = data + right_child(index);
		verify(data,right_child(index));
	}

}
*/

vptree_t make_vp_tree(vpoint_t * data, int length)
{
    vptree_t tree;
    tree.size = length;
    tree.memsize = pow2ceil(length);

    // allocate all space for the tree in one malloc
    tree.node_data = (node_t *) calloc( tree.memsize, sizeof(node_t) );

    //printf("%d size %d memsize\n", tree.size, tree.memsize);
    build(&tree, 1, data, 0, length-1);
    //verify(tree.node_data, 1);

    return tree;
}

vptree_t make_vp_python(double * lat, double * lon, int length)
{
    vptree_t tree;
    tree.size = length;
    tree.memsize = pow2ceil(length);

    // allocate all space for the tree in one malloc
    tree.node_data = (node_t *) calloc( tree.memsize, sizeof(node_t) );

    //printf("%d size %d memsize\n", tree.size, tree.memsize);
    vpoint_t * data = (vpoint_t *) malloc( length * sizeof(vpoint_t) );

    int i;
    for(i=0; i<length; i++) {
        data[i].lat = lat[i];
        data[i].lon = lon[i];
    }
    build(&tree, 1, data, 0, length-1);
    free(data);
    return tree;
}

long long radq(vptree_t tree, int pi, vpoint_t p, double radius)
{
    long long count = 0;
    node_t * node = tree.node_data + pi;

    if(distance(p, node->point) < radius) {
        count++;
    }

    if(!has_rchild(*node)) return count;
    if(!has_lchild(*node)) {
        count += radq(tree, right_child(pi), p, radius);
        return count;
    }

    if(distance(p, node->point) + radius < node->r) {
        count += radq(tree, left_child(pi), p, radius);
    } else if(distance(p, node->point) - radius >= node->r) {
        count += radq(tree, right_child(pi), p, radius);
    } else {
        count += radq(tree, left_child(pi), p, radius);
        count += radq(tree, right_child(pi), p, radius);
    }
    return count;
}

long long pair_count(vptree_t tree, double * lat, double * lon,
                    int length, double radius, int num_threads)
{
    long long count = 0;
    int i;
    vpoint_t q;
    for(i=0; i<length; i++) {
        q.lat = lat[i];
        q.lon = lon[i];
        count += radq(tree, 1, q, radius);
    }

    return count;
}
