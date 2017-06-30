#include "kdtree.h"
#include <stdio.h>

int nodesize(void) {
    return sizeof(node_t);
}

static inline unsigned long long convert(double x)
{
    /* could do return *((long long *) &x) here, but the following does not
     * raise a strict-aliasing warning from gcc */ 
    union u {
        double d;
        long long l;
    };
    return ((union u *) &x)->d;
}

unsigned long long node_hash(node_t * n)
{
    return convert(n->x) + convert(n->y) + convert(n->z);
}

unsigned long long _hash(node_t * data, int i)
{
    node_t n = *(data+i);
    unsigned long long sum = node_hash(data+i);

    if(n.has_lchild) {
        if(n.has_rchild) {
            sum += _hash(data,right_child(i));
        }
        sum += 2 * _hash(data,left_child(i));
    }

    return sum;
}

unsigned long long hash(kdtree_t t)
{
    return _hash(t.node_data,1);
}

static void _verify(node_t * data, int index, enum dim d) {
	d=d%3;

    node_t *parent, *lc, *rc;
    char *errmsg = "node %p should not be %s child of %p\n";
    int comp = 1;

    parent = data + index;

	if (parent->has_lchild) {
        lc = data + left_child(index);
		switch(d) {
        case X:
            comp = parent->x < lc->x;
            break;
        case Y:
            comp = parent->y < lc->y;
            break;
        case Z:
            comp = parent->z < lc->z;
            break;
		}
        if(comp) {
            fprintf(stderr,errmsg,lc,"left",parent);
        }
		_verify(data,left_child(index),d+1);
	}
	if (parent->has_rchild) {
        rc = data + right_child(index);
		switch(d) {
        case X:
            comp = parent->x > rc->x;
            break;
        case Y:
            comp = parent->y > rc->y;
            break;
        case Z:
            comp = parent->z > rc->z;
            break;
		}
        if(comp) {
            fprintf(stderr,errmsg,rc,"right",parent);
        }
		_verify(data,right_child(index),d+1);
	}
}

void verify_tree(kdtree_t t, enum dim d)
{
    _verify(t.node_data,1,d);
}

static int _count(kdtree_t t, int index)
{
    node_t * p = t.node_data + index;
    if(!(p->has_rchild)) return (p->has_lchild) + 1;
    return (_count(t, left_child(index)) + _count(t, right_child(index)) + 1);
}

int count_tree(kdtree_t t)
{
    return _count(t,1);
}
